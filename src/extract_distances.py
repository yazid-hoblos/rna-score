"""
RNA Distance Extractor Script

This script processes RNA structures from PDB or mmCIF files (local or remote) to compute
distance distributions between selected atoms. It supports both intra-chain and inter-chain
interactions, various atomic selection modes, and parallel processing.

Parameters / Command-line options:

  -h, --help
      Show this help message and exit.

  --pdb PDB
      Single RNA PDB ID (remote) or filename (local). Default: None

  --list LIST
      Text file containing a list of RNA PDB IDs and optional chain IDs to process.
      Format: <PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]
      Example:
          1EHZ A
          1Y26 B C
          2OEU
      Default: None

  --folder FOLDER
      Directory containing local PDB or mmCIF files to process in batch.
      Default: None

  --chains CHAINS [CHAINS ...]
      Specific chains to process for a single PDB or files in a folder. Default: all chains

  --format {pdb,mmcif}
      Input file format. Default: pdb

  --all-models
      Process all NMR models instead of just the first model. Default: False

  --cores CORES
      Number of CPU cores to use for parallel processing. Default: None (uses all available cores)

  --atom-mode ATOM_MODE [ATOM_MODE ...]
      Atom selection mode:
        - Single atom name (e.g., "C3'")
        - "centroid" for nucleotide centroid
        - "all" for all atoms
        - list of atom names (e.g., "C3'" "P" "O")
      Default: ["C3'"]

  --dist-mode {intra,inter}
      Interaction mode:
        - "intra" computes distances within the same chain (sequence-separated by --seq-sep)
        - "inter" computes distances between different chains
      Default: intra

  --cutoff CUTOFF
      Maximum distance (Å) to consider an interaction. Default: 20.0

  --seq-sep SEQ_SEP
      Minimum sequence separation for intra-chain interactions. Ignored for inter-chain. Default: 4

  --bin-width BIN_WIDTH
      Width of histogram bins (applies if --method=histogram). Default: 1.0

  --method {histogram,kde}
      Output type:
        - histogram: counts of interactions in bins
        - kde: raw distances (can be used for kernel density estimation)
      Default: histogram

  --out-dir OUT_DIR
      Directory for output files. Default: dist_data

  --save-detailed
      Export detailed CSV log of interactions including coordinates, atom names, chain IDs, and distances. Default: False
"""

import argparse
import os
import numpy as np
import pandas as pd
import time
from concurrent.futures import ProcessPoolExecutor
from itertools import repeat
import glob
from scipy.spatial import cKDTree
from collections import defaultdict
from utils.structure_io import FastParser, OnlineFetcher

# Global configuration removed; config is passed explicitly per worker

def parse_input_file(filename):
    """
    Parse a list of PDB IDs and optional chains from a text file.
    Each line: <PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]
    Returns a list of dictionaries: {'id': str, 'chains': list or None, 'is_local': bool}
    """
    targets = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts:
                continue
            pdb_id = parts[0]
            chains = parts[1:] if len(parts) > 1 else None
            is_file = os.path.exists(pdb_id)
            targets.append({
                'id': pdb_id,
                'chains': chains,
                'is_local': is_file
            })
    return targets

def worker_process(target, cfg):
    """
    Worker function to fetch, parse, and compute distances for a single target.
    Returns (results_dict, detailed_df) or (None, None) if processing fails.
    """
    try:
        ident = target['id']
        file_id = os.path.basename(ident).split('.')[0] # Use filename as ID for local files

        # Fetch file content
        if target['is_local']:
            try:
                with open(ident, 'r') as f:
                    content = f.read().splitlines()
                # Infer format for local files if not explicitly set
                if ident.lower().endswith(('.cif', '.mmcif')):
                    file_format = 'mmcif'
                else:
                    file_format = 'pdb'

            except:
                return None, None
        else:
            content = OnlineFetcher.get_content(ident, file_format=cfg['fmt'])
            file_format = cfg['fmt']
            file_id = ident # Use PDB ID for remote files
        
        if not content:
            return None, None

        # Parse structure
        parser = FastParser(atom_mode=cfg['atom_mode'])
        data = parser.parse(content, file_format=file_format, pdb_name=file_id,
                            chains_filter=target['chains'], use_all_models=cfg['all_models'])
        if data is None:
            return None, None

        # Compute distances
        dists, detailed_df = compute_distances(
            data,
            cutoff=cfg['cutoff'],
            seq_separation=cfg['seq_sep'],
            dist_mode=cfg['dist_mode'],
            detailed_output=cfg['detailed']
        )
        
        return dists, detailed_df
    except Exception:
        return None, None


def compute_distances(data, cutoff=20.0, seq_separation=4, dist_mode='intra', detailed_output=False):
    """
    Computes pairwise distances between atoms in an RNA structure.

    Parameters:
    -----------
    data : dict
        Parsed RNA structure data from FastParser.
    cutoff : float
        Maximum distance between atoms to be considered a contact (in Ångströms).
    seq_separation : int
        Minimum sequence separation for intra-chain distance computations.
    dist_mode : str
        Either 'intra' (within the same chain) or 'inter' (between different chains).
    detailed_output : bool, optional
        Whether to return a detailed pandas DataFrame with all distances.

    Returns:
    --------
    tuple:
        results : dict
            Dictionary mapping residue pair types to lists of distances.
        detailed_df : pandas.DataFrame or None
            DataFrame with detailed distance information if requested.
    """
    if data is None:
        return {}, None
    
    coords = data['coords']
    chains = data['chain_ids']
    res_ids = data['res_ids']
    res_names = data['res_names']
    
    atom_names = data['atom_names']
    alt_locs = data['alt_locs']
    b_factors = data['b_factors']
    model_ids = data['model_ids']
    pdb_name = data['pdb_name']
    
    # Build a KDTree for fast neighbor search
    tree = cKDTree(coords, compact_nodes=True, balanced_tree=True)
    pairs = tree.query_pairs(cutoff)
    
    if not pairs:
        return {}, None
    
    pairs_np = np.array(list(pairs))
    idx_i = pairs_np[:, 0]
    idx_j = pairs_np[:, 1]
    
    same_model = (model_ids[idx_i] == model_ids[idx_j])
    same_chain = (chains[idx_i] == chains[idx_j])
    
    if dist_mode == 'intra':
        seq_diff = np.abs(res_ids[idx_i] - res_ids[idx_j])
        mask = same_model & same_chain & (seq_diff >= seq_separation)
    else:
        mask = same_model & (~same_chain)

    final_i = idx_i[mask]
    final_j = idx_j[mask]
    
    diffs = coords[final_i] - coords[final_j]
    dists = np.linalg.norm(diffs, axis=1)
    
    types_i = res_names[final_i]
    types_j = res_names[final_j]
    
    results = defaultdict(list)
    pair_keys = []

    # Organize distances by residue type
    for t1, t2, d in zip(types_i, types_j, dists):
        key = t1 + t2 if t1 < t2 else t2 + t1
        pair_keys.append(key)
        results[key].append(d)
        results['XX'].append(d)
        
    detailed_df = None
    if detailed_output and len(dists) > 0:
        detailed_df = pd.DataFrame({
            'PDB': pdb_name,
            'Model': model_ids[final_i],
            'Chain1': chains[final_i],
            'Res1': types_i,
            'ResID1': res_ids[final_i],
            'Atom1': atom_names[final_i],
            'AltLoc1': alt_locs[final_i],
            'B_Factor1': b_factors[final_i],
            'Chain2': chains[final_j],
            'Res2': types_j,
            'ResID2': res_ids[final_j],
            'Atom2': atom_names[final_j],
            'AltLoc2': alt_locs[final_j],
            'B_Factor2': b_factors[final_j],
            'Distance': dists,
            'Type': pair_keys
        })
        
    return results, detailed_df


def process_dataset_parallel(targets, config, max_workers=None):
    """
    Process a list of targets in parallel using ProcessPoolExecutor.
    Returns final_results (dict of distance lists) and all_detailed_dfs (list of pandas DataFrames)
    """
    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']
    final_results = {k: [] for k in keys}
    all_detailed_dfs = []

    atom_mode_print = config['atom_mode']
    if isinstance(atom_mode_print, list):
        atom_mode_print = "+".join(atom_mode_print)

    print(f"Starting Parallel Processing (Atom Mode: {atom_mode_print})...")
    start_time = time.time()
    
    # Determine max_workers safely
    if max_workers is None or max_workers <= 0:
        max_workers = os.cpu_count() or 1

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(worker_process, targets, repeat(config))
        for i, (dists, df) in enumerate(results):
            if dists is None:
                continue
            for k, v in dists.items():
                if k in final_results:
                    final_results[k].extend(v)
            if df is not None:
                all_detailed_dfs.append(df)
            if (i + 1) % 10 == 0:
                print(f"  ...processed {i+1}/{len(targets)} structures")

    elapsed = time.time() - start_time
    print(f"Batch processing complete in {elapsed:.2f} seconds.")
    return final_results, all_detailed_dfs

def main():
    parser = argparse.ArgumentParser(description="RNA Distance Extractor")
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb', type=str, help="Single RNA PDB ID or filename.")
    group.add_argument('--list', type=str, help="Text file containing a list of RNA PDB IDs and optional chain IDs.")
    group.add_argument('--folder', type=str, help="Directory containing local PDB or mmCIF files to process in batch.")
    
    parser.add_argument('--chains', nargs='+', help="Specific chains to process. Default: all chains")
    parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help="Input file format. Default: pdb")
    parser.add_argument('--all-models', action='store_true', help="Process all NMR models. Default: False")
    parser.add_argument('--cores', type=int, default=None, help="Number of CPU cores to use. Default: all available cores")
    parser.add_argument('--atom-mode', nargs='+', default=["C3'"], help="Atom selection mode (e.g., C3', centroid, all, or list). Default: C3'")
    parser.add_argument('--dist-mode', choices=['intra', 'inter'], default='intra', help="Interaction mode. Default: intra")
    parser.add_argument('--cutoff', type=float, default=20.0, help="Maximum distance (Å) to consider an interaction. Default: 20.0")
    parser.add_argument('--seq-sep', type=int, default=4, help="Minimum sequence separation for intra-chain interactions. Default: 4")
    parser.add_argument('--bin-width', type=float, default=1.0, help="Width of histogram bins (applies if method=histogram). Default: 1.0")
    parser.add_argument('--method', choices=['histogram', 'kde'], default='histogram', help="Output type: counts (histogram) or raw distances (kde). Default: histogram")
    parser.add_argument('--out-dir', type=str, default='dist_data', help="Directory for output files. Default: dist_data")
    parser.add_argument('--save-detailed', action='store_true', help="Export detailed CSV log of interactions. Default: False")
    
    args = parser.parse_args()
    
    # Normalize atom mode
    mode_arg = args.atom_mode
    if len(mode_arg) == 1:
        mode_arg = mode_arg[0]

    config = {
        'atom_mode': mode_arg,
        'dist_mode': args.dist_mode,
        'cutoff': args.cutoff,
        'seq_sep': args.seq_sep,
        'fmt': args.format,
        'all_models': args.all_models,
        'detailed': args.save_detailed
    }

    targets = []
    if args.pdb:
        is_local = os.path.exists(args.pdb)
        targets.append({'id': args.pdb, 'chains': args.chains, 'is_local': is_local})
        # Single PDB processing uses one core unless it's a local file.
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=1 if not is_local else args.cores)
    
    elif args.list:
        targets = parse_input_file(args.list)
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=args.cores)

    elif args.folder:
        print(f"Searching for files in: {args.folder}")
        targets = []
        
        # Define extensions to search based on format preference
        exts = []
        # Look for both common PDB and mmCIF extensions regardless of explicit format argument
        exts.extend(['pdb', 'ent', 'cif', 'mmcif'])

        # Use glob to find all matching files efficiently
        for ext in exts:
            pattern = os.path.join(args.folder, f"*.{ext}")
            for filepath in glob.glob(pattern):
                targets.append({
                    'id': filepath,
                    'chains': args.chains, # Apply common chain filter to all files
                    'is_local': True
                })
        
        if not targets:
            print(f"No files found in {args.folder} matching {exts} extensions.")
            return

        print(f"Found {len(targets)} files to process. Starting batch processing...")
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=args.cores)
        
    else:
        # Should not be reached due to argparse mutual exclusivity, but for safety
        return 

    os.makedirs(args.out_dir, exist_ok=True)
    print(f"\nWriting results to '{args.out_dir}/' (Method: {args.method.upper()})...")

    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']

    if args.method == 'histogram':
        bin_edges = np.arange(0, args.cutoff + args.bin_width, args.bin_width)
        for k in keys:
            vals = final_results[k]
            hist = np.histogram(vals, bins=bin_edges)[0] if vals else np.zeros(len(bin_edges)-1, dtype=int)
            np.savetxt(f"{args.out_dir}/{k}_histogram.txt", hist, fmt='%d')
    elif args.method == 'kde':
        for k in keys:
            vals = final_results[k]
            file_path = f"{args.out_dir}/{k}_kde_raw.txt"
            if vals:
                np.savetxt(file_path, vals, fmt='%.3f')
            else:
                open(file_path, 'w').close()

    if args.save_detailed and all_detailed_dfs:
        print("Writing detailed interaction log...")
        final_df = pd.concat(all_detailed_dfs, ignore_index=True)
        csv_path = os.path.join(args.out_dir, "detailed_interactions.csv")
        final_df.to_csv(csv_path, index=False)
        print(f"Log saved to {csv_path}")

    print("Pipeline completed successfully.")

if __name__ == "__main__":
    import multiprocessing
    multiprocessing.freeze_support()
    main()