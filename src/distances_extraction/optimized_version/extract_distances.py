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

  --chains CHAINS [CHAINS ...]
      Specific chains to process for a single PDB. Default: all chains

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
        - list of atom names
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
from core import FastParser, OnlineFetcher, DistanceComputer

# Global configuration shared with worker processes
_WORKER_CONFIG = {}

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

def worker_process(target):
    """
    Worker function to fetch, parse, and compute distances for a single target.
    Returns (results_dict, detailed_df) or (None, None) if processing fails.
    """
    try:
        cfg = _WORKER_CONFIG
        ident = target['id']

        # Fetch file content
        if target['is_local']:
            try:
                with open(ident, 'r') as f:
                    content = f.read().splitlines()
            except:
                return None, None
        else:
            content = OnlineFetcher.get_content(ident, file_format=cfg['fmt'])
        
        if not content:
            return None, None

        # Parse structure
        parser = FastParser(atom_mode=cfg['atom_mode'])
        data = parser.parse(content, file_format=cfg['fmt'], pdb_name=ident,
                            chains_filter=target['chains'], use_all_models=cfg['all_models'])
        if data is None:
            return None, None

        # Compute distances
        computer = DistanceComputer(cutoff=cfg['cutoff'], seq_separation=cfg['seq_sep'])
        dists, detailed_df = computer.compute(data, dist_mode=cfg['dist_mode'], detailed_output=cfg['detailed'])
        
        return dists, detailed_df
    except Exception:
        return None, None

def process_dataset_parallel(targets, config, max_workers=None):
    """
    Process a list of targets in parallel using ProcessPoolExecutor.
    Returns final_results (dict of distance lists) and all_detailed_dfs (list of pandas DataFrames)
    """
    global _WORKER_CONFIG
    _WORKER_CONFIG = config
    
    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']
    final_results = {k: [] for k in keys}
    all_detailed_dfs = []

    atom_mode_print = config['atom_mode']
    if isinstance(atom_mode_print, list):
        atom_mode_print = "+".join(atom_mode_print)

    print(f"Starting Parallel Processing (Atom Mode: {atom_mode_print})...")
    start_time = time.time()

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(worker_process, targets))
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
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=1)
    elif args.list:
        targets = parse_input_file(args.list)
        final_results, all_detailed_dfs = process_dataset_parallel(targets, config, max_workers=args.cores)

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
