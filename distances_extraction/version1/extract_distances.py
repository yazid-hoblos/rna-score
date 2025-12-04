"""
RNA DISTANCE EXTRACTOR
======================

This script calculates atomic distances from a dataset of PDB structures.
It utilizes the 'core' library to process files and generates raw data 
files suitable for statistical analysis (Histograms or KDE).

Main Capabilities:
------------------
1. Process multiple PDBs from a list file.
2. Filter for specific chains, models, or atom types.
3. Calculate Intrachain or Interchain distances.
4. Output data as binned counts (Histograms) or raw distance lists (KDE).

Arguments:
----------
Input:
  --pdb <ID>            Single PDB ID (e.g., 1EHZ) or local filename.
  --list <FILE>         Text file containing a list of PDB IDs and optional chain selections.

Filters:
  --chains <LIST>       Specific chains to process (only for single --pdb input).
  --format <FMT>        Input format: 'pdb' (default) or 'mmcif'.
  --all-models          Process all NMR models (Default: Model 1 only).

Parameters:
  --atom-mode <ARGS>    Which atoms to use?
                        - "C3'" (default)
                        - "centroid" (Geometric center)
                        - "all" (All atoms)
                        - <NAME> (e.g., "P")
                        - <LIST> (e.g., "P C3'")
  --dist-mode <MODE>    "intra" (folding) or "inter" (complexes).
  --cutoff <FLOAT>      Max distance in Angstroms (default: 20.0).
  --seq-sep <INT>       Min sequence separation for intrachain (default: 4).
  --bin-width <FLOAT>   Width of histogram bins (default: 1.0).

Output:
  --method <METHOD>     "histogram" (counts) or "kde" (raw values).
  --out-dir <DIR>       Output folder (default: "dist_data").
  --save-detailed       Save a CSV with detailed atom-pair info.

Usage Examples:
---------------
  python extract_distances.py --list training_set.txt --out-dir results
  python extract_distances.py --pdb 1EHZ --chains A --atom-mode all
"""

import argparse
import os
import numpy as np
import pandas as pd
from core import FastParser, OnlineFetcher, DistanceComputer

def parse_input_file(filename):
    """
    Reads the input file containing PDB IDs and chain selections.
    """
    targets = []
    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if not parts: continue
            
            pdb_id = parts[0]
            chains = parts[1:] if len(parts) > 1 else None
            is_file = os.path.exists(pdb_id)
            
            targets.append({
                'id': pdb_id,
                'chains': chains,
                'is_local': is_file
            })
    return targets

def process_dataset(targets, atom_mode="C3'", dist_mode='intra', cutoff=20.0, 
                   seq_sep=4, file_format='pdb', all_models=False, detailed=False):
    """
    Iterates through the dataset and processes each structure.
    Returns aggregated results.
    """
    parser_engine = FastParser(atom_mode=atom_mode)
    computer = DistanceComputer(cutoff=cutoff, seq_separation=seq_sep)
    
    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']
    final_results = {k: [] for k in keys}
    all_detailed_dfs = []

    # Handle atom mode display string
    mode_str = atom_mode if isinstance(atom_mode, str) else "+".join(atom_mode)
    print(f"Starting Processing (Mode: {mode_str}, Interaction: {dist_mode})")

    count_total = 0

    for t in targets:
        ident = t['id']
        chains = t['chains']
        is_local = t['is_local']
        
        print(f"Processing {ident}...", end=" ", flush=True)
        
        content = None
        current_fmt = file_format
        
        if is_local:
            if ident.endswith('.cif'): current_fmt = 'mmcif'
            try:
                with open(ident, 'r') as f:
                    content = f.read().splitlines()
            except Exception as e:
                print(f"Failed to read file: {e}")
                continue
        else:
            content = OnlineFetcher.get_content(ident, file_format=file_format)
            
        if not content:
            print("Skipped (Empty/Download Failed).")
            continue

        # Parse
        data = parser_engine.parse(content, file_format=current_fmt, pdb_name=ident, chains_filter=chains, use_all_models=all_models)
        
        if data is None:
            print("No valid RNA atoms found.")
            continue
            
        # Compute
        dists, detailed_df = computer.compute(data, dist_mode=dist_mode, detailed_output=detailed)
        
        count = 0
        for k, v in dists.items():
            if k in final_results:
                final_results[k].extend(v)
                count += len(v)
        
        if detailed_df is not None:
            all_detailed_dfs.append(detailed_df)
            
        print(f"Done. ({count} interactions)")
        count_total += count

    return final_results, all_detailed_dfs

def main():
    parser = argparse.ArgumentParser(description="RNA Distance Extractor")
    
    # Configuration
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--pdb', type=str, help="Single PDB ID or filename.")
    group.add_argument('--list', type=str, help="Input dataset file.")
    
    parser.add_argument('--chains', nargs='+', help="Specific chains to process.")
    parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help="File format.")
    parser.add_argument('--all-models', action='store_true', help="Process all NMR models.")
    
    parser.add_argument('--atom-mode', nargs='+', default=["C3'"], 
                        help="Atomic representation (e.g. 'C3'', 'all', 'P', or list).")
    
    parser.add_argument('--dist-mode', choices=['intra', 'inter'], default='intra', help="Interaction mode.")
    parser.add_argument('--cutoff', type=float, default=20.0, help="Max distance (Angstroms).")
    parser.add_argument('--seq-sep', type=int, default=4, help="Min separation for intrachain (default 4).")
    parser.add_argument('--bin-width', type=float, default=1.0, help="Histogram bin width.")
    
    parser.add_argument('--method', choices=['histogram', 'kde'], default='histogram', help="Output format.")
    parser.add_argument('--out-dir', type=str, default='dist_data', help="Output directory.")
    parser.add_argument('--save-detailed', action='store_true', help="Save detailed CSV.")
    
    args = parser.parse_args()

    # Normalize atom_mode
    mode_arg = args.atom_mode
    if len(mode_arg) == 1:
        mode_arg = mode_arg[0]

    # Run
    targets = []
    if args.pdb:
        is_local = os.path.exists(args.pdb)
        targets.append({'id': args.pdb, 'chains': args.chains, 'is_local': is_local})
    elif args.list:
        targets = parse_input_file(args.list)

    final_results, all_detailed_dfs = process_dataset(
        targets=targets,
        atom_mode=mode_arg,
        dist_mode=args.dist_mode,
        cutoff=args.cutoff,
        seq_sep=args.seq_sep,
        file_format=args.format,
        all_models=args.all_models,
        detailed=args.save_detailed
    )

    # Save
    os.makedirs(args.out_dir, exist_ok=True)
    print(f"\nWriting results to '{args.out_dir}/' (Format: {args.method.upper()})...")
    
    keys = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'CU', 'GG', 'GU', 'UU', 'XX']

    if args.method == 'histogram':
        bin_edges = np.arange(0, args.cutoff + args.bin_width, args.bin_width)
        for k in keys:
            vals = final_results[k]
            if not vals:
                hist = np.zeros(len(bin_edges) - 1, dtype=int)
            else:
                hist, _ = np.histogram(vals, bins=bin_edges)
            np.savetxt(f"{args.out_dir}/{k}_histogram.txt", hist, fmt='%d')
            
    elif args.method == 'kde':
        for k in keys:
            vals = final_results[k]
            if not vals:
                open(f"{args.out_dir}/{k}_kde_raw.txt", 'w').close()
            else:
                np.savetxt(f"{args.out_dir}/{k}_kde_raw.txt", vals, fmt='%.3f')

    if args.save_detailed and all_detailed_dfs:
        print("Writing detailed interaction log...")
        final_df = pd.concat(all_detailed_dfs, ignore_index=True)
        csv_path = os.path.join(args.out_dir, "detailed_interactions.csv")
        final_df.to_csv(csv_path, index=False)
        print(f"Log saved to {csv_path}")

    print("Pipeline completed successfully.")

if __name__ == "__main__":
    main()