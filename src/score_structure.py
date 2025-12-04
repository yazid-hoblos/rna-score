#!/usr/bin/env python3
"""
RNA Structure Scoring Script

Evaluates predicted RNA structures using trained scoring tables.
Computes total Gibbs free energy estimate by summing interpolated scores.

Usage:
  python score_structure.py --pdb structure.pdb --tables training_output
"""

import argparse
import os
import numpy as np
from core import FastParser, DistanceComputer


class ScoringTables:
    """Container for trained scoring tables with interpolation."""
    
    def __init__(self, tables_dir):
        """Load all scoring tables from directory."""
        self.pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'UC', 'UG', 'UU']
        self.tables = {}
        
        for pair in self.pairs:
            score_file = os.path.join(tables_dir, f'score_{pair}.csv')
            if not os.path.exists(score_file):
                continue
            
            data = np.loadtxt(score_file, delimiter=',', skiprows=1)
            self.tables[pair] = {
                'distance_min': data[:, 0],
                'distance_max': data[:, 1],
                'distance_mid': data[:, 2],
                'score': data[:, 3]
            }
    
    def get_score(self, pair, distance):
        """
        Get interpolated score for a base pair at a given distance.
        
        Args:
            pair: Base pair type (e.g., 'AA', 'AU')
            distance: Distance in Angstroms
        
        Returns:
            Interpolated score value
        """
        if pair not in self.tables:
            return 0.0
        
        table = self.tables[pair]
        distances = table['distance_mid']
        scores = table['score']
        
        # Linear interpolation
        score = np.interp(distance, distances, scores, left=scores[0], right=scores[-1])
        return float(score)
    
    def canonical_pair_key(self, res1, res2):
        """Convert residue pair to canonical key (alphabetically sorted)."""
        res1 = res1.strip().upper()
        res2 = res2.strip().upper()
        
        # Convert modified bases to standard
        base_map = {'A': 'A', 'U': 'U', 'C': 'C', 'G': 'G', 'T': 'U'}
        res1 = base_map.get(res1, res1)
        res2 = base_map.get(res2, res2)
        
        if res1 > res2:
            res1, res2 = res2, res1
        
        return f"{res1}{res2}"


def parse_structure(filepath, file_format='pdb'):
    """Parse structure file and extract distance data."""
    try:
        with open(filepath, 'r') as f:
            content = f.read().splitlines()
    except Exception as e:
        print(f"ERROR: Could not read file {filepath}: {e}")
        return None
    
    # Parse structure
    parser = FastParser(atom_mode="C3'")
    data = parser.parse(
        content,
        file_format=file_format,
        pdb_name=os.path.basename(filepath),
        chains_filter=None,
        use_all_models=False
    )
    
    return data


def compute_structure_score(data, scoring_tables, cutoff=20.0, seq_sep=4):
    """
    Compute total score for a structure.
    
    Args:
        data: Parsed structure data
        scoring_tables: ScoringTables object
        cutoff: Maximum distance to consider
        seq_sep: Minimum sequence separation
    
    Returns:
        Dictionary with total score and detailed interactions
    """
    # Compute distances
    computer = DistanceComputer(cutoff=cutoff, seq_separation=seq_sep)
    dists, detailed_df = computer.compute(data, dist_mode='intra', detailed_output=True)
    
    total_score = 0.0
    interaction_scores = []
    
    # Score each interaction
    if detailed_df is not None and len(detailed_df) > 0:
        for _, row in detailed_df.iterrows():
            res1 = row['Res1']
            res2 = row['Res2']
            distance = row['Distance']
            
            pair_key = scoring_tables.canonical_pair_key(res1, res2)
            score = scoring_tables.get_score(pair_key, distance)
            
            total_score += score
            interaction_scores.append({
                'res1': res1,
                'res2': res2,
                'chain1': row['Chain1'],
                'chain2': row['Chain2'],
                'seq1': row['ResID1'],
                'seq2': row['ResID2'],
                'distance': distance,
                'pair': pair_key,
                'score': score
            })
    
    return {
        'total_score': total_score,
        'n_interactions': len(interaction_scores),
        'interactions': interaction_scores
    }


def score_single_structure(filepath, scoring_tables, file_format, cutoff, seq_sep):
    """Score a single structure and return results."""
    data = parse_structure(filepath, file_format)
    if data is None:
        return None
    
    result = compute_structure_score(data, scoring_tables, cutoff, seq_sep)
    result['filename'] = os.path.basename(filepath)
    return result


def main():
    parser = argparse.ArgumentParser(
        description="Score RNA structure(s) using trained scoring function"
    )
    
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        '--pdb',
        type=str,
        help='Single PDB or mmCIF file to score'
    )
    group.add_argument(
        '--folder',
        type=str,
        help='Directory containing PDB or mmCIF files to score in batch'
    )
    group.add_argument(
        '--list',
        type=str,
        help='Text file with list of PDB/mmCIF file paths to score'
    )
    
    parser.add_argument(
        '--tables',
        type=str,
        default='training_output',
        help='Directory containing trained score tables'
    )
    parser.add_argument(
        '--format',
        choices=['pdb', 'mmcif'],
        default='pdb',
        help='Input file format (default: pdb)'
    )
    parser.add_argument(
        '--cutoff',
        type=float,
        default=20.0,
        help='Distance cutoff for scoring (default: 20.0 Å)'
    )
    parser.add_argument(
        '--seq-sep',
        type=int,
        default=4,
        help='Minimum sequence separation for interactions (default: 4)'
    )
    parser.add_argument(
        '--detailed',
        default=False,
        action='store_true',
        help='Print detailed per-interaction scores'
    )
    parser.add_argument(
        '--output',
        type=str,
        help='Save scores to CSV file (required for batch mode)'
    )
    args = parser.parse_args()
    
    # Load scoring tables
    print(f"Loading scoring tables from {args.tables}...")
    scoring_tables = ScoringTables(args.tables)
    
    if not scoring_tables.tables:
        print("ERROR: No scoring tables found!")
        return 1
    
    print(f"  Loaded {len(scoring_tables.tables)} base pair tables")
    
    # Collect files to score
    files_to_score = []
    
    if args.pdb:
        files_to_score.append(args.pdb)
    
    elif args.folder:
        import glob
        exts = ['pdb', 'ent', 'cif', 'mmcif']
        for ext in exts:
            pattern = os.path.join(args.folder, f"*.{ext}")
            files_to_score.extend(glob.glob(pattern))
        
        if not files_to_score:
            print(f"ERROR: No structure files found in {args.folder}")
            return 1
    
    elif args.list:
        with open(args.list, 'r') as f:
            for line in f:
                filepath = line.strip()
                if filepath and os.path.exists(filepath):
                    files_to_score.append(filepath)
        
        if not files_to_score:
            print(f"ERROR: No valid files found in list {args.list}")
            return 1
    
    # Batch mode requires output file
    if len(files_to_score) > 1 and not args.output:
        print("ERROR: Batch mode requires --output to save results")
        return 1
    
    # Score structures
    all_results = []
    
    for idx, filepath in enumerate(files_to_score, 1):
        print(f"\n[{idx}/{len(files_to_score)}] Scoring: {os.path.basename(filepath)}...")
        
        result = score_single_structure(filepath, scoring_tables, args.format, args.cutoff, args.seq_sep)
        
        if result is None:
            print(f"  WARNING: Failed to process {filepath}, skipping...")
            continue
        
        all_results.append(result)
        
        # Print individual results for single file mode
        if len(files_to_score) == 1:
            print(f"\n{'='*60}")
            print(f"SCORING RESULTS")
            print(f"{'='*60}")
            print(f"Structure: {result['filename']}")
            print(f"Number of interactions: {result['n_interactions']}")
            print(f"Total Gibbs free energy estimate: {result['total_score']:.3f}")
            print(f"{'='*60}")
            
            # Print detailed scores if requested
            if args.detailed and result['interactions']:
                print(f"\nDetailed interaction scores:")
                print(f"{'Res1':<6} {'Res2':<6} {'Chain1':<7} {'Chain2':<7} {'Seq1':<6} {'Seq2':<6} {'Dist(Å)':<8} {'Pair':<6} {'Score':<8}")
                print("-" * 70)
                for inter in result['interactions'][:20]:  # Show first 20
                    print(f"{inter['res1']:<6} {inter['res2']:<6} {inter['chain1']:<7} {inter['chain2']:<7} "
                          f"{inter['seq1']:<6} {inter['seq2']:<6} {inter['distance']:<8.2f} "
                          f"{inter['pair']:<6} {inter['score']:<8.3f}")
                
                if len(result['interactions']) > 20:
                    print(f"... ({len(result['interactions']) - 20} more interactions)")
        else:
            # Brief output for batch mode
            print(f"  Score: {result['total_score']:.3f} ({result['n_interactions']} interactions)")
    
    # Save results
    if args.output and all_results:
        import pandas as pd
        
        if len(files_to_score) == 1:
            # Single file: save detailed interactions
            df = pd.DataFrame(all_results[0]['interactions'])
            df.to_csv(args.output, index=False)
            print(f"\nDetailed scores saved to: {args.output}")
        else:
            # Batch mode: save summary
            summary_data = []
            for r in all_results:
                summary_data.append({
                    'Filename': r['filename'],
                    'Total_Score': r['total_score'],
                    'N_Interactions': r['n_interactions']
                })
            df = pd.DataFrame(summary_data)
            df.to_csv(args.output, index=False)
            print(f"\n✓ Scored {len(all_results)}/{len(files_to_score)} structures")
            print(f"Summary saved to: {args.output}")
    print(f"Number of interactions: {result['n_interactions']}")
    print(f"Total Gibbs free energy estimate: {result['total_score']:.3f}")
    print(f"{'='*60}")
    
    # Print detailed scores if requested
    if args.detailed and result['interactions']:
        print(f"\nDetailed interaction scores:")
        print(f"{'Res1':<6} {'Res2':<6} {'Chain1':<7} {'Chain2':<7} {'Seq1':<6} {'Seq2':<6} {'Dist(Å)':<8} {'Pair':<6} {'Score':<8}")
        print("-" * 70)
        for inter in result['interactions'][:20]:  # Show first 20
            print(f"{inter['res1']:<6} {inter['res2']:<6} {inter['chain1']:<7} {inter['chain2']:<7} "
                  f"{inter['seq1']:<6} {inter['seq2']:<6} {inter['distance']:<8.2f} "
                  f"{inter['pair']:<6} {inter['score']:<8.3f}")
        
        if len(result['interactions']) > 20:
            print(f"... ({len(result['interactions']) - 20} more interactions)")
    
    return 0


if __name__ == '__main__':
    exit(main())
