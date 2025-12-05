#!/usr/bin/env python3
"""
RNA Scoring Function Training Script

Computes observed and reference frequencies from distance histograms,
then calculates pseudo-energy scores using the log-ratio formula.

Usage:
  python train.py --input-dir dist_data --output-dir training_output
"""

import argparse
import os
import numpy as np
import json
from pathlib import Path
import scipy.stats


def load_histogram(filepath):
    """Load a histogram file and return counts as numpy array."""
    if not os.path.exists(filepath):
        return None
    return np.loadtxt(filepath, dtype=float)


def load_kde_raw(filepath):
    """Load raw distances from KDE file."""
    if not os.path.exists(filepath):
        return None
    return np.loadtxt(filepath, dtype=float)


def compute_frequency(counts, pseudocount=1e-6):
    """
    Compute frequency from counts with pseudocount smoothing.
    
    Args:
        counts: Array of counts per bin
        pseudocount: Small value to avoid division by zero
    
    Returns:
        Frequency array (normalized counts)
    """
    total = np.sum(counts)
    if total == 0:
        return np.zeros_like(counts, dtype=float)
    
    # Add pseudocount to avoid zero frequencies
    smoothed_counts = counts + pseudocount
    smoothed_total = total + pseudocount * len(counts)
    
    return smoothed_counts / smoothed_total


def compute_kde(counts, bin_mid, bandwidth=None):
    """
    Compute KDE from histogram counts and bin midpoints.
    Returns density evaluated at bin midpoints.
    """
    # Reconstruct sample points from histogram
    samples = np.repeat(bin_mid, counts.astype(int))
    if len(samples) == 0:
        return np.zeros_like(bin_mid)
    kde = scipy.stats.gaussian_kde(samples, bw_method=bandwidth)
    density = kde.evaluate(bin_mid)
    return density / np.sum(density)  # Normalize


def compute_kde_from_raw(raw_distances, bin_mid, bandwidth=None):
    """
    Compute KDE density from raw distance samples evaluated at bin midpoints.
    
    Args:
        raw_distances: Array of raw distance values
        bin_mid: Bin midpoints where density is evaluated
        bandwidth: Bandwidth method for KDE
    
    Returns:
        Normalized density array
    """
    if len(raw_distances) == 0:
        return np.zeros_like(bin_mid)
    try:
        kde = scipy.stats.gaussian_kde(raw_distances, bw_method=bandwidth)
        density = kde.evaluate(bin_mid)
        return density / np.sum(density)  # Normalize
    except:
        return np.zeros_like(bin_mid)


def compute_score(obs_freq, ref_freq, max_score=10.0):
    """
    Compute pseudo-energy score: -log(f_obs / f_ref)
    
    Args:
        obs_freq: Observed frequency
        ref_freq: Reference frequency
        max_score: Maximum allowed score value
    
    Returns:
        Score array, capped at max_score
    """
    # Avoid division by zero
    # safe_obs = np.maximum(obs_freq, 1e-10)
    # safe_ref = np.maximum(ref_freq, 1e-10)
    
    score = -np.log(obs_freq / ref_freq)
    print("Debug info:")
    # print(obs_freq)
    # print(ref_freq)
    # print(obs_freq / ref_freq)
    # print(score)
    
    # Cap at maximum score
    score = np.minimum(score, max_score)
    score = np.nan_to_num(score, nan=max_score)
    print(score)
    
    return score


def main():
    parser = argparse.ArgumentParser(
        description="Train RNA scoring function from distance histograms"
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='dist_data',
        help='Directory containing histogram files or KDE data from extract_distances.py'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='training_output',
        help='Directory for output score and frequency tables'
    )
    parser.add_argument(
        '--max-score',
        type=float,
        default=10.0,
        help='Maximum score value (default: 10.0)'
    )
    parser.add_argument(
        '--pseudocount',
        type=float,
        default=1e-6,
        help='Pseudocount for smoothing (default: 1e-6)'
    )
    parser.add_argument(
        '--cutoff',
        type=float,
        default=20.0,
        help='Maximum distance in Angstroms (default: 20.0)'
    )
    # TODO
    # parser.add_argument(
    #     '--min-distance',
    #     type=float,
    #     default=2.0,
    #     help='Minimum distance in Angstroms (default: 2.0)'
    # )
    
    parser.add_argument(
        '--bin-width',
        type=float,
        default=1.0,
        help='Bin width in Angstroms (default: 1.0)'
    )
    parser.add_argument(
        '--method',
        choices=['histogram', 'kde'],
        default='histogram',
        help='Training method: histogram (default) or kde'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Base pairs to process
    pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'CU', 'GU', 'UU']
    
    # Load reference data (XX - all pairs)
    if args.method == 'kde':
        ref_file = os.path.join(args.input_dir, 'XX_kde_raw.txt')
        if not os.path.exists(ref_file):
            print(f"ERROR: KDE reference file not found: {ref_file}")
            print("Run extract_distances.py with --method kde first to generate raw distance files.")
            return 1
        ref_raw = load_kde_raw(ref_file)
        # Generate distance bin information from raw data
        n_bins = int(args.cutoff / args.bin_width)
        bin_edges = np.linspace(0, args.cutoff, n_bins + 1)
        bin_min = bin_edges[:-1]
        bin_max = bin_edges[1:]
        bin_mid = (bin_min + bin_max) / 2
        # Compute reference frequency using KDE
        ref_freq = compute_kde_from_raw(ref_raw, bin_mid)
    else:
        ref_file = os.path.join(args.input_dir, 'XX_histogram.txt')
        if not os.path.exists(ref_file):
            print(f"ERROR: Reference file not found: {ref_file}")
            print("Run extract_distances.py with --method histogram first to generate histograms.")
            return 1
        ref_counts = load_histogram(ref_file)
        # Generate distance bin information
        n_bins = len(ref_counts)
        # add warning if distance < 2 Å detected
        bin_edges = np.linspace(2, args.cutoff, n_bins + 1)

        if np.any(bin_edges < 2):
            print(f"WARNING: Detected distances < 2 Å in reference histogram. Check input data.")

        bin_min = bin_edges[:-1]
        bin_max = bin_edges[1:]
        bin_mid = (bin_min + bin_max) / 2
        # Compute reference frequency from histogram
        ref_freq = compute_frequency(ref_counts, args.pseudocount)
    
    print(f"Training scoring function...")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    print(f"  Number of bins: {n_bins}")
    print(f"  Distance range: 0-{args.cutoff} Å")
    print(f"  Max score: {args.max_score}")
    
    # Save reference frequency table
    if args.method == 'kde':
        ref_table = np.column_stack([bin_min, bin_max, bin_mid, ref_freq])
        header = 'Distance_Min,Distance_Max,Distance_Mid,Frequency'
    else:
        ref_table = np.column_stack([bin_min, bin_max, bin_mid, ref_counts, ref_freq])
        header = 'Distance_Min,Distance_Max,Distance_Mid,Count,Frequency'
    
    ref_path = os.path.join(args.output_dir, 'freq_XX.csv')
    np.savetxt(
        ref_path,
        ref_table,
        delimiter=',',
        header=header,
        comments='',
        fmt='%.6f'
    )
    
    # Process each base pair
    processed_pairs = []
    for pair in pairs:
        if args.method == 'kde':
            data_file = os.path.join(args.input_dir, f'{pair}_kde_raw.txt')
            if not os.path.exists(data_file):
                print(f"  WARNING: KDE raw file not found for {pair}, skipping...")
                continue
            obs_raw = load_kde_raw(data_file)
            if obs_raw is None or len(obs_raw) == 0:
                print(f"  WARNING: {pair} has no data, skipping...")
                continue
            obs_freq = compute_kde_from_raw(obs_raw, bin_mid)
        else:
            hist_file = os.path.join(args.input_dir, f'{pair}_histogram.txt')
            if not os.path.exists(hist_file):
                print(f"  WARNING: Histogram not found for {pair}, skipping...\"")
                continue
            obs_counts = load_histogram(hist_file)
            if obs_counts is None or len(obs_counts) != n_bins:
                print(f"  WARNING: {pair} has different bin count, skipping...")
                continue
            obs_freq = compute_frequency(obs_counts, args.pseudocount)
        
        # Compute score
        score = compute_score(obs_freq, ref_freq, args.max_score)
        
        # Save frequency table
        freq_table = np.column_stack([bin_min, bin_max, bin_mid, obs_freq])
        freq_path = os.path.join(args.output_dir, f'freq_{pair}.csv')
        np.savetxt(
            freq_path,
            freq_table,
            delimiter=',',
            header='Distance_Min,Distance_Max,Distance_Mid,Frequency',
            comments='',
            fmt='%.6f'
        )
        
        # Save score table
        score_table = np.column_stack([bin_min, bin_max, bin_mid, score])
        score_path = os.path.join(args.output_dir, f'score_{pair}.csv')
        np.savetxt(
            score_path,
            score_table,
            delimiter=',',
            header='Distance_Min,Distance_Max,Distance_Mid,Score',
            comments='',
            fmt='%.6f'
        )
        
        processed_pairs.append(pair)
        print(f"  ✓ Processed {pair}")
    
    # Save metadata
    metadata = {
        'input_dir': args.input_dir,
        'output_dir': args.output_dir,
        'max_score': args.max_score,
        'pseudocount': args.pseudocount,
        'cutoff': args.cutoff,
        'bin_width': args.bin_width,
        'n_bins': int(n_bins),
        'pairs': processed_pairs
    }
    
    metadata_path = os.path.join(args.output_dir, 'metadata.json')
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\n✓ Training complete!")
    print(f"  Processed {len(processed_pairs)} base pairs")
    print(f"  Frequency tables: freq_*.csv")
    print(f"  Score tables: score_*.csv")
    print(f"  Metadata: metadata.json")
    
    return 0


if __name__ == '__main__':
    exit(main())
