#!/usr/bin/env python3
"""
RNA Scoring Profile Plotter

Visualizes pseudo-energy scores as a function of distance for each base pair.

Usage:
  python plot_scores.py --input-dir training_output --output-dir plots
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


def load_score_table(filepath):
    """Load a score table CSV file."""
    if not os.path.exists(filepath):
        return None
    
    data = np.loadtxt(filepath, delimiter=',', skiprows=1)
    return {
        'distance_min': data[:, 0],
        'distance_max': data[:, 1],
        'distance_mid': data[:, 2],
        'score': data[:, 3]
    }


def plot_single_profile(pair, score_data, output_path):
    """Plot score profile for a single base pair."""
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.plot(score_data['distance_mid'], score_data['score'], 'o-', linewidth=2, markersize=4)
    ax.set_xlabel('Distance (Å)', fontsize=12)
    ax.set_ylabel('Pseudo-energy Score', fontsize=12)
    ax.set_title(f'Scoring Profile: {pair}', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(y=0, color='gray', linestyle='--', linewidth=1, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def plot_combined_profiles(pairs_data, output_path):
    """Plot all base pair profiles in a single figure."""
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()
    
    for idx, (pair, score_data) in enumerate(pairs_data.items()):
        ax = axes[idx]
        ax.plot(score_data['distance_mid'], score_data['score'], 'o-', linewidth=2, markersize=3)
        ax.set_xlabel('Distance (Å)', fontsize=10)
        ax.set_ylabel('Score', fontsize=10)
        ax.set_title(pair, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3)
        ax.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()


def main():
    parser = argparse.ArgumentParser(
        description="Plot RNA scoring profiles from trained tables"
    )
    parser.add_argument(
        '--input-dir',
        type=str,
        default='training_output',
        help='Directory containing score tables from train.py'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='plots',
        help='Directory for output plot files'
    )
    parser.add_argument(
        '--combined',
        action='store_true',
        help='Generate a combined plot with all profiles'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Base pairs
    pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'UC', 'UG', 'UU']
    
    print(f"Plotting scoring profiles...")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    
    # Load all score tables
    pairs_data = {}
    for pair in pairs:
        score_file = os.path.join(args.input_dir, f'score_{pair}.csv')
        score_data = load_score_table(score_file)
        
        if score_data is None:
            print(f"  WARNING: Score table not found for {pair}, skipping...")
            continue
        
        pairs_data[pair] = score_data
        
        # Plot individual profile
        output_path = os.path.join(args.output_dir, f'{pair}_profile.png')
        plot_single_profile(pair, score_data, output_path)
        print(f"  ✓ Plotted {pair}")
    
    # Generate combined plot if requested
    if args.combined and pairs_data:
        combined_path = os.path.join(args.output_dir, 'all_profiles.png')
        plot_combined_profiles(pairs_data, combined_path)
        print(f"  ✓ Generated combined plot")
    
    print(f"\n✓ Plotting complete!")
    print(f"  Generated {len(pairs_data)} individual plots")
    if args.combined:
        print(f"  Combined plot: all_profiles.png")
    
    return 0


if __name__ == '__main__':
    exit(main())
