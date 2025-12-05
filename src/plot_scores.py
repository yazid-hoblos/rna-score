#!/usr/bin/env python3
"""
RNA Scoring Profile Plotter

Visualizes pseudo-energy scores as a function of distance for each base pair.

Usage:
  python plot_scores.py --input-dir training_output --output-dir plots
"""

import argparse
import os
from pathlib import Path
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


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
    sns.set_theme(style="whitegrid")
    x = score_data['distance_mid']
    y = score_data['score']
    cmap = plt.get_cmap('plasma')
    norm = Normalize(vmin=np.min(y), vmax=np.max(y))

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    segment_values = (y[:-1] + y[1:]) / 2

    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=3)
    lc.set_array(segment_values)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.add_collection(lc)
    sc = ax.scatter(x, y, c=y, cmap=cmap, norm=norm, s=40, zorder=3)

    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y) - 0.1*np.ptp(y), np.max(y) + 0.1*np.ptp(y))
    ax.set_xlabel('Distance (Å)', fontsize=12)
    ax.set_ylabel('Pseudo-energy Score', fontsize=12)
    ax.set_title(f'Scoring Profile: {pair}', fontsize=14, fontweight='bold')
    ax.axhline(0, linestyle='--', color='gray', alpha=0.5, linewidth=1)
    ax.grid(alpha=0.3)

    cbar = fig.colorbar(sc, ax=ax)
    cbar.set_label('Score', fontsize=12)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    

def plot_single_profile_plotly(pair, score_data, output_path):
    x = score_data['distance_mid']
    y = score_data['score']
    y_norm = (y - np.min(y)) / (np.max(y) - np.min(y))
    segment_colors = sample_colorscale('Plasma', y_norm)
    traces = []
    for i in range(len(x)-1):
        trace = go.Scatter(
            x=x[i:i+2],
            y=y[i:i+2],
            mode='lines+markers',
            line=dict(color=segment_colors[i], width=3),
            marker=dict(color=segment_colors[i], size=6),
            showlegend=False,
            hoverinfo='x+y'
        )
        traces.append(trace)
    fig = go.Figure(data=traces)
    fig.update_layout(
        title=f'Scoring Profile: {pair}',
        xaxis_title='Distance (Å)',
        yaxis_title='Pseudo-energy Score',
        template='plotly_white',
        hovermode='closest'
    )
    fig.write_html(output_path)
    return fig


def plot_combined_profiles(pairs_data, output_path):
    """Plot all base pair profiles in a single figure."""
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()
    
    for idx, (pair, score_data) in enumerate(pairs_data.items()):
        ax = axes[idx]
        x = score_data['distance_mid']
        y = score_data['score']
        ax.plot(x, y, 'o-', linewidth=2, markersize=3)
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
    parser.add_argument(
        '--plotly',
        action='store_true',
        help='Use Plotly for interactive plots'
    )
    
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(os.path.join(args.output_dir, 'png'), exist_ok=True)
    
    pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'UC', 'UG', 'UU']
    
    print(f"Plotting scoring profiles...")
    print(f"  Input directory: {args.input_dir}")
    print(f"  Output directory: {args.output_dir}")
    
    pairs_data = {}
    for pair in pairs:
        score_file = os.path.join(args.input_dir, f'score_{pair}.csv')
        score_data = load_score_table(score_file)
        
        if score_data is None:
            print(f"  WARNING: Score table not found for {pair}, skipping...")
            continue
        
        pairs_data[pair] = score_data
        
        output_path = os.path.join(args.output_dir, 'png', f'{pair}_profile.png')
        plot_single_profile(pair, score_data, output_path)
        print(f"  ✓ Plotted {pair}")
    
    
    if args.combined and pairs_data:
        combined_path = os.path.join(args.output_dir,'png', 'all_profiles.png')
        plot_combined_profiles(pairs_data, combined_path)
        print(f"  ✓ Generated combined plot")

    if args.plotly:
        os.makedirs(os.path.join(args.output_dir, 'html'), exist_ok=True)
        for pair, score_data in pairs_data.items():
            output_path = os.path.join(args.output_dir, 'html', f'{pair}_profile.html')
            plot_single_profile_plotly(pair, score_data, output_path)
            print(f"  ✓ Plotted {pair} (Plotly)")
    
    print(f"\n✓ Plotting complete!")
    print(f"  Generated {len(pairs_data)} individual plots")
    if args.combined:
        print(f"  Combined plot: all_profiles.png")
    
    return 0


if __name__ == '__main__':
    exit(main())
