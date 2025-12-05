#!/usr/bin/env python3

import argparse
import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from plotly.colors import sample_colorscale
import plotly.graph_objects as go

def _ensure_dir(path: str):
    """Create directory if it does not exist."""
    Path(path).mkdir(parents=True, exist_ok=True)


def load_score_table(filepath):
    """Load a score table CSV file (handles KDE or non-KDE)."""
    if not os.path.exists(filepath):
        return None

    data = np.loadtxt(filepath, delimiter=',', skiprows=1)

    if 'kde' in filepath:
        return {'distance': data[:, 0], 'score': data[:, 1]}

    return {
        'distance_min': data[:, 0],
        'distance_max': data[:, 1],
        'distance': data[:, 2],
        'score': data[:, 3]
    }


def _compute_segments(x, y):
    """Compute segment coordinates and mid-segment color values."""
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    segment_values = (y[:-1] + y[1:]) / 2
    return segments, segment_values


def _compute_normalized_values(y):
    """Normalize y for color mapping (Plotly)."""
    return (y - np.min(y)) / (np.max(y) - np.min(y))


def plot_single_profile(pair, score_data, output_path):
    sns.set_theme(style="whitegrid")

    x, y = score_data['distance'], score_data['score']
    segments, segment_values = _compute_segments(x, y)

    cmap = plt.get_cmap('plasma')
    norm = Normalize(vmin=np.min(y), vmax=np.max(y))

    fig, ax = plt.subplots(figsize=(8, 5))

    lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=3)
    lc.set_array(segment_values)
    ax.add_collection(lc)

    sc = ax.scatter(x, y, c=y, cmap=cmap, norm=norm, s=40, zorder=3)

    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(y) - 0.1*np.ptp(y), np.max(y) + 0.1*np.ptp(y))

    ax.axhline(0, linestyle="--", color="gray", alpha=0.5)
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Pseudo-energy Score")
    ax.set_title(f"Scoring Profile: {pair}", fontsize=14, fontweight="bold")

    fig.colorbar(sc, ax=ax, label="Score")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_combined_profiles(pairs_data, output_path, cols=5):
    sns.set_theme(style="whitegrid")

    n_pairs = len(pairs_data)
    rows = int(np.ceil(n_pairs / cols))

    fig, axes = plt.subplots(rows, cols, figsize=(6 * cols, 4 * rows), squeeze=False)
    axes = axes.flatten()

    for ax, (pair, score_data) in zip(axes, pairs_data.items()):
        x, y = score_data['distance'], score_data['score']
        segments, segment_values = _compute_segments(x, y)

        cmap = plt.get_cmap('plasma')
        norm = Normalize(vmin=np.min(y), vmax=np.max(y))

        lc = LineCollection(segments, cmap=cmap, norm=norm, linewidth=2)
        lc.set_array(segment_values)
        ax.add_collection(lc)

        ax.scatter(x, y, c=y, cmap=cmap, norm=norm, s=20)

        ax.axhline(0, linestyle="--", color="gray", alpha=0.5)
        ax.set_title(f"Pair: {pair}", fontsize=12, fontweight="bold")
        ax.set_xlabel("Distance (Å)")
        ax.set_ylabel("Pseudo-energy Score")

    for ax in axes[n_pairs:]:
        ax.axis("off")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()


def plot_single_profile_plotly(pair, score_data, output_path):
    x, y = score_data["distance"], score_data["score"]

    y_norm = _compute_normalized_values(y)
    segment_colors = sample_colorscale("Plasma", y_norm)

    fig = go.Figure()
    for i in range(len(x) - 1):
        fig.add_trace(go.Scatter(
            x=x[i:i+2],
            y=y[i:i+2],
            mode="lines+markers",
            line=dict(color=segment_colors[i], width=3),
            marker=dict(color=segment_colors[i], size=6),
            hoverinfo="x+y",
            showlegend=False
        ))

    fig.update_layout(
        title=f"Scoring Profile: {pair}",
        xaxis_title="Distance (Å)",
        yaxis_title="Pseudo-energy Score",
        template="plotly_white",
        hovermode="closest"
    )

    fig.write_html(output_path)


def plot_combined_profiles_plotly(pairs_data, output_path):
    fig = go.Figure()
    trace_map = {}

    for idx, (pair, score) in enumerate(pairs_data.items()):
        x, y = score['distance'], score['score']
        y_norm = _compute_normalized_values(y)
        colors = sample_colorscale("Plasma", y_norm)

        trace_map[pair] = []
        visible = (idx == 0)

        for i in range(len(x) - 1):
            fig.add_trace(go.Scatter(
                x=x[i:i+2],
                y=y[i:i+2],
                mode="lines+markers",
                name=pair if i == 0 else None,
                line=dict(color=colors[i], width=3),
                marker=dict(color=colors[i], size=6),
                visible=visible
            ))
            trace_map[pair].append(len(fig.data) - 1)

    buttons = []
    all_visible = [True] * len(fig.data)

    for pair, indices in trace_map.items():
        visibility = [False] * len(fig.data)
        for i in indices:
            visibility[i] = True

        buttons.append({
            "label": pair,
            "method": "update",
            "args": [{"visible": visibility}, {"title": f"Scoring Profile: {pair}"}]
        })

    buttons.append({"label": "All", "method": "update", "args": [{"visible": all_visible},
                                                                 {"title": "Combined Scoring Profiles"}]})

    fig.update_layout(
        title=f"Scoring Profile: {list(pairs_data.keys())[0]}",
        xaxis_title="Distance (Å)",
        yaxis_title="Pseudo-energy Score",
        template="plotly_white",
        updatemenus=[{"buttons": buttons, "x": 1.15, "y": 0.95}],
        showlegend=False
    )

    fig.write_html(output_path)


def run_workflow(args):
    """Main workflow for generating RNA scoring profile plots."""

    _ensure_dir(args.output_dir)
    _ensure_dir(os.path.join(args.output_dir, "png"))

    pairs = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'CU', 'GU', 'UU']
    pairs_data = {}

    print(f"Plotting scoring profiles from {args.input_dir} …")

    for pair in pairs:
        path = os.path.join(args.input_dir, f"score_{pair}.csv")
        data = load_score_table(path)

        png_dir = os.path.join(args.output_dir, "png")
        _ensure_dir(png_dir)
        html_dir = os.path.join(args.output_dir, "html")
        _ensure_dir(html_dir)

        if data is None:
            print(f" [!WARNING] missing: {pair}")
            continue

        pairs_data[pair] = data
        out_path = os.path.join(png_dir, f"{pair}_profile.png")
        plot_single_profile(pair, data, out_path)
        print(f"  ✓ {pair}")

    for pair, data in pairs_data.items():
        out_html = os.path.join(html_dir, f"{pair}_profile.html")
        plot_single_profile_plotly(pair, data, out_html)

    if args.combined and pairs_data:
        combined_png = os.path.join(args.output_dir, "png", "all_profiles.png")
        plot_combined_profiles(pairs_data, combined_png)
        print("  ✓ Combined plot saved.")

        if args.combined:
            combined_html = os.path.join(html_dir, "all_profiles.html")
            plot_combined_profiles_plotly(pairs_data, combined_html)

        print("  ✓ Plotly HTML plots generated.")

    print("\n✓ Complete, plots saved in:", args.output_dir)
    print('\tplotly interactive plots in "' + os.path.join(args.output_dir, "html") + '"')
    print('\tstatic png plots in "' + os.path.join(args.output_dir, "png") + '"')
    return 0


def main():
    parser = argparse.ArgumentParser(description="Plot RNA scoring profiles from trained tables")
    parser.add_argument("--input-dir", type=str, default="training_output")
    parser.add_argument("--output-dir", type=str, default="plots")
    parser.add_argument("--combined", action="store_true")
    parser.add_argument("--plotly", action="store_true")
    args = parser.parse_args()

    return run_workflow(args)


if __name__ == "__main__":
    exit(main())
