#!/usr/bin/env python3

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
# from matplotlib.cm import get_cmap

def _ensure_dir(path):
    os.makedirs(path, exist_ok=True)

PAIRS = ['AA', 'AC', 'AG', 'AU', 'CC', 'CG', 'GG', 'UC', 'UG', 'UU']

plasma = plt.get_cmap("plasma")
PAIR_COLORS = {pair: plasma(i / (len(PAIRS)-1)) for i, pair in enumerate(PAIRS)}

def load_freq_table(filepath):
    if not os.path.exists(filepath):
        return None
    data = np.loadtxt(filepath, delimiter=",", skiprows=1)
    return {"value": data[:, 0], "freq": data[:, 1]}

def plot_single_hist(pair, data, outpath):
    sns.set_theme(style="whitegrid")
    x = data["value"]
    f = data["freq"]
    
    # -- orangy pinkish color of plasma colormap
    current_color = plasma(0.75)

    plt.figure(figsize=(7, 5))
    plt.bar(x, f, color=current_color, alpha=0.85)
    plt.title(f"Frequency Histogram: {pair}", fontsize=14)
    plt.xlabel("Value", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()

def plot_combined_hist(pairs_data, outpath):
    sns.set_theme(style="whitegrid")
    plt.figure(figsize=(10, 6))

    for pair, data in pairs_data.items():
        plt.plot(
            data["value"],
            data["freq"],
            label=pair,
            color=PAIR_COLORS[pair],
            linewidth=2
        )

    plt.title("Combined Histograms")
    plt.xlabel("Value")
    plt.ylabel("Frequency")
    plt.legend()
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()

def rgba_to_plotly(rgba):
    r, g, b, a = rgba
    return f"rgba({int(r*255)}, {int(g*255)}, {int(b*255)}, {a})"

def plot_single_hist_plotly(pair, data, outpath):
    current_color=plasma(0.75)

    fig = go.Figure()
    fig.add_bar(
        x=data["value"],
        y=data["freq"],
        marker=dict(color=rgba_to_plotly(current_color), opacity=0.85),
        name=pair
    )
    fig.update_layout(
        title=f"Frequency Histogram: {pair}",
        xaxis_title="Value",
        yaxis_title="Frequency",
        template="plotly_white"
    )
    fig.write_html(outpath)

def plot_combined_hist_plotly(pairs_data, outpath):
    fig = go.Figure()

    for pair, data in pairs_data.items():
        fig.add_scatter(
            x=data["value"],
            y=data["freq"],
            mode="lines",
            name=pair,
            line=dict(color=rgba_to_plotly(PAIR_COLORS[pair]), width=3)
        )

    fig.update_layout(
        title="Combined Histograms",
        xaxis_title="Value",
        yaxis_title="Frequency",
        template="plotly_white"
    )
    fig.write_html(outpath)

def run_workflow(args):
    print(f"Plotting scores distributions from {args.input_dir} …")
    _ensure_dir(args.output_dir)
    _ensure_dir(os.path.join(args.output_dir, "html"))
    _ensure_dir(os.path.join(args.output_dir, "png"))

    pairs_data = {}
    for pair in PAIRS:
        filepath = os.path.join(args.input_dir, f"freq_{pair}.csv")
        data = load_freq_table(filepath)
        if data is None:
            print(f"  [!WARNING] missing file for {pair}")
            continue

        pairs_data[pair] = data
        out_png = os.path.join(args.output_dir, "png", f"{pair}_dist.png")
        plot_single_hist(pair, data, out_png)
        print(f" ✓ {pair}")

    for pair, data in pairs_data.items():
        out_html = os.path.join(args.output_dir, "html", f"{pair}_dist.html")
        plot_single_hist_plotly(pair, data, out_html)

    if args.combined and pairs_data:
        out_png = os.path.join(args.output_dir, "png", "all_dist.png")
        plot_combined_hist(pairs_data, out_png)
        print(" ✓ Combined histogram saved.")

        if args.combined:
            out_html = os.path.join(args.output_dir, "html", "all_dist.html")
            plot_combined_hist_plotly(pairs_data, out_html)

        print(" ✓ Plotly HTML histograms generated.")

    print("Done")
    print('✓ Complete - distrbutions found in "' + args.output_dir + '"')
    print('\tplotly intreactive plots in "' + os.path.join(args.output_dir, "html") + '"')
    print('\tstatic png plots in "' + os.path.join(args.output_dir, "png") + '"')
    return 0

def main():
    parser = argparse.ArgumentParser(description="Plot histogram frequency tables")
    parser.add_argument("--input-dir", type=str, default="freq_tables")
    parser.add_argument("--output-dir", type=str, default="hist_plots")
    parser.add_argument("--combined", action="store_true")
    parser.add_argument("--plotly", action="store_true")
    args = parser.parse_args()
    return run_workflow(args)

if __name__ == "__main__":
    exit(main())

