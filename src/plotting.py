#!/usr/bin/env python3
"""
Visualize frequency distributions and scoring functions.

The script expects the CSV outputs of training.py (freq_*.csv, score_*.csv).
It generates per-pair plots for the frequency distribution, the scoring
profile, and a combined interaction plot overlaying both.

Usage example:
    python3 src/plotting.py --input-dir data/examples/training_output \\
        --output-dir data/examples/plots
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

import matplotlib.pyplot as plt
import pandas as pd


def find_pairs(input_dir: Path, only: Optional[Iterable[str]]) -> List[str]:
    """Return residue-pair keys present in the input directory."""
    if only:
        return [pair.upper() for pair in only]
    pairs = []
    for freq_file in sorted(input_dir.glob("freq_*.csv")):
        key = freq_file.stem.split("_", 1)[1]
        if key == "XX":
            continue
        pairs.append(key.upper())
    return pairs


def plot_frequency(df: pd.DataFrame, pair: str, output_dir: Path) -> None:
    """Plot frequency distribution for a residue pair."""
    width = float((df["Distance_Max"] - df["Distance_Min"]).iloc[0]) * 0.9
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(df["Distance_Mid"], df["Frequency"], width=width, color="#2a9d8f", edgecolor="black")
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Frequency")
    ax.set_title(f"Frequency distribution – {pair}")
    fig.tight_layout()
    output_path = output_dir / f"frequency_{pair}.png"
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_score(df: pd.DataFrame, pair: str, output_dir: Path) -> None:
    """Plot scoring function for a residue pair."""
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(df["Distance_Mid"], df["Score"], marker="o", color="#e76f51")
    ax.set_xlabel("Distance (Å)")
    ax.set_ylabel("Score")
    ax.set_title(f"Scoring function – {pair}")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    output_path = output_dir / f"scoring_{pair}.png"
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def plot_interaction(
    freq_df: pd.DataFrame,
    score_df: pd.DataFrame,
    pair: str,
    output_dir: Path,
) -> None:
    """Plot combined frequency and score interaction profile."""
    width = float((freq_df["Distance_Max"] - freq_df["Distance_Min"]).iloc[0]) * 0.9
    fig, ax1 = plt.subplots(figsize=(9, 4))
    ax1.bar(
        freq_df["Distance_Mid"],
        freq_df["Frequency"],
        width=width,
        color="#90caf9",
        edgecolor="black",
        alpha=0.6,
    )
    ax1.set_xlabel("Distance (Å)")
    ax1.set_ylabel("Frequency", color="#1e88e5")
    ax1.tick_params(axis="y", labelcolor="#1e88e5")

    ax2 = ax1.twinx()
    ax2.plot(score_df["Distance_Mid"], score_df["Score"], marker="s", color="#d81b60")
    ax2.set_ylabel("Score", color="#d81b60")
    ax2.tick_params(axis="y", labelcolor="#d81b60")

    fig.suptitle(f"Interaction profile – {pair}")
    fig.tight_layout()
    output_path = output_dir / f"interaction_{pair}.png"
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint for plotting frequency and scoring profiles."""
    parser = argparse.ArgumentParser(description="Plot frequency and scoring profiles.")
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("training_output"),
        help="Directory with freq_*.csv and score_*.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("plots"),
        help="Where to store generated images",
    )
    parser.add_argument("--pairs", nargs="*", help="Subset of pairs to plot (e.g., AU CG GG)")
    args = parser.parse_args(argv)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pairs = find_pairs(args.input_dir, args.pairs)
    if not pairs:
        raise SystemExit("No pair files found in input directory.")

    for pair in pairs:
        freq_path = args.input_dir / f"freq_{pair}.csv"
        score_path = args.input_dir / f"score_{pair}.csv"
        if not freq_path.exists() or not score_path.exists():
            print(f"Skipping {pair}: missing freq or score file.")
            continue
        freq_df = pd.read_csv(freq_path)
        score_df = pd.read_csv(score_path)
        plot_frequency(freq_df, pair, args.output_dir)
        plot_score(score_df, pair, args.output_dir)
        plot_interaction(freq_df, score_df, pair, args.output_dir)
        print(f"Generated plots for {pair}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
