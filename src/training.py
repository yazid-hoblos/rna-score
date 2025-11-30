#!/usr/bin/env python3
"""
Build residue-pair scoring tables from prepared distance data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rna_utils import PAIR_CANONICAL, canonical_pair_key

PAIR_KEYS = list(PAIR_CANONICAL.values())


def build_bins(max_distance: float, bin_size: float) -> np.ndarray:
    """Return bin edges from 0 to max_distance inclusive."""
    edges = np.arange(0, max_distance + bin_size, bin_size, dtype=float)
    if edges[-1] < max_distance:
        edges = np.append(edges, max_distance)
    return edges


def histogram(distances: Iterable[float], bins: np.ndarray) -> np.ndarray:
    """Count distances into provided bins."""
    distances_arr = np.fromiter(distances, dtype=float)
    if distances_arr.size == 0:
        return np.zeros(len(bins) - 1, dtype=float)
    counts, _ = np.histogram(distances_arr, bins=bins)
    return counts.astype(float)


def compute_tables(
    df: pd.DataFrame, bins: np.ndarray, max_score: float, pseudocount: float
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    """Compute frequency and scoring tables for each residue pair and reference."""
    # pylint: disable=too-many-locals
    distances_ref: List[float] = []
    pair_distance_map: Dict[str, List[float]] = {key: [] for key in PAIR_KEYS}
    for _, row in df.iterrows():
        key = canonical_pair_key(row["Residue_Type_1"], row["Residue_Type_2"])
        if not key:
            continue
        distance = float(row["Distance"])
        if distance < 0 or distance > bins[-1]:
            continue
        pair_distance_map[key].append(distance)
        distances_ref.append(distance)

    ref_counts = histogram(distances_ref, bins)
    ref_total = ref_counts.sum()
    ref_freq = (ref_counts + pseudocount) / (ref_total + pseudocount * len(ref_counts))

    pair_tables: Dict[str, pd.DataFrame] = {}
    bin_left = bins[:-1]
    bin_right = bins[1:]
    bin_mid = bin_left + (bin_right - bin_left) / 2

    for key, dists in pair_distance_map.items():
        counts = histogram(dists, bins)
        total = counts.sum()
        if total > 0:
            denom = total + pseudocount * len(counts)
            freq = (counts + pseudocount) / denom
        else:
            freq = np.zeros_like(counts)
        safe_ref = np.where(ref_freq <= 0, pseudocount, ref_freq)
        safe_freq = np.where(freq <= 0, pseudocount, freq)
        score = -np.log(np.divide(safe_freq, safe_ref))
        score = np.clip(score, None, max_score)
        pair_tables[key] = pd.DataFrame(
            {
                "Distance_Min": bin_left,
                "Distance_Max": bin_right,
                "Distance_Mid": bin_mid,
                "Count": counts,
                "Frequency": freq,
                "Score": score,
            }
        )

    ref_table = pd.DataFrame(
        {
            "Distance_Min": bin_left,
            "Distance_Max": bin_right,
            "Distance_Mid": bin_mid,
            "Count": ref_counts,
            "Frequency": ref_freq,
        }
    )
    return pair_tables, ref_table


def save_tables(
    pair_tables: Dict[str, pd.DataFrame],
    ref_table: pd.DataFrame,
    output_dir: Path,
) -> None:
    """Persist frequency/score tables to CSV files in the output directory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for key, table in pair_tables.items():
        freq_path = output_dir / f"freq_{key}.csv"
        score_path = output_dir / f"score_{key}.csv"
        table[["Distance_Min", "Distance_Max", "Distance_Mid", "Frequency"]].to_csv(
            freq_path,
            index=False,
        )
        table[["Distance_Min", "Distance_Max", "Distance_Mid", "Score"]].to_csv(
            score_path,
            index=False,
        )
    ref_table.to_csv(output_dir / "freq_XX.csv", index=False)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint for training scoring tables."""
    parser = argparse.ArgumentParser(description="Train residue-pair scoring tables.")
    parser.add_argument("--input", required=True, type=Path, help="CSV from data_preparation.py")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("training_output"),
        help="Where to store tables",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=20.0,
        help="Max distance used for bins",
    )
    parser.add_argument("--bin-size", type=float, default=1.0, help="Bin width in Å")
    parser.add_argument(
        "--pseudocount",
        type=float,
        default=1e-4,
        help="Additive smoothing to avoid zeros",
    )
    parser.add_argument(
        "--max-score",
        type=float,
        default=10.0,
        help="Upper bound for pseudo-energy",
    )
    args = parser.parse_args(argv)

    df = pd.read_csv(args.input)
    if df.empty:
        raise SystemExit("Input CSV is empty; run data_preparation first.")

    bins = build_bins(args.max_distance, args.bin_size)
    pair_tables, ref_table = compute_tables(df, bins, args.max_score, args.pseudocount)
    save_tables(pair_tables, ref_table, args.output_dir)

    metadata = {
        "input": str(args.input),
        "output_dir": str(args.output_dir),
        "max_distance": args.max_distance,
        "bin_size": args.bin_size,
        "pairs": PAIR_KEYS,
    }
    (args.output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2))
    print(f"Wrote frequency and score tables to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
