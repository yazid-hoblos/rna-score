#!/usr/bin/env python3
"""
Build KDE-based residue-pair scoring tables using R's density() via rpy2.
Reads *_kde_raw.txt from extract_distances.py and writes freq/score CSVs.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

# Ordered residue pair keys (including the XX reference)
PAIR_KEYS = ["AA", "AC", "AG", "AU", "CC", "CG", "CU", "GG", "GU", "UU", "XX"]


def load_raw_vectors(raw_dir: Path) -> Dict[str, np.ndarray]:
    """Load raw distance vectors dumped by extract_distances.py --method kde."""
    data: Dict[str, np.ndarray] = {}
    for key in PAIR_KEYS:
        path = raw_dir / f"{key}_kde_raw.txt"
        if path.exists() and path.stat().st_size > 0:
            data[key] = np.loadtxt(path, dtype=float, ndmin=1)
        else:
            data[key] = np.array([], dtype=float)
    return data


def _r_density(
    samples: np.ndarray,
    grid: np.ndarray,
    kernel: str,
    bandwidth: str | float | None,
    adjust: float,
    cut: float,
    weights: np.ndarray | None,
) -> Tuple[np.ndarray, np.ndarray]:
    # Use R's density() for flexible kernels/bandwidth rules
    from rpy2.robjects import FloatVector
    from rpy2.robjects.packages import importr

    stats = importr("stats")
    kwargs = {
        "kernel": kernel,
        "n": len(grid),
        "from": float(grid[0]),
        "to": float(grid[-1]),
        "cut": cut,
        "adjust": adjust,
    }
    if bandwidth is not None:
        kwargs["bw"] = bandwidth
    if weights is not None:
        kwargs["weights"] = FloatVector(weights.tolist())

    res = stats.density(FloatVector(samples.tolist()), **kwargs)
    x = np.array(res.rx2("x"), dtype=float)
    y = np.array(res.rx2("y"), dtype=float)

    # Normalize to guard against tiny numeric drift
    area = np.trapz(y, x)
    if area > 0:
        y = y / area
    return x, y


def _py_density(samples: np.ndarray, grid: np.ndarray, bandwidth: float) -> Tuple[np.ndarray, np.ndarray]:
    # SciPy fallback (Gaussian kernel only)
    kde = gaussian_kde(samples, bw_method=bandwidth)
    y = kde(grid)
    area = np.trapz(y, grid)
    if area > 0:
        y = y / area
    return grid, y


def compute_density(
    samples: np.ndarray,
    grid: np.ndarray,
    backend: str,
    kernel: str,
    bandwidth: str | float | None,
    adjust: float,
    cut: float,
) -> Tuple[np.ndarray, np.ndarray]:
    samples = samples.astype(float)
    if samples.size == 0:
        return grid, np.zeros_like(grid)

    if backend == "r":
        return _r_density(samples, grid, kernel, bandwidth, adjust, cut, None)

    # SciPy does not support the R kernel zoo; bandwidth is either rule name or numeric
    bw = bandwidth if isinstance(bandwidth, (int, float)) else "scott"
    return _py_density(samples, grid, float(adjust) * (1.0 if bw == "scott" else float(bw)))


def build_tables(
    raw_dir: Path,
    grid: np.ndarray,
    backend: str,
    kernel: str,
    bandwidth: str | float | None,
    adjust: float,
    cut: float,
    min_density: float,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
    # Load raw distances and compute reference density
    dist_map = load_raw_vectors(raw_dir)
    ref_x, ref_y = compute_density(dist_map["XX"], grid, backend, kernel, bandwidth, adjust, cut)
    ref_safe = np.clip(ref_y, min_density, None)

    pair_tables: Dict[str, pd.DataFrame] = {}
    for key in PAIR_KEYS:
        if key == "XX":
            continue
        x, y = compute_density(dist_map[key], grid, backend, kernel, bandwidth, adjust, cut)
        safe_y = np.clip(y, min_density, None)
        score = -np.log(safe_y / ref_safe)
        pair_tables[key] = pd.DataFrame(
            {
                "Distance": x,
                "Pair_Density": y,
                "Ref_Density": ref_y,
                "Score": score,
            }
        )

    ref_table = pd.DataFrame({"Distance": ref_x, "Density": ref_y})
    return pair_tables, ref_table


def save_tables(pair_tables: Dict[str, pd.DataFrame], ref_table: pd.DataFrame, output_dir: Path) -> None:
    """Persist frequency and score tables to CSV."""
    output_dir.mkdir(parents=True, exist_ok=True)
    for key, table in pair_tables.items():
        table[["Distance", "Pair_Density", "Ref_Density"]].to_csv(output_dir / f"freq_{key}.csv", index=False)
        table[["Distance", "Score"]].to_csv(output_dir / f"score_{key}.csv", index=False)
    ref_table.to_csv(output_dir / "freq_XX.csv", index=False)


def pick_backend(force_python: bool) -> str:
    # Prefer R + rpy2 when available; otherwise fall back to SciPy
    if force_python:
        return "python"
    try:
        import rpy2  # noqa: F401
    except ImportError:
        return "python"
    return "r"


def main() -> int:
    parser = argparse.ArgumentParser(description="Train KDE-based scoring tables from raw distance files.")
    parser.add_argument("--raw-dir", type=Path, required=True, help="Directory with *_kde_raw.txt from extract_distances.py")
    parser.add_argument("--output-dir", type=Path, default=Path("kde_output"), help="Where to write freq_*/score_* CSVs")
    parser.add_argument("--max-distance", type=float, default=20.0, help="Upper bound of the grid in Angstrom")
    parser.add_argument("--grid-step", type=float, default=0.1, help="Grid spacing for evaluating the density")
    parser.add_argument("--kernel", default="gaussian", help="R density() kernel (gaussian, epanechnikov, rect, triangular, cosine)")
    parser.add_argument(
        "--bandwidth",
        default="SJ-dpi",
        help="R density() bw method (nrd0, nrd, ucv, bcv, SJ-ste, SJ-dpi) or numeric bandwidth for SciPy fallback",
    )
    parser.add_argument("--adjust", type=float, default=1.0, help="Bandwidth multiplier passed to density()")
    parser.add_argument("--cut", type=float, default=3.0, help="How far the grid extends beyond data range (R density cut)")
    parser.add_argument("--min-density", type=float, default=1e-8, help="Floor to avoid log(0) in scoring")
    parser.add_argument("--python-backend", action="store_true", help="Force SciPy backend instead of R/rpy2")
    args = parser.parse_args()

    grid = np.arange(0.0, args.max_distance + args.grid_step, args.grid_step, dtype=float)
    backend = pick_backend(args.python_backend)

    pair_tables, ref_table = build_tables(
        raw_dir=args.raw_dir,
        grid=grid,
        backend=backend,
        kernel=args.kernel,
        bandwidth=args.bandwidth,
        adjust=args.adjust,
        cut=args.cut,
        min_density=args.min_density,
    )
    save_tables(pair_tables, ref_table, args.output_dir)

    meta = {
        "raw_dir": str(args.raw_dir),
        "output_dir": str(args.output_dir),
        "max_distance": args.max_distance,
        "grid_step": args.grid_step,
        "kernel": args.kernel,
        "bandwidth": args.bandwidth,
        "adjust": args.adjust,
        "cut": args.cut,
        "backend": backend,
    }
    (args.output_dir / "metadata.json").write_text(json.dumps(meta, indent=2))
    print(f"Wrote KDE tables to {args.output_dir} using {backend} backend")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
