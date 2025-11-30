#!/usr/bin/env python3
"""
Score a predicted RNA structure using pre-trained residue-pair tables.
"""
from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from data_preparation import PDBEntry, load_structure
from rna_utils import ATOM_NAME_ALIASES, canonical_pair_key, extract_chain_residues


def load_score_tables(score_dir: Path) -> Dict[str, pd.DataFrame]:
    """Load all score_*.csv tables from a directory."""
    tables: Dict[str, pd.DataFrame] = {}
    for score_file in score_dir.glob("score_*.csv"):
        pair = score_file.stem.split("_", 1)[1].upper()
        if pair == "XX":
            continue
        tables[pair] = pd.read_csv(score_file)
    if not tables:
        raise FileNotFoundError(f"No score_*.csv files found in {score_dir}")
    return tables


def interpolate_score(distance: float, table: pd.DataFrame) -> float:
    """Interpolate a score for a specific distance using mid-bin values."""
    mids = table["Distance_Mid"].to_numpy()
    scores = table["Score"].to_numpy()
    return float(np.interp(distance, mids, scores, left=scores[0], right=scores[-1]))


def collect_distances(
    structure,
    chain_ids: Iterable[str],
    atom_choice: str,
    max_distance: float,
    min_separation: int,
) -> List[Tuple[str, str, float]]:
    """Collect all qualifying residue-residue distances for scoring."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    target_names = ATOM_NAME_ALIASES[atom_choice]
    distances: List[Tuple[str, str, float]] = []
    model = structure[0]
    for chain_id in chain_ids:
        if chain_id not in model:
            logging.warning("Chain %s not found; skipping", chain_id)
            continue
        residues = extract_chain_residues(model[chain_id], target_names)
        for i, (_, res_i, coord_i) in enumerate(residues):
            for _, res_j, coord_j in residues[i + min_separation :]:
                dist = float(np.linalg.norm(coord_i - coord_j))
                if dist > max_distance:
                    continue
                distances.append((res_i, res_j, dist))
    return distances


def score_structure(
    structure,
    tables: Dict[str, pd.DataFrame],
    atom_choice: str,
    max_distance: float,
    min_separation: int,
    chains: Optional[List[str]] = None,
) -> float:
    """Compute total pseudo-energy score for a structure."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments
    model = structure[0]
    chain_ids = chains if chains else [c.id for c in model]
    distances = collect_distances(structure, chain_ids, atom_choice, max_distance, min_separation)
    total_score = 0.0
    for res1, res2, dist in distances:
        key = canonical_pair_key(res1, res2)
        if not key or key not in tables:
            continue
        total_score += interpolate_score(dist, tables[key])
    return total_score


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint to score structures with pre-trained tables."""
    parser = argparse.ArgumentParser(description="Score a predicted RNA structure.")
    parser.add_argument("--structure", required=True, help="PDB ID or path to PDB/mmCIF file")
    parser.add_argument(
        "--scores-dir",
        type=Path,
        default=Path("training_output"),
        help="Directory with score_*.csv files",
    )
    parser.add_argument(
        "--atom",
        choices=list(ATOM_NAME_ALIASES.keys()),
        default="C3'",
        help="Atom used for scoring",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=20.0,
        help="Maximum distance considered (Å)",
    )
    parser.add_argument(
        "--min-separation",
        type=int,
        default=4,
        help="Minimum sequence separation (i, i+N)",
    )
    parser.add_argument("--chains", nargs="*", help="Chain IDs to score; default: all chains")
    parser.add_argument("--output", type=Path, help="Optional CSV path to append score result")
    parser.add_argument("--log-level", default="INFO", help="Logging level")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(levelname)s: %(message)s",
    )
    tables = load_score_tables(args.scores_dir)

    struct_path = Path(args.structure)
    entry = PDBEntry(identifier=args.structure, path=struct_path if struct_path.exists() else None)
    try:
        structure = load_structure(entry, cache_dir=Path(".cache/pdb"), prefer_format="auto")
    except (OSError, ValueError) as exc:
        logging.error("Failed to load structure: %s", exc)
        return 1

    available_chains = [c.id for c in structure[0]]
    requested_chains = args.chains if args.chains else None
    if requested_chains:
        requested_chains = [str(c) for c in requested_chains]
        missing = [c for c in requested_chains if c not in available_chains]
        if len(missing) == len(requested_chains):
            logging.error(
                "None of the requested chains %s were found. Available: %s",
                ",".join(requested_chains),
                ",".join(available_chains),
            )
            return 1
        if missing:
            logging.warning(
                "Skipping missing chains %s. Available chains: %s",
                ",".join(missing),
                ",".join(available_chains),
            )

    total = score_structure(
        structure,
        tables,
        atom_choice=args.atom,
        max_distance=args.max_distance,
        min_separation=args.min_separation,
        chains=requested_chains,
    )
    pdb_id = Path(args.structure).stem.upper()
    logging.info("Total score for %s: %.4f", pdb_id, total)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        is_new = not args.output.exists()
        with args.output.open("a") as f:
            if is_new:
                f.write("PDB_ID,Total_Score\n")
            f.write(f"{pdb_id},{total:.6f}\n")
        logging.info("Appended result to %s", args.output)
    else:
        print(f"{pdb_id},{total:.6f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
