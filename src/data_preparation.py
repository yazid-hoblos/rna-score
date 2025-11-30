#!/usr/bin/env python3
"""
Prepare RNA structural data by extracting inter-residue distances.

The script downloads or reads PDB/mmCIF files, keeps only C3'/C5' atoms,
computes intrachain distances separated by at least `--min-separation`,
filters out distances above `--max-distance`, and writes a CSV with one row
per residue pair.
"""
from __future__ import annotations

import argparse
import csv
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np
import requests
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure

from rna_utils import ATOM_NAME_ALIASES, extract_chain_residues


@dataclass
class PDBEntry:
    """Identifier for a PDB entry with an optional chain or local path."""

    identifier: str
    chain: Optional[str] = None
    path: Optional[Path] = None


def parse_input_list(path: Path) -> List[PDBEntry]:
    """Parse CSV or whitespace list of PDB IDs with optional chain column."""
    entries: List[PDBEntry] = []
    text = path.read_text().splitlines()
    has_comma = any("," in line for line in text)
    for raw in text:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts: Sequence[str]
        if has_comma:
            parts = [p.strip() for p in line.split(",") if p.strip()]
        else:
            parts = line.split()
        if not parts:
            continue
        identifier = parts[0]
        chain = parts[1] if len(parts) > 1 else None
        path_candidate = Path(identifier)
        entry_path = path_candidate if path_candidate.exists() else None
        entries.append(PDBEntry(identifier=identifier, chain=chain, path=entry_path))
    return entries


def download_structure(pdb_id: str, cache_dir: Optional[Path], prefer_format: str) -> Path:
    """Download a PDB/mmCIF file and return the local path (cached if present)."""
    pdb_id = pdb_id.lower()
    fmt = prefer_format.lower()
    if fmt not in {"pdb", "cif", "mmcif", "auto"}:
        raise ValueError(f"Unsupported format {prefer_format}")
    ext = "cif" if fmt in {"cif", "mmcif", "auto"} else "pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.{ext}"
    if cache_dir:
        cache_dir.mkdir(parents=True, exist_ok=True)
        local_path = cache_dir / f"{pdb_id}.{ext}"
        if local_path.exists():
            return local_path
    else:
        local_path = Path(f"{pdb_id}.{ext}")
    logging.info("Downloading %s to %s", url, local_path)
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    local_path.write_bytes(response.content)
    return local_path


def parse_structure(path: Path) -> Structure:
    """Parse a local structure file into a Bio.PDB Structure object."""
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(QUIET=True)
    elif suffix == ".pdb":
        parser = PDBParser(QUIET=True)
    else:
        raise ValueError(f"Unsupported structure format for {path}")
    structure = parser.get_structure(path.stem, str(path))
    if structure is None or not isinstance(structure, Structure):
        raise ValueError(f"Failed to parse structure from {path}")
    return structure


def load_structure(entry: PDBEntry, cache_dir: Optional[Path], prefer_format: str) -> Structure:
    """Load a structure either from a local path or by downloading from RCSB."""
    if entry.path:
        return parse_structure(entry.path)
    path = download_structure(entry.identifier, cache_dir, prefer_format)
    return parse_structure(path)


def collect_chain_records(
    structure: Structure,
    chain_id: str,
    target_names: Iterable[str],
    max_distance: float,
    min_separation: int,
    pdb_id: str,
) -> List[Tuple]:
    """Return distance rows for a single chain."""
    # pylint: disable=too-many-arguments,too-many-positional-arguments,too-many-locals
    residues = extract_chain_residues(structure[0][chain_id], target_names)

    records: List[Tuple] = []
    for i, (res_i_idx, res_i_type, coord_i) in enumerate(residues):
        for res_j_idx, res_j_type, coord_j in residues[i + min_separation :]:
            distance = float(np.linalg.norm(coord_i - coord_j))
            if distance > max_distance:
                continue
            records.append(
                (
                    pdb_id,
                    res_i_type,
                    chain_id,
                    res_i_idx,
                    coord_i[0],
                    coord_i[1],
                    coord_i[2],
                    res_j_type,
                    chain_id,
                    res_j_idx,
                    coord_j[0],
                    coord_j[1],
                    coord_j[2],
                    distance,
                )
            )
    return records


def extract_distances(
    structure: Structure,
    entry: PDBEntry,
    atom_choice: str,
    max_distance: float,
    min_separation: int,
) -> List[Tuple]:
    """Extract all qualifying distance rows for an entry across one or more chains."""
    target_names = ATOM_NAME_ALIASES[atom_choice]
    pdb_id = entry.identifier.upper()
    chain_ids = [entry.chain] if entry.chain else [c.id for c in structure[0]]
    records: List[Tuple] = []
    for chain_id in chain_ids:
        try:
            records.extend(
                collect_chain_records(
                    structure,
                    chain_id,
                    target_names,
                    max_distance,
                    min_separation,
                    pdb_id,
                )
            )
        except KeyError:
            logging.warning("Chain %s not found in %s", chain_id, pdb_id)
    return records


def write_output(rows: List[Tuple], output_path: Path) -> None:
    """Write the aggregated rows to CSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "PDB_ID",
        "Residue_Type_1",
        "Chain_ID_1",
        "Residue_Seq_Num_1",
        "X_Coord_1",
        "Y_Coord_1",
        "Z_Coord_1",
        "Residue_Type_2",
        "Chain_ID_2",
        "Residue_Seq_Num_2",
        "X_Coord_2",
        "Y_Coord_2",
        "Z_Coord_2",
        "Distance",
    ]
    with output_path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entrypoint for preparing distance rows from structures."""
    parser = argparse.ArgumentParser(description="Prepare RNA distance dataset.")
    parser.add_argument(
        "--ids-file",
        required=True,
        type=Path,
        help="Text/CSV file containing PDB IDs (optionally followed by chain).",
    )
    parser.add_argument(
        "--output",
        default=Path("prepared_distances.csv"),
        type=Path,
        help="Destination CSV path.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path(".cache/pdb"),
        help="Directory to cache downloaded structures.",
    )
    parser.add_argument(
        "--atom",
        choices=list(ATOM_NAME_ALIASES.keys()),
        default="C3'",
        help="Atom to use for distance calculations.",
    )
    parser.add_argument(
        "--max-distance",
        type=float,
        default=20.0,
        help="Maximum distance (Å) to retain.",
    )
    parser.add_argument(
        "--min-separation",
        type=int,
        default=4,
        help="Minimum sequence separation between residues (i, i+N).",
    )
    parser.add_argument(
        "--structure-format",
        choices=["auto", "pdb", "cif", "mmcif"],
        default="auto",
        help="Preferred download format when fetching from RCSB.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (DEBUG, INFO, WARNING, ERROR).",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(levelname)s: %(message)s",
    )

    entries = parse_input_list(args.ids_file)
    if not entries:
        logging.error("No entries parsed from %s", args.ids_file)
        return 1

    all_rows: List[Tuple] = []
    for entry in entries:
        try:
            structure = load_structure(entry, args.cache_dir, args.structure_format)
        except (requests.RequestException, OSError, ValueError) as exc:
            logging.error("Failed to load %s: %s", entry.identifier, exc)
            continue
        rows = extract_distances(
            structure,
            entry,
            args.atom,
            args.max_distance,
            args.min_separation,
        )
        logging.info(
            "Processed %s (chains: %s) -> %d distance rows",
            entry.identifier,
            entry.chain if entry.chain else "all",
            len(rows),
        )
        all_rows.extend(rows)

    if not all_rows:
        logging.warning("No distance rows produced; check inputs or parameters.")
    write_output(all_rows, args.output)
    logging.info("Wrote %d rows to %s", len(all_rows), args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
