"""
Shared utilities for RNA residue handling and chain extraction.
"""
from __future__ import annotations

from typing import Iterable, List, Optional, Tuple

import numpy as np

# Residue normalization and atom aliases
NUCLEOTIDE_MAP = {
    "A": "A",
    "ADE": "A",
    "G": "G",
    "GUA": "G",
    "C": "C",
    "CYT": "C",
    "U": "U",
    "URA": "U",
    "T": "U",  # occasionally T is used in RNA structures
}

ATOM_NAME_ALIASES = {
    "C3'": {"C3'", "C3*"},
    "C5'": {"C5'", "C5*"},
}

# Symmetric pair labels
PAIR_CANONICAL = {
    "AA": "AA",
    "AC": "AC",
    "AG": "AG",
    "AU": "AU",
    "CC": "CC",
    "CG": "CG",
    "CU": "UC",
    "GG": "GG",
    "GU": "UG",
    "UU": "UU",
}


def standardize_residue_name(resname: str) -> Optional[str]:
    """Return canonical nucleotide single-letter code or None if not RNA."""
    return NUCLEOTIDE_MAP.get(resname.strip().upper())


def canonical_pair_key(res1: str, res2: str) -> Optional[str]:
    """Return symmetric residue-pair key (e.g., AU/UA -> AU)."""
    r1 = standardize_residue_name(res1)
    r2 = standardize_residue_name(res2)
    if not r1 or not r2:
        return None
    return PAIR_CANONICAL.get("".join(sorted([r1, r2])))


def find_atom_coordinate(residue, target_names: Iterable[str]) -> Optional[np.ndarray]:
    """Locate the coordinate for the first matching atom name on a residue."""
    for atom in residue.get_atoms():
        if atom.get_name() in target_names:
            return atom.get_coord()
    return None


def extract_chain_residues(chain, target_names: Iterable[str]) -> List[Tuple[int, str, np.ndarray]]:
    """Collect sorted (resseq, res_type, coord) for matching residues in a chain."""
    residues: List[Tuple[int, str, np.ndarray]] = []
    for residue in chain:
        het_flag, resseq, _icode = residue.id
        if het_flag != " ":
            continue
        res_type = standardize_residue_name(residue.get_resname())
        if not res_type:
            continue
        coord = find_atom_coordinate(residue, target_names)
        if coord is None:
            continue
        residues.append((resseq, res_type, coord))
    residues.sort(key=lambda x: x[0])
    return residues
