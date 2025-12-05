"""
RNA Structure I/O Module

This module provides tools to fetch and parse RNA structures from PDB and mmCIF files.

Main Components:
----------------
1. OnlineFetcher:
   - Downloads PDB or mmCIF files from the RCSB Protein Data Bank (in-memory).
   - Returns the file content as a list of lines.

2. FastParser:
   - Efficiently parses PDB or mmCIF files.
   - Supports filtering by chain, atom type, and models.
   - Supports extracting specific atoms, all atoms, or centroids per residue.
   - Returns structured arrays containing coordinates, residue info, atom info, 
     chain IDs, model IDs, and B-factors.

Workflow Example:
-----------------
1. Fetch RNA structure:
    lines = OnlineFetcher.get_content('1ABC', file_format='pdb')

2. Parse structure:
    parser = FastParser(atom_mode="C3'")
    data = parser.parse(lines, file_format='pdb', pdb_name='1ABC')

Dependencies:
-------------
- numpy
- urllib
- collections
"""

import numpy as np
from collections import defaultdict
import urllib.request
import re


class OnlineFetcher:
    """
    Class for fetching PDB or mmCIF files from the RCSB Protein Data Bank online repository.
    """
    BASE_URL = "https://files.rcsb.org/download/"

    @staticmethod
    def get_content(pdb_id, file_format='pdb'):
        """
        Fetches the content of a PDB or mmCIF file for a given PDB ID.

        Parameters:
        -----------
        pdb_id : str
            The PDB identifier of the structure.
        file_format : str, optional
            Either 'pdb' or 'mmcif' to specify the file format. Default is 'pdb'.

        Returns:
        --------
        list[str] or None
            List of lines from the file if successful; None if download fails.
        """
        ext = "cif" if file_format == 'mmcif' else "pdb"
        url = f"{OnlineFetcher.BASE_URL}{pdb_id.upper()}.{ext}"
        
        try:
            with urllib.request.urlopen(url) as response:
                return response.read().decode('utf-8').splitlines()
        except Exception as e:
            return None


class FastParser:
    """
    Efficient parser for extracting atomic coordinates and related information from PDB or mmCIF files for RNA.

    Attributes:
    -----------
    atom_mode : str or set
        Defines which atoms to retain during parsing. Options:
        - specific atom name (e.g., "C3'")
        - 'centroid': compute centroid of all atoms per residue
        - 'all': keep all atoms
        - list of atom names
    allowed_bases : set
        Set of RNA nucleotide bases to keep ('A', 'U', 'C', 'G').
    """
    @staticmethod
    def _normalize_atom_name(atom_name):
        """
        Normalize atom names for consistent matching (removes quotes, whitespace, asterisks).
        """
        if not atom_name:
            return atom_name
        # Note: The original normalization included .replace("'", "'").
        # Standard PDB atom names like C3' use a prime/asterisk. This ensures it's handled.
        normalized = atom_name.strip().replace('"', '').replace("'", "'").replace('*', "'")
        return normalized
    
    def __init__(self, atom_mode=["C3'"]):
        self.atom_mode = atom_mode
        self.allowed_bases = {'A', 'U', 'C', 'G'}
        if isinstance(atom_mode, list):
            self.target_atoms = set([self._normalize_atom_name(a) for a in atom_mode])
            self.use_centroid = False
            self.keep_all = False
        elif atom_mode == "centroid":
            self.use_centroid = True
            self.keep_all = False
            self.target_atoms = None
        elif atom_mode == "all":
            self.keep_all = True
            self.use_centroid = False
            self.target_atoms = None
        else:
            self.target_atoms = {self._normalize_atom_name(atom_mode)}
            self.use_centroid = False
            self.keep_all = False

    def parse(self, lines, file_format='pdb', pdb_name='unknown', chains_filter=None, use_all_models=False):
        """
        Parse the structure file and extract atom coordinates.

        Parameters:
        -----------
        lines : list[str]
            Lines of the PDB or mmCIF file content.
        file_format : str, optional
            'pdb' or 'mmcif'. Default is 'pdb'.
        pdb_name : str, optional
            Name/ID of the structure for reference. Default is 'unknown'.
        chains_filter : list or None, optional
            List of chain IDs to include. If None, all chains are included.
        use_all_models : bool, optional
            Whether to include all NMR models. Default is False (first model only).

        Returns:
        --------
        dict or None
            Dictionary with parsed atom data, or None if parsing fails.
        """
        if file_format == 'pdb':
            return self._parse_pdb(lines, pdb_name, chains_filter, use_all_models)
        elif file_format == 'mmcif':
            return self._parse_cif(lines, pdb_name, chains_filter, use_all_models)
        else:
            return None

    def _parse_pdb(self, lines, pdb_name, chains_filter, use_all_models):
        """
        Parses a PDB file line by line.

        Returns:
        --------
        dict or None
            Dictionary containing parsed atom data or None if no atoms are valid.
        """
        coords = []
        res_names, res_ids, chain_ids = [], [], []
        atom_names, alt_locs, b_factors = [], [], []
        model_ids = []
        
        current_model = 1
        centroid_buffer = defaultdict(list)
        
        for line in lines:
            if line.startswith("MODEL"):
                try:
                    current_model = int(line.split()[1])
                except:
                    pass
                if not use_all_models and current_model > 1:
                    break
                continue
            
            if line.startswith("ENDMDL") and not use_all_models:
                break
            
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            
            atom_name = self._normalize_atom_name(line[12:16].strip())
            alt_loc = line[16].strip()
            res_name = line[17:20].strip()
            chain = line[21].strip()
            
            # Exclude hydrogen atoms (common convention is for atom name to start with 'H')
            if atom_name.startswith('H'):
                continue

            if res_name not in self.allowed_bases:
                continue
            if chains_filter and chain not in chains_filter:
                continue
            
            try:
                res_id = int(line[22:26].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                try:
                    b_factor = float(line[60:66].strip())
                except:
                    b_factor = 0.0
            except:
                continue
            
            key = (current_model, chain, res_id, res_name)
            
            if self.use_centroid:
                centroid_buffer[key].append([x, y, z, atom_name, alt_loc, b_factor])
            elif self.keep_all or atom_name in self.target_atoms:
                coords.append([x, y, z])
                res_names.append(res_name)
                res_ids.append(res_id)
                chain_ids.append(chain)
                atom_names.append(atom_name)
                alt_locs.append(alt_loc if alt_loc else '')
                b_factors.append(b_factor)
                model_ids.append(current_model)
        
        if self.use_centroid:
            for (model, chain, res_id, res_name), atoms in centroid_buffer.items():
                arr = np.array([a[:3] for a in atoms], dtype=float)
                centroid = arr.mean(axis=0)
                coords.append(centroid.tolist())
                res_names.append(res_name)
                res_ids.append(res_id)
                chain_ids.append(chain)
                atom_names.append('CENTROID')
                alt_locs.append('')
                b_factors.append(0.0)
                model_ids.append(model)

        if not coords:
            return None

        # Return structured data
        return {
            'pdb_name': pdb_name,
            'coords': np.array(coords, dtype=np.float32),
            'res_names': np.array(res_names, dtype='<U3'),
            'res_ids': np.array(res_ids, dtype=np.int32),
            'chain_ids': np.array(chain_ids, dtype='<U4'),
            'atom_names': np.array(atom_names, dtype='<U4'),
            'alt_locs': np.array(alt_locs, dtype='<U1'),
            'b_factors': np.array(b_factors, dtype=np.float32),
            'model_ids': np.array(model_ids, dtype=np.int32)
        }

    def _parse_cif(self, lines, pdb_name, chains_filter, use_all_models):
        """
        Parses a mmCIF file line by line.

        Returns:
        --------
        dict or None
            Dictionary containing parsed atom data or None if no atoms are valid.
        """
        coords = []
        res_names, res_ids, chain_ids = [], [], []
        atom_names, alt_locs, b_factors = [], [], []
        model_ids = []

        centroid_buffer = defaultdict(list)

        in_atom_site = False
        headers = []
        indices = {}

        # Debugging removed
        for line in lines:
            if line.startswith("loop_"):
                in_atom_site = False
                headers = []
                continue

            if line.startswith("_atom_site."):
                in_atom_site = True
                headers.append(line.strip())
                continue

            if in_atom_site:
                if line.startswith("_"):
                    in_atom_site = False
                    continue

                # Initialize header indices mapping
                if not indices and headers:
                    for i, h in enumerate(headers):
                        indices[h] = i
                    req = ['_atom_site.group_PDB', '_atom_site.label_comp_id',
                           '_atom_site.auth_asym_id', '_atom_site.auth_seq_id',
                           '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']
                    if not all(r in indices for r in req):
                        return None

                parts = line.split()
                if len(parts) < len(headers):
                    continue

                group = parts[indices['_atom_site.group_PDB']]
                if group not in ('ATOM', 'HETATM'):
                    continue

                res_name = parts[indices['_atom_site.label_comp_id']].strip()
                if res_name not in self.allowed_bases:
                    continue

                chain = parts[indices['_atom_site.auth_asym_id']]
                if chains_filter and chain not in chains_filter:
                    continue

                try:
                    res_id = int(parts[indices['_atom_site.auth_seq_id']])
                    x = float(parts[indices['_atom_site.Cartn_x']])
                    y = float(parts[indices['_atom_site.Cartn_y']])
                    z = float(parts[indices['_atom_site.Cartn_z']])

                    atom_name_raw = parts[indices.get('_atom_site.label_atom_id', indices.get('_atom_site.auth_atom_id', 0))]
                    atom_name = self._normalize_atom_name(atom_name_raw)
                    alt_loc = parts[indices.get('_atom_site.label_alt_id', 0)] if '_atom_site.label_alt_id' in indices else '.'
                    b_factor = float(parts[indices['_atom_site.B_iso_or_equiv']]) if '_atom_site.B_iso_or_equiv' in indices else 0.0
                    model_num = int(parts[indices['_atom_site.pdbx_PDB_model_num']]) if '_atom_site.pdbx_PDB_model_num' in indices else 1

                    # Exclude hydrogen atoms (common convention is for atom name to start with 'H')
                    if atom_name.startswith('H'):
                        continue

                    if alt_loc == '.':
                        alt_loc = ''

                    if not use_all_models and model_num > 1:
                        continue

                except (ValueError, IndexError):
                    continue
                key = (model_num, chain, res_id, res_name)
                if self.use_centroid:
                    centroid_buffer[key].append([x, y, z, atom_name, alt_loc, b_factor])
                elif self.keep_all or atom_name in self.target_atoms:
                    coords.append([x, y, z])
                    res_names.append(res_name)
                    res_ids.append(res_id)
                    chain_ids.append(chain)
                    atom_names.append(atom_name)
                    alt_locs.append(alt_loc)
                    b_factors.append(b_factor)
                    model_ids.append(model_num)

        if self.use_centroid:
            for (model, chain, res_id, res_name), atoms in centroid_buffer.items():
                arr = np.array([a[:3] for a in atoms], dtype=float)
                centroid = arr.mean(axis=0)
                coords.append(centroid.tolist())
                res_names.append(res_name)
                res_ids.append(res_id)
                chain_ids.append(chain)
                atom_names.append('CENTROID')
                alt_locs.append('')
                b_factors.append(0.0)
                model_ids.append(model)

        if not coords:
            return None

        return {
            'pdb_name': pdb_name,
            'coords': np.array(coords, dtype=np.float32),
            'res_names': np.array(res_names, dtype='<U3'),
            'res_ids': np.array(res_ids, dtype=np.int32),
            'chain_ids': np.array(chain_ids, dtype='<U4'),
            'atom_names': np.array(atom_names, dtype='<U4'),
            'alt_locs': np.array(alt_locs, dtype='<U1'),
            'b_factors': np.array(b_factors, dtype=np.float32),
            'model_ids': np.array(model_ids, dtype=np.int32)
        }