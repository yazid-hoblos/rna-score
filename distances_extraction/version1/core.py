"""
RNA STRUCTURE ANALYSIS CORE LIBRARY
===================================

A specialized library for high-throughput processing of RNA 3D structures. 
This module provides optimized tools for parsing PDB/mmCIF files and computing 
pairwise atomic distances using spatial indexing.

It is designed as a SHARED DEPENDENCY:
1. Training Phase: Used to extract raw statistics from the database.
2. Scoring Phase: Used to parse puzzle structures and compute energy.

Key Features:
    - **Vectorized Parsing**: Converts PDB/mmCIF text directly into NumPy arrays 
      (Structure-of-Arrays). This is significantly faster than standard parsers 
      because it avoids the overhead of creating Python objects for every atom, 
      while still retaining all essential metadata (B-factors, names, etc.).
    - **Spatial Optimization**: Utilizes KD-Trees for O(N log N) neighbor searching, 
      enabling efficient processing of large macromolecular assemblies.
    - **Flexible Atom Selection**: Supports extraction of specific atom (e.g., C3', P),
      centroid, subsets, or all-atom representations.
    - **In-Memory Fetching**: Capable of streaming structure data directly from RCSB PDB.

Classes:
    - OnlineFetcher: Handles retrieval of structure data from remote repositories.
    - FastParser: Parses raw PDB/mmCIF content into structured NumPy arrays.
    - DistanceComputer: Performs geometric calculations and filtering.

Usage:
    >>> from core import FastParser, DistanceComputer
    >>> parser = FastParser(atom_mode="C3'")
    >>> data = parser.parse(content, file_format='pdb', pdb_name='1EHZ')
    >>> computer = DistanceComputer(cutoff=20.0)
    >>> distances, _ = computer.compute(data, dist_mode='intra')
"""

import numpy as np
from scipy.spatial import cKDTree
from collections import defaultdict
import urllib.request
import os
import sys
import pandas as pd

class OnlineFetcher:
    """
    Provides static methods to retrieve structural data from remote databases.
    """
    BASE_URL = "https://files.rcsb.org/download/"

    @staticmethod
    def get_content(pdb_id, file_format='pdb'):
        """
        Retrieves structure content from RCSB PDB.

        Args:
            pdb_id (str): The 4-character PDB identifier (e.g., '1EHZ').
            file_format (str): Format specifier ('pdb' or 'mmcif').

        Returns:
            list[str] or None: File content split by lines, or None if request fails.
        """
        ext = "cif" if file_format == 'mmcif' else "pdb"
        url = f"{OnlineFetcher.BASE_URL}{pdb_id.upper()}.{ext}"
        
        try:
            with urllib.request.urlopen(url) as response:
                return response.read().decode('utf-8').splitlines()
        except Exception as e:
            print(f"Error fetching {pdb_id} from web: {e}")
            return None

class FastParser:
    """
    Optimized parser for macromolecular structure files.
    
    Implements a 'Structure-of-Arrays' approach, storing atomic data in parallel 
    NumPy arrays rather than object hierarchies. This design choice is optimized 
    for vectorized numerical operations and low memory overhead.
    """
    def __init__(self, atom_mode="C3'"):
        """
        Initialize the parser configuration.

        Args:
            atom_mode (str or list): Specifies which atoms to extract.
                - 'C3'': Standard RNA backbone representation.
                - 'centroid': Computes geometric center of residue.
                - 'all': Extracts all atoms.
                - str (e.g., 'P'): Extracts specific atom type.
                - list (e.g., ['P', "C3'"]): Extracts subset of atoms (Coarse Graining).
        """
        self.atom_mode = atom_mode
        self.allowed_bases = {'A', 'U', 'C', 'G'}
        
        # Optimize lookup if atom_mode is a list
        if isinstance(self.atom_mode, list):
            self.atom_mode = set(self.atom_mode)

    def parse(self, content_lines, file_format, pdb_name, chains_filter=None, use_all_models=False):
        """
        Parses structure content into a dictionary of arrays.

        Args:
            content_lines (list): Raw text lines of the file.
            file_format (str): 'pdb' or 'mmcif'.
            pdb_name (str): Identifier for the structure.
            chains_filter (list, optional): List of chain IDs to include. Defaults to None (all).
            use_all_models (bool, optional): If True, processes all models (NMR). Defaults to False.

        Returns:
            dict: Structured data containing 'coords', 'res_names', 'chain_ids', etc.
        """
        if file_format == 'mmcif':
            return self._parse_cif(content_lines, pdb_name, chains_filter, use_all_models)
        else:
            return self._parse_pdb(content_lines, pdb_name, chains_filter, use_all_models)

    def _should_keep_atom(self, atom_name):
        """Internal filter for atom selection strategy."""
        if self.atom_mode == 'all' or self.atom_mode == 'centroid':
            return True
        if isinstance(self.atom_mode, set):
            return atom_name in self.atom_mode
        return atom_name == self.atom_mode

    def _parse_pdb(self, lines, pdb_name, chains_filter, use_all_models):
        """Internal handler for PDB format parsing."""
        coords = []
        res_names, res_ids, chain_ids = [], [], []
        atom_names, alt_locs, b_factors = [], [], []
        model_ids = []
        
        centroid_buffer = defaultdict(list) 

        current_model = 1
        model_found = False
        
        for line in lines:
            if line.startswith("MODEL"):
                try:
                    current_model = int(line.split()[1])
                except:
                    current_model = 1
                if model_found and not use_all_models:
                    break 
                model_found = True
                
            if line.startswith("ENDMDL"):
                if not use_all_models:
                    break

            if not line.startswith("ATOM"):
                continue

            # Strict column-based parsing (PDB standard)
            alt_loc = line[16]
            if alt_loc != ' ' and alt_loc != 'A':
                continue
                
            res_name = line[17:20].strip()
            if res_name not in self.allowed_bases:
                continue
                
            chain_id = line[21]
            if chains_filter and chain_id not in chains_filter:
                continue

            atom_name = line[12:16].strip()
            
            if not self._should_keep_atom(atom_name):
                continue

            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                res_seq = int(line[22:26])
                bfactor = float(line[60:66])
            except ValueError:
                continue 

            def append_data(c, rn, ri, ci, an, al, bf, mid):
                coords.append(c)
                res_names.append(rn)
                res_ids.append(ri)
                chain_ids.append(ci)
                atom_names.append(an)
                alt_locs.append(al)
                b_factors.append(bf)
                model_ids.append(mid)

            if self.atom_mode == "centroid":
                key = (chain_id, res_seq, res_name, current_model)
                centroid_buffer[key].append([x, y, z])
            else:
                append_data([x, y, z], res_name, res_seq, chain_id, atom_name, alt_loc, bfactor, current_model)

        # Centroid Calculation
        if self.atom_mode == "centroid":
            for (chain, seq, rname, mid), atom_coords in centroid_buffer.items():
                coords.append(np.mean(atom_coords, axis=0))
                res_names.append(rname)
                res_ids.append(seq)
                chain_ids.append(chain)
                atom_names.append("CTR")
                alt_locs.append(".") 
                b_factors.append(0.0)
                model_ids.append(mid)

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

    def _parse_cif(self, lines, pdb_name, chains_filter, use_all_models):
        """Internal handler for mmCIF format parsing."""
        coords = []
        res_names, res_ids, chain_ids = [], [], []
        atom_names, alt_locs, b_factors = [], [], []
        model_ids = []
        
        centroid_buffer = defaultdict(list)
        
        in_atom_site = False
        headers = []
        indices = {}
        
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
                
                # Dynamically map column names to indices
                if not indices and headers:
                    for i, h in enumerate(headers):
                        indices[h] = i
                    req = ['_atom_site.group_PDB', '_atom_site.label_comp_id', 
                           '_atom_site.auth_asym_id', '_atom_site.auth_seq_id', 
                           '_atom_site.Cartn_x', '_atom_site.Cartn_y', '_atom_site.Cartn_z']
                    if not all(r in indices for r in req):
                        return None 

                parts = line.split()
                if len(parts) < len(indices): continue 

                if parts[indices['_atom_site.group_PDB']] != 'ATOM': continue
                
                model_num = 1
                if '_atom_site.pdbx_PDB_model_num' in indices:
                    model_num = int(parts[indices['_atom_site.pdbx_PDB_model_num']])
                    if not use_all_models and model_num > 1:
                        break 
                
                alt = '.'
                if '_atom_site.label_alt_id' in indices:
                    alt = parts[indices['_atom_site.label_alt_id']]
                    if alt != '.' and alt != 'A': continue

                res_name = parts[indices['_atom_site.label_comp_id']]
                if res_name not in self.allowed_bases: continue
                
                chain_id = parts[indices['_atom_site.auth_asym_id']]
                if chains_filter and chain_id not in chains_filter: continue
                
                atom_name = "UNK"
                if '_atom_site.label_atom_id' in indices:
                    atom_name = parts[indices['_atom_site.label_atom_id']].replace('"', '')
                
                if not self._should_keep_atom(atom_name):
                    continue
                
                bfactor = 0.0
                if '_atom_site.B_iso_or_equiv' in indices:
                    try:
                        bfactor = float(parts[indices['_atom_site.B_iso_or_equiv']])
                    except:
                        pass

                try:
                    x = float(parts[indices['_atom_site.Cartn_x']])
                    y = float(parts[indices['_atom_site.Cartn_y']])
                    z = float(parts[indices['_atom_site.Cartn_z']])
                    res_seq = int(parts[indices['_atom_site.auth_seq_id']])
                except ValueError:
                    continue

                def append_data(c, rn, ri, ci, an, al, bf, mid):
                    coords.append(c)
                    res_names.append(rn)
                    res_ids.append(ri)
                    chain_ids.append(ci)
                    atom_names.append(an)
                    alt_locs.append(al)
                    b_factors.append(bf)
                    model_ids.append(mid)

                if self.atom_mode == "centroid":
                    key = (chain_id, res_seq, res_name, model_num)
                    centroid_buffer[key].append([x, y, z])
                else:
                    append_data([x, y, z], res_name, res_seq, chain_id, atom_name, alt, bfactor, model_num)

        if self.atom_mode == "centroid":
            for (chain, seq, rname, mid), atom_coords in centroid_buffer.items():
                coords.append(np.mean(atom_coords, axis=0))
                res_names.append(rname)
                res_ids.append(seq)
                chain_ids.append(chain)
                atom_names.append("CTR")
                alt_locs.append(".")
                b_factors.append(0.0)
                model_ids.append(mid)

        if not coords: return None

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

class DistanceComputer:
    """
    Computes Euclidean distances between atoms using efficient spatial algorithms.
    """
    def __init__(self, cutoff=20.0, seq_separation=4):
        """
        Args:
            cutoff (float): Max distance to consider (default 20A).
            seq_separation (int): Minimum sequence index difference |i - j| for intrachain pairs.
                                  Default 4 enforces the rule "separated by at least 3 positions".
        """
        self.cutoff = cutoff
        self.seq_separation = seq_separation

    def compute(self, data, dist_mode, detailed_output=False):
        """
        Computes pairwise distances using KD-Tree spatial hashing.

        Returns:
            tuple: (distances_dict, detailed_dataframe or None)
        """
        if data is None: return {}, None
        
        coords = data['coords']
        chains = data['chain_ids']
        res_ids = data['res_ids']
        res_names = data['res_names']
        
        atom_names = data['atom_names']
        alt_locs = data['alt_locs']
        b_factors = data['b_factors']
        model_ids = data['model_ids']
        pdb_name = data['pdb_name']
        
        # 1. Spatial Indexing
        tree = cKDTree(coords)
        
        # 2. Neighbor Query (Radius Search)
        pairs = tree.query_pairs(self.cutoff)
        
        if not pairs: return {}, None
        
        pairs_np = np.array(list(pairs))
        idx_i = pairs_np[:, 0]
        idx_j = pairs_np[:, 1]
        
        # 3. Vectorized Filtering
        same_model = (model_ids[idx_i] == model_ids[idx_j])
        same_chain = (chains[idx_i] == chains[idx_j])
        
        if dist_mode == 'intra':
            seq_diff = np.abs(res_ids[idx_i] - res_ids[idx_j])
            mask = same_model & same_chain & (seq_diff >= self.seq_separation)
        else:
            mask = same_model & (~same_chain)

        final_i = idx_i[mask]
        final_j = idx_j[mask]
        
        # 4. Metric Computation
        diffs = coords[final_i] - coords[final_j]
        dists = np.linalg.norm(diffs, axis=1)
        
        types_i = res_names[final_i]
        types_j = res_names[final_j]
        
        results = defaultdict(list)
        pair_keys = []

        for t1, t2, d in zip(types_i, types_j, dists):
            key = t1 + t2 if t1 < t2 else t2 + t1
            pair_keys.append(key)
            results[key].append(d)
            results['XX'].append(d)
            
        detailed_df = None
        if detailed_output and len(dists) > 0:
            detailed_df = pd.DataFrame({
                'PDB': pdb_name,
                'Model': model_ids[final_i],
                'Chain1': chains[final_i],
                'Res1': types_i,
                'ResID1': res_ids[final_i],
                'Atom1': atom_names[final_i],
                'AltLoc1': alt_locs[final_i],
                'B_Factor1': b_factors[final_i],
                'Chain2': chains[final_j],
                'Res2': types_j,
                'ResID2': res_ids[final_j],
                'Atom2': atom_names[final_j],
                'AltLoc2': alt_locs[final_j],
                'B_Factor2': b_factors[final_j],
                'Distance': dists,
                'Type': pair_keys
            })
            
        return results, detailed_df