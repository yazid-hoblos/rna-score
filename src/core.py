"""
RNA Structure Parser and Distance Calculator

This module provides tools to fetch, parse, and analyze atomic coordinates 
from PDB and mmCIF files for RNA structures.

Main Components:
----------------
1. OnlineFetcher:
   - Downloads PDB or mmCIF files from the RCSB Protein Data Bank.
   - Returns the file content as a list of lines.

2. FastParser:
   - Efficiently parses PDB or mmCIF files.
   - Supports filtering by chain, atom type, and models.
   - Supports extracting either specific atoms, all atoms, or centroids per residue.
   - Returns structured arrays containing coordinates, residue info, atom info, 
     chain IDs, model IDs, and B-factors.

3. DistanceComputer:
   - Computes pairwise distances between atoms in RNA structures.
   - Supports intra-chain (with optional sequence separation) and inter-chain modes.
   - Can return a dictionary of distances grouped by residue pairs and optionally 
     a detailed pandas DataFrame for all contacts.

Workflow Example:
-----------------
1. Fetch RNA structure:
    lines = OnlineFetcher.get_content('1ABC', file_format='pdb')

2. Parse structure:
    parser = FastParser(atom_mode="C3'")
    data = parser.parse(lines, file_format='pdb', pdb_name='1ABC')

3. Compute distances:
    dist_comp = DistanceComputer(cutoff=20.0, seq_separation=4)
    results, detailed_df = dist_comp.compute(data, dist_mode='intra', detailed_output=True)

Dependencies:
-------------
- numpy
- scipy
- pandas
- urllib
- collections
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
    def __init__(self, atom_mode="C3'"):
        self.atom_mode = atom_mode
        self.allowed_bases = {'A', 'U', 'C', 'G'}
        
        if isinstance(self.atom_mode, list):
            self.atom_mode = set(self.atom_mode)

    def parse(self, content_lines, file_format, pdb_name, chains_filter=None, use_all_models=False):
        """
        Parses the file content depending on its format.

        Parameters:
        -----------
        content_lines : list[str]
            Lines of the PDB or mmCIF file.
        file_format : str
            'pdb' or 'mmcif'
        pdb_name : str
            Identifier for the RNA structure.
        chains_filter : list[str], optional
            List of chain IDs to keep. Default is None (keep all chains).
        use_all_models : bool, optional
            Whether to parse all models in the structure. Default is False.

        Returns:
        --------
        dict or None
            Dictionary containing parsed data arrays, or None if no valid atoms were found.
        """
        if file_format == 'mmcif':
            return self._parse_cif(content_lines, pdb_name, chains_filter, use_all_models)
        else:
            return self._parse_pdb(content_lines, pdb_name, chains_filter, use_all_models)

    def _should_keep_atom(self, atom_name):
        """
        Checks if an atom should be kept based on the atom_mode setting.
        """
        if self.atom_mode == 'all' or self.atom_mode == 'centroid':
            return True
        if isinstance(self.atom_mode, set):
            return atom_name in self.atom_mode
        return atom_name == self.atom_mode

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
        
        centroid_buffer = defaultdict(list)  # Used if atom_mode is 'centroid'

        current_model = 1
        model_found = False
        
        for line in lines:
            # Handle multi-model structures
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

            # Skip alternative locations except 'A' or blank
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

            # Helper function to append parsed data
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

        # Compute centroids if needed
        if self.atom_mode == "centroid":
            for (chain, seq, rname, mid), atom_coords in centroid_buffer.items():
                center = np.mean(atom_coords, axis=0)
                coords.append(center)
                res_names.append(rname)
                res_ids.append(seq)
                chain_ids.append(chain)
                atom_names.append("CTR")
                alt_locs.append(".") 
                b_factors.append(0.0)
                model_ids.append(mid)

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

                # Helper function to append parsed data
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

        # Compute centroids if needed
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
    Computes pairwise distances between atoms in an RNA structure, either intra-chain or inter-chain,
    with optional sequence separation filtering.
    """
    def __init__(self, cutoff=20.0, seq_separation=4):
        """
        Parameters:
        -----------
        cutoff : float
            Maximum distance between atoms to be considered a contact (in Ångströms).
        seq_separation : int
            Minimum sequence separation for intra-chain distance computations.
        """
        self.cutoff = cutoff
        self.seq_separation = seq_separation

    def compute(self, data, dist_mode, detailed_output=False):
        """
        Computes pairwise distances based on the specified mode.

        Parameters:
        -----------
        data : dict
            Parsed RNA structure data from FastParser.
        dist_mode : str
            Either 'intra' (within the same chain) or 'inter' (between different chains).
        detailed_output : bool, optional
            Whether to return a detailed pandas DataFrame with all distances.

        Returns:
        --------
        tuple:
            results : dict
                Dictionary mapping residue pair types to lists of distances.
            detailed_df : pandas.DataFrame or None
                DataFrame with detailed distance information if requested.
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
        
        # Build a KDTree for fast neighbor search
        tree = cKDTree(coords, compact_nodes=True, balanced_tree=True)
        pairs = tree.query_pairs(self.cutoff)
        
        if not pairs: return {}, None
        
        pairs_np = np.array(list(pairs))
        idx_i = pairs_np[:, 0]
        idx_j = pairs_np[:, 1]
        
        same_model = (model_ids[idx_i] == model_ids[idx_j])
        same_chain = (chains[idx_i] == chains[idx_j])
        
        if dist_mode == 'intra':
            seq_diff = np.abs(res_ids[idx_i] - res_ids[idx_j])
            mask = same_model & same_chain & (seq_diff >= self.seq_separation)
        else:
            mask = same_model & (~same_chain)

        final_i = idx_i[mask]
        final_j = idx_j[mask]
        
        diffs = coords[final_i] - coords[final_j]
        dists = np.linalg.norm(diffs, axis=1)
        
        types_i = res_names[final_i]
        types_j = res_names[final_j]
        
        results = defaultdict(list)
        pair_keys = []

        # Organize distances by residue type
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
