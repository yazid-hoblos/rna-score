# structural-RNA-project - M2 GENIOMHE 

Goal: Creation of an objective function for the RNA folding problem

## Installation

```bash
pip install -r requirements.txt
```

## RNA Structures Setup

To download RNA structures from PDB in both PDB and mmCIF formats

```bash
python3 scripts/access_rna_structures.py -n 50 --rna-only
```
The script uses the JSON queries in the `queries/` directory.

### Options Summary

- `-n, --number`: Maximum number of structures to download (default: 10)
- `--all`: Download all available RNA structures (overrides -n)
- `--rna-only`: Download only pure RNA structures (no protein/DNA); sorted by resolution (best first)
- `-f, --formats`: File formats to download: pdb, cif (default: both)
- `-o, --output`: Output directory (default: rna_structures)
- `-w, --workers`: Number of parallel downloads (default: 5)
- `--info`: Show information about structures
- `--list-only`: Only list PDB IDs without downloading

### Output Structure

```
rna_structures/
├── pdb/           # PDB format files
│   ├── 1a1t.pdb
│   ├── 1a34.pdb
│   └── ...
├── mmcif/         # mmCIF format files
│   ├── 1a1t.cif
│   ├── 1a34.cif
│   └── ...
├── downloaded_ids.txt    # List of downloaded PDB IDs
└── failed_downloads.txt  # List of failed downloads (if any)
```

## PDB Files Validation

To validate the downloaded PDB files:

```bash
python3 scripts/validate_pdb_files.py -i rna_structures/pdb/
```

This generates (in `rna_structures/` by default):
- `validation_report.txt`: Summary of validation results
- `valid_pdbs_ids.txt`: List of valid PDB file IDs
- `invalid_pdbs_ids.txt`: List of invalid PDB file IDs
- `validated/`: Directory containing only the valid PDB files