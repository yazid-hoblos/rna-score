# rna-score: RNA Scoring Library

**Goal:** Creation of an objective function for the RNA folding problem

This project develops a scoring function to evaluate predicted RNA tertiary structures based on interatomic distance distributions.

**Supervised by:** Professor Guillaume Postic  
**Team:** Rayane Adams, Joelle Assy, Yazid Hoblos, Denys Buryi, Raul Duran De Alba

---

## Overview

rna-score provides tools to download, process, and score RNA 3D structures using statistical models of interatomic distances. The library is available as a Python CLI tool and is also accessible via a web interface:

🌐 **Try it online:** [https://rna-score.onrender.com/](https://rna-score.onrender.com/)

---

## Features Supported for Distances Computation

| Component           | Description                                                         | CLI Option                                   | Details                                                                                                                                                                           |
|---------------------|---------------------------------------------------------------------|-----------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| **Input Source**    | Select the structure(s) to process, either remote or local.         | `--pdb`, `--list`, `--folder`                 | Choose **one**: a PDB ID, a local file, a list file (`<ID> [CHAIN ...]`), or a directory of structures.                                                                          |
| **Input Format**    | Specify the format used for parsing the structure.                  | `--format {pdb, mmcif}`                       | Default: **pdb**. Automatically detected for local files.                                                                                                                         |
| **Atomic Selection**| Choose how the structure is represented for distance calculations.| `--atom-mode`                                 | Options: **"C3'"** (default), **centroid**, **all**, or multiple atom names (e.g., `"P" "C4"`).                                                                                     |
| **Interaction Mode**| Determine whether distances are measured within or between chains.  | `--dist-mode {intra, inter}`                  | **intra** (default): within one chain. **inter**: between distinct chains.                                                                                                       |
| **Sequence Separation** | Minimum offset for intra-chain contacts.                        | `--seq-sep SEQ_SEP`                           | Default: **4** residues. Ignored in *inter* mode. Distances considered from *i* to *i+4*.                                                                                         |
| **Distance Cutoff** | Maximum atom–atom distance (Å) counted as a contact.                | `--cutoff CUTOFF`                             | Default: **20.0 Å**.                                                                                                                                                              |
| **Output Type**     | Determines the type of distance distribution produced.              | `--method {histogram, kde}`                   | **histogram** (default): binned counts. **kde**: raw distances for kernel density estimation.                                                                                     |
| **Parallelization** | Control how many CPU cores to use.                                  | `--cores CORES`                               | Default: all available cores.                                                                                                                                                     |
| **NMR Models**      | Whether to process all models in NMR structures.                    | `--all-models`                                | Default: only the first model is used.                                                                                                                                            |
| **Detailed Log**    | Save a CSV file with full information for every measured distance.  | `--save-detailed`                             | Saves: **PDB**, **Model**, **Chain IDs**, **Residue IDs**, **Atom Names**, **B-factors**, **AltLocs**, **Distance**, and **Pair Type**.                                          |
| **Output Directory**| Location where results will be written.                             | `--out-dir OUT_DIR`                           | Default: `dist_data/`.                                                                                                                                                            |

## Code Structure

- **`src/`**  
	Main Python package containing:
	- `rna_score/cli.py` – Command-line interface entry point
	- `access_rna_structures.py` – Downloading RNA structures
	- `extract_distances.py` – Distance extraction from structures
	- `kde_training.py`, `train.py` – Training scoring tables (histogram/KDE)
	- `score_structures.py` – Scoring new structures
	- `plot_distributions.py`, `plot_scores.py` – Visualization utilities
	- `utils/` – Helper functions (e.g., structure I/O, validation)

- **`tests/`**  
	Unit and integration tests

- **`requirements.txt`**  
	Python dependencies

- **`setup.py`**  
	Installation script

---

## CLI Usage

Install (editable mode for development):

```bash
pip install -r requirements.txt
pip install -e .
```

### 1. Download RNA Structures

```bash
rna-score access -n 50 --rna-only -f cif -o data/rna_structures --workers 4
```
*Add `--validate` to filter out invalid downloaded files.*

### 2. Extract Interatomic Distances

```bash
rna-score extract --folder rna_structures/mmcif --format mmcif --out-dir dist_data
```

### 3. Train Scoring Tables

```bash
rna-score train --input-dir dist_data --output-dir training_output --method histogram
```

### 4. Score Structures

```bash
rna-score score --folder rna_structures/mmcif --tables training_output --format mmcif --output scores.csv
```

### 5. Plot Results

```bash
rna-score plot --input-dir training_output --output-dir plots --combined
```

### 6. Full Workflow (all steps in one command)

```bash
# add pdb ids and chains for scoring
cat <<EOF > tests/scoring_list.txt
1EHZ A
1Y26 B C
EOF

rna-score workflow --train-folder data/rna_structures/mmcif --score-list tests/scoring_list.txt --output-dir tests/workflow_output --format mmcif --method histogram
```

This runs extraction, training, scoring, and plotting in a single step. See `rna-score workflow --help` for all options.

*Each subcommand supports `--help` / `-h` for details.*

## Web Interface

You can also use rna-score directly in your browser:  
👉 [https://rna-score.onrender.com/](https://rna-score.onrender.com/)

---