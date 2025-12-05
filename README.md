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
rna-score access -n 50 --rna-only -f cif -o rna_structures --workers 4
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

*Each subcommand supports `--help` / `-h` for details.*

---

## Web Interface

You can also use rna-score directly in your browser:  
👉 [https://rna-score.onrender.com/](https://rna-score.onrender.com/)

---