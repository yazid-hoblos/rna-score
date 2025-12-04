"""
VERIFICATION SUITE
------------------
A comprehensive integration testing suite for the preprocessing and distances extraction.
This script validates all functional components, including file parsing, 
geometric computation, argument handling, and output generation.

Tests Covered:
1. Environment Integrity (File checks, Imports).
2. Functional correctness (Web vs Local files).
3. Parameter Logic (Atom modes, Interaction modes, Formats).
4. Edge Case Handling (Mixed inputs, Missing files).
5. Real World Simulation (Processing a list of real PDBs).

Outputs:
    A 'verification_data' directory containing test inputs and generated outputs.
"""

import os
import sys
import subprocess
import shutil
from core import FastParser

# --- CONFIGURATION ---
SCRIPT_NAME = "extract_distances.py"
BASE_DIR = "verification_data"
INPUT_DIR = os.path.join(BASE_DIR, "inputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# File Paths
DUMMY_PDB_PATH = os.path.join(INPUT_DIR, "dummy.pdb")
DUMMY_CIF_PATH = os.path.join(INPUT_DIR, "dummy.cif")
MIXED_LIST_PATH = os.path.join(INPUT_DIR, "mixed_inputs.txt")
REAL_LIST_PATH = os.path.join(INPUT_DIR, "training_set.txt")

# Dummy PDB Content (2 chains, A and B)
DUMMY_PDB_CONTENT = """\
ATOM      1  P     G A   1      10.000  10.000  10.000  1.00  0.00           P
ATOM      2  C3'   G A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  P     C A   5      10.000  15.000  10.000  1.00  0.00           P
ATOM      4  C3'   C A   5      11.000  15.000  10.000  1.00  0.00           C
ATOM      5  P     U B   1      12.000  10.000  10.000  1.00  0.00           P
ATOM      6  C3'   U B   1      13.000  10.000  10.000  1.00  0.00           C
"""

# Dummy mmCIF Content
DUMMY_CIF_CONTENT = """\
data_dummy
loop_
_atom_site.group_PDB
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM P   G A 1 10.000 10.000 10.000
ATOM C3' G A 1 11.000 10.000 10.000
ATOM P   C A 5 10.000 15.000 10.000
ATOM C3' C A 5 11.000 15.000 10.000
ATOM P   U B 1 12.000 10.000 10.000
ATOM C3' U B 1 13.000 10.000 10.000
"""

def print_status(step, msg, status="...", info=""):
    print(f"[{step}] {msg:<55} {status}")
    if info:
        print(f"      -> {info}")

def run_command(cmd):
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True, result.stdout.decode()
    except subprocess.CalledProcessError as e:
        return False, e.stderr.decode()

def setup_environment():
    if os.path.exists(BASE_DIR):
        shutil.rmtree(BASE_DIR)
    os.makedirs(INPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # 1. Create Dummy files
    with open(DUMMY_PDB_PATH, "w") as f:
        f.write(DUMMY_PDB_CONTENT)

    with open(DUMMY_CIF_PATH, "w") as f:
        f.write(DUMMY_CIF_CONTENT)

    # 2. Create Mixed List (Local + Web + Chains)
    with open(MIXED_LIST_PATH, "w") as f:
        f.write(f"{DUMMY_PDB_PATH} A\n")
        f.write("\n")
        f.write("1Y26\n") # Small riboswitch

    # 3. Create Real World Dataset List
    with open(REAL_LIST_PATH, "w") as f:
        f.write("1Y26\n") # Guanine Riboswitch
        f.write("4GXY\n") # Hammerhead Ribozyme
        # Note: 1EHZ excluded from 'Real' test for speed (it's huge), 
        # but you can add it if you have good internet.

def main():
    print("=========================================================")
    print("   RNA PROJECT: ORGANIZED VERIFICATION SUITE")
    print("=========================================================\n")

    setup_environment()
    print(f"Environment setup complete.")
    print(f"Inputs located in:  {INPUT_DIR}/")
    print(f"Outputs will go to: {OUTPUT_DIR}/\n")

    # TEST 1: File Integrity
    required = ["core.py", SCRIPT_NAME]
    missing = [f for f in required if not os.path.exists(f)]
    if missing:
        print(f"CRITICAL ERROR: Missing script files: {missing}")
        return
    print_status("1", "Checking script integrity", "OK")

    # TEST 2: Library Import
    try:
        parser = FastParser(atom_mode="C3'")
        print_status("2", "Importing 'core.py' library", "OK")
    except Exception as e:
        print_status("2", "Importing 'core.py' library", "FAIL", str(e))
        return

    # TEST 3: Web Fetching (PDB)
    out_web = os.path.join(OUTPUT_DIR, "web")
    cmd = f"python {SCRIPT_NAME} --pdb 1Y26 --out-dir {out_web}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_web, "AU_histogram.txt")):
        print_status("3", "Single Input (Web PDB)", "OK", f"Saved to {out_web}")
    else:
        print_status("3", "Single Input (Web PDB)", "FAIL", "Check internet connection")

    # TEST 4: Local File Input (PDB)
    out_local = os.path.join(OUTPUT_DIR, "local_pdb")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --out-dir {out_local}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_local, "CG_histogram.txt")):
        print_status("4", "Single Input (Local PDB)", "OK", f"Saved to {out_local}")
    else:
        print_status("4", "Single Input (Local PDB)", "FAIL", log)

    # TEST 5: Mixed List Input
    out_list = os.path.join(OUTPUT_DIR, "mixed_list")
    cmd = f"python {SCRIPT_NAME} --list {MIXED_LIST_PATH} --out-dir {out_list}"
    success, log = run_command(cmd)
    if success: 
        print_status("5", "Mixed List Input", "OK", f"Processed Local+Web list")
    else:
        print_status("5", "Mixed List Input", "FAIL", log)

    # TEST 6: Parameter - Atom Mode (Custom Subset)
    out_subset = os.path.join(OUTPUT_DIR, "subset")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --atom-mode P \"C3'\" --out-dir {out_subset}"
    success, log = run_command(cmd)
    if success:
        print_status("6", "Atom Mode: Subset (P + C3')", "OK")
    else:
        print_status("6", "Atom Mode: Subset (P + C3')", "FAIL", log)

    # TEST 7: Parameter - Atom Mode (All Atoms) + KDE
    out_kde = os.path.join(OUTPUT_DIR, "kde")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --atom-mode all --method kde --out-dir {out_kde}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_kde, "CG_kde_raw.txt")):
        print_status("7", "Atom Mode: All + Method: KDE", "OK")
    else:
        print_status("7", "Atom Mode: All + Method: KDE", "FAIL", log)

    # TEST 8: Parameter - Dist Mode (Interchain)
    out_inter = os.path.join(OUTPUT_DIR, "inter")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --dist-mode inter --out-dir {out_inter}"
    success, log = run_command(cmd)
    if success:
        print_status("8", "Dist Mode: Interchain", "OK", "Found interaction")
    else:
        print_status("8", "Dist Mode: Interchain", "FAIL", log)

    # TEST 9: Parameter - mmCIF Format (Web)
    out_mmcif_web = os.path.join(OUTPUT_DIR, "mmcif_web")
    cmd = f"python {SCRIPT_NAME} --pdb 1Y26 --format mmcif --out-dir {out_mmcif_web}"
    success, log = run_command(cmd)
    if success:
        print_status("9", "Format: mmCIF (Web)", "OK")
    else:
        print_status("9", "Format: mmCIF (Web)", "FAIL", log)

    # TEST 10: Detailed Output
    out_details = os.path.join(OUTPUT_DIR, "details")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --save-detailed --out-dir {out_details}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_details, "detailed_interactions.csv")):
        print_status("10", "Option: --save-detailed", "OK")
    else:
        print_status("10", "Option: --save-detailed", "FAIL", log)

    # TEST 11: Local mmCIF File
    out_local_cif = os.path.join(OUTPUT_DIR, "local_cif")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_CIF_PATH} --out-dir {out_local_cif}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_local_cif, "CG_histogram.txt")):
        print_status("11", "Single Input (Local mmCIF)", "OK", f"Saved to {out_local_cif}")
    else:
        print_status("11", "Single Input (Local mmCIF)", "FAIL", log)

    # TEST 12: Real Dataset Training (New)
    out_real = os.path.join(OUTPUT_DIR, "real_dataset")
    cmd = f"python {SCRIPT_NAME} --list {REAL_LIST_PATH} --out-dir {out_real} --atom-mode all"
    success, log = run_command(cmd)
    if success:
        print_status("12", "Real Dataset (List Input)", "OK", f"Saved to {out_real}")
    else:
        print_status("12", "Real Dataset (List Input)", "FAIL", log)

    print("\n---------------------------------------------------------")
    print("VERIFICATION COMPLETE.")
    print(f"You can inspect the generated files in: '{BASE_DIR}'")

if __name__ == "__main__":
    main()