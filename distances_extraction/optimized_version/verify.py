"""
RNA Distance Extraction Verification Suite

This script sets up a controlled environment to test and verify the functionality
of the RNA distance extraction pipeline (`extract_distances.py`). It creates dummy
PDB and mmCIF files, a mixed input list, and runs multiple test cases including:
- Local and web PDB files
- Single and multiple structures
- Various atom modes
- Distance computation modes (intra- and inter-chain)
- Histogram and raw distance outputs (kde)
- Detailed CSV logs

It reports success/failure for each test step and organizes outputs in a dedicated directory.

Structure
---------
BASE_DIR/
    inputs/
        dummy.pdb      # dummy RNA structure in PDB format
        dummy.cif      # dummy RNA structure in mmCIF format
        mixed_inputs.txt  # list of local and web structures
    outputs/
        web/
        local_pdb/
        mixed_list/
        subset/
        kde/
        inter/
        mmcif_web/
        details/
        local_cif/
"""

import os
import sys
import subprocess
import shutil
from core import FastParser

# Script and directory configuration
SCRIPT_NAME = "extract_distances.py"
BASE_DIR = "verification_data"
INPUT_DIR = os.path.join(BASE_DIR, "inputs")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

DUMMY_PDB_PATH = os.path.join(INPUT_DIR, "dummy.pdb")
DUMMY_CIF_PATH = os.path.join(INPUT_DIR, "dummy.cif")
MIXED_LIST_PATH = os.path.join(INPUT_DIR, "mixed_inputs.txt")

# Dummy RNA PDB content
DUMMY_PDB_CONTENT = """\
ATOM      1  P     G A   1      10.000  10.000  10.000  1.00  0.00           P
ATOM      2  C3'   G A   1      11.000  10.000  10.000  1.00  0.00           C
ATOM      3  P     C A   5      10.000  15.000  10.000  1.00  0.00           P
ATOM      4  C3'   C A   5      11.000  15.000  10.000  1.00  0.00           C
ATOM      5  P     U B   1      12.000  10.000  10.000  1.00  0.00           P
ATOM      6  C3'   U B   1      13.000  10.000  10.000  1.00  0.00           C
"""

# Dummy RNA mmCIF content
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
    """
    Print formatted status messages for verification steps.

    Parameters
    ----------
    step : str
        Step number or identifier.
    msg : str
        Description of the step.
    status : str, optional
        Status indicator (default: "...").
    info : str, optional
        Additional information to display (default: "").
    """
    print(f"[{step}] {msg:<55} {status}")
    if info:
        print(f"      -> {info}")

def run_command(cmd):
    """
    Run a shell command and capture output or errors.

    Parameters
    ----------
    cmd : str
        Command to execute.

    Returns
    -------
    tuple
        - success : bool, True if command ran successfully
        - output : str, captured stdout or stderr
    """
    try:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return True, result.stdout.decode()
    except subprocess.CalledProcessError as e:
        return False, e.stderr.decode()

def setup_environment():
    """
    Create input/output directories and populate dummy RNA files for testing.
    """
    if os.path.exists(BASE_DIR):
        shutil.rmtree(BASE_DIR)
    os.makedirs(INPUT_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    with open(DUMMY_PDB_PATH, "w") as f:
        f.write(DUMMY_PDB_CONTENT)

    with open(DUMMY_CIF_PATH, "w") as f:
        f.write(DUMMY_CIF_CONTENT)

    with open(MIXED_LIST_PATH, "w") as f:
        f.write(f"{DUMMY_PDB_PATH} A\n")
        f.write("\n")
        f.write("1EHZ A\n")
        f.write("1Y26\n")

def main():
    """
    Run the verification suite, testing multiple configurations of the RNA distance pipeline.

    Test steps include:
    1. Script integrity
    2. Library import
    3. Web PDB input
    4. Local PDB input
    5. Mixed input list
    6. Atom-mode subset selection
    7. Atom-mode "all" + method KDE
    8. Inter-chain distances
    9. Web mmCIF input
    10. Save detailed interactions
    11. Local mmCIF input

    Outputs are organized under `verification_data/outputs/`.
    """
    print("=========================================================")
    print("   RNA PROJECT: VERIFICATION SUITE")
    print("=========================================================\n")

    setup_environment()
    print(f"Environment setup complete.")
    print(f"Inputs located in:  {INPUT_DIR}/")
    print(f"Outputs will go to: {OUTPUT_DIR}/\n")

    # Check for required files
    required = ["core.py", SCRIPT_NAME]
    missing = [f for f in required if not os.path.exists(f)]
    if missing:
        print(f"CRITICAL ERROR: Missing script files: {missing}")
        return
    print_status("1", "Checking script integrity", "OK")

    # Check library import
    try:
        parser = FastParser(atom_mode="C3'")
        print_status("2", "Importing 'core.py' library", "OK")
    except Exception as e:
        print_status("2", "Importing 'core.py' library", "FAIL", str(e))
        return

    # Run multiple test cases
    out_web = os.path.join(OUTPUT_DIR, "web")
    cmd = f"python {SCRIPT_NAME} --pdb 1EHZ --chains A --out-dir {out_web}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_web, "AU_histogram.txt")):
        print_status("3", "Single Input (Web PDB)", "OK", f"Saved to {out_web}")
    else:
        print_status("3", "Single Input (Web PDB)", "FAIL", "Check internet connection")

    out_local = os.path.join(OUTPUT_DIR, "local_pdb")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --out-dir {out_local}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_local, "CG_histogram.txt")):
        print_status("4", "Single Input (Local PDB)", "OK", f"Saved to {out_local}")
    else:
        print_status("4", "Single Input (Local PDB)", "FAIL", log)

    out_list = os.path.join(OUTPUT_DIR, "mixed_list")
    cmd = f"python {SCRIPT_NAME} --list {MIXED_LIST_PATH} --out-dir {out_list}"
    success, log = run_command(cmd)
    if success:
        print_status("5", "Mixed List Input", "OK", f"Processed Local+Web list")
    else:
        print_status("5", "Mixed List Input", "FAIL", log)

    out_subset = os.path.join(OUTPUT_DIR, "subset")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --atom-mode P \"C3'\" --out-dir {out_subset}"
    success, log = run_command(cmd)
    if success:
        print_status("6", "Atom Mode: Subset (P + C3')", "OK", f"Saved to {out_subset}")
    else:
        print_status("6", "Atom Mode: Subset (P + C3')", "FAIL", log)

    out_kde = os.path.join(OUTPUT_DIR, "kde")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --atom-mode all --method kde --out-dir {out_kde}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_kde, "CG_kde_raw.txt")):
        print_status("7", "Atom Mode: All + Method: KDE", "OK")
    else:
        print_status("7", "Atom Mode: All + Method: KDE", "FAIL", log)

    out_inter = os.path.join(OUTPUT_DIR, "inter")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --dist-mode inter --out-dir {out_inter}"
    success, log = run_command(cmd)
    if success:
        print_status("8", "Dist Mode: Interchain", "OK", "Found interaction")
    else:
        print_status("8", "Dist Mode: Interchain", "FAIL", log)

    out_mmcif_web = os.path.join(OUTPUT_DIR, "mmcif_web")
    cmd = f"python {SCRIPT_NAME} --pdb 1EHZ --format mmcif --chains A --out-dir {out_mmcif_web}"
    success, log = run_command(cmd)
    if success:
        print_status("9", "Format: mmCIF (Web)", "OK")
    else:
        print_status("9", "Format: mmCIF (Web)", "FAIL", log)

    out_details = os.path.join(OUTPUT_DIR, "details")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_PDB_PATH} --save-detailed --out-dir {out_details}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_details, "detailed_interactions.csv")):
        print_status("10", "Option: --save-detailed", "OK")
    else:
        print_status("10", "Option: --save-detailed", "FAIL", log)

    out_local_cif = os.path.join(OUTPUT_DIR, "local_cif")
    cmd = f"python {SCRIPT_NAME} --pdb {DUMMY_CIF_PATH} --out-dir {out_local_cif}"
    success, log = run_command(cmd)
    if success and os.path.exists(os.path.join(out_local_cif, "CG_histogram.txt")):
        print_status("11", "Single Input (Local mmCIF)", "OK", f"Saved to {out_local_cif}")
    else:
        print_status("11", "Single Input (Local mmCIF)", "FAIL", log)

    print("\n---------------------------------------------------------")
    print("VERIFICATION COMPLETE.")
    print(f"You can inspect the generated files in: '{BASE_DIR}'")

if __name__ == "__main__":
    main()
