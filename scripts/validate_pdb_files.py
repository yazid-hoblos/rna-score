#!/usr/bin/env python3
"""
Validate and filter RNA structures before training.
Checks for C3' atoms, chain length, and structure quality.
"""

import sys
from pathlib import Path
from typing import Dict, List, Tuple
import argparse


def parse_pdb_for_validation(pdb_file: Path) -> Dict:
    """
    Parse PDB file and extract validation information.
    
    Returns:
        Dictionary with validation metrics
    """
    info = {
        "filename": pdb_file.name,
        "has_c3": False,
        "chain_count": 0,
        "rna_residues": 0,
        "c3_atoms": 0,
        "chains": {},
        "resolution": None,
        "valid": False,
        "issues": []
    }
    
    chains = {}
    c3_count = 0
    
    try:
        with open(pdb_file, 'r') as f:
            for line in f:
                # Get resolution
                if line.startswith("REMARK   2 RESOLUTION"):
                    try:
                        res_str = line[25:].strip().split()[0]
                        if res_str != "NOT":
                            info["resolution"] = float(res_str)
                    except:
                        pass
                
                # Parse ATOM/HETATM records
                if line.startswith(("ATOM", "HETATM")):
                    atom_name = line[12:16].strip()
                    chain_id = line[21]
                    res_name = line[17:20].strip()
                    
                    # Check if RNA residue
                    if res_name in ["A", "U", "G", "C", "ADE", "URA", "GUA", "CYT", 
                                   "DA", "DU", "DG", "DC", "RA", "RU", "RG", "RC"]:
                        if chain_id not in chains:
                            chains[chain_id] = {"residues": set(), "has_c3": False}
                        
                        res_num = line[22:26].strip()
                        chains[chain_id]["residues"].add(res_num)
                        
                        # Check for C3' atom
                        if atom_name == "C3'" or atom_name == "C3*":
                            chains[chain_id]["has_c3"] = True
                            c3_count += 1
                            info["has_c3"] = True
    
    except Exception as e:
        info["issues"].append(f"Error reading file: {e}")
        return info
    
    # Compile chain information
    for chain_id, chain_data in chains.items():
        residue_count = len(chain_data["residues"])
        info["chains"][chain_id] = {
            "length": residue_count,
            "has_c3": chain_data["has_c3"]
        }
        info["rna_residues"] += residue_count
    
    info["chain_count"] = len(chains)
    info["c3_atoms"] = c3_count
    
    # Validation checks
    if info["chain_count"] == 0:
        info["issues"].append("No RNA chains found")
    if not info["has_c3"]:
        info["issues"].append("No C3' atoms found")
    if info["rna_residues"] < 10:
        info["issues"].append(f"Too few RNA residues ({info['rna_residues']} < 10)")
    
    # Mark as valid if passes all checks
    info["valid"] = len(info["issues"]) == 0
    
    return info


def validate_dataset(input_dir: Path, output_dir: Path = None, 
                     min_residues: int = 10, max_resolution: float = None) -> Tuple[List, List]:
    """
    Validate all PDB files in a directory.
    
    Args:
        input_dir: Directory containing PDB files
        output_dir: Optional directory to copy valid files
        min_residues: Minimum number of RNA residues required
        max_resolution: Maximum resolution (Angstroms)
    
    Returns:
        Tuple of (valid_files, invalid_files)
    """
    pdb_files = list(input_dir.glob("*.pdb"))
    
    if not pdb_files:
        print(f"No PDB files found in {input_dir}")
        return [], []
    
    print(f"Validating {len(pdb_files)} PDB files...")
    print(f"Criteria: min_residues={min_residues}, max_resolution={max_resolution}")
    
    valid_files = []
    invalid_files = []
    
    # Statistics
    stats = {
        "total": len(pdb_files),
        "valid": 0,
        "no_c3": 0,
        "too_small": 0,
        "bad_resolution": 0,
        "no_rna": 0,
        "error": 0
    }
    
    for pdb_file in pdb_files:
        info = parse_pdb_for_validation(pdb_file)
        
        # Additional validation based on parameters
        if info["valid"]:
            if info["rna_residues"] < min_residues:
                info["valid"] = False
                info["issues"].append(f"Too few residues ({info['rna_residues']} < {min_residues})")
                stats["too_small"] += 1
            
            if max_resolution and info["resolution"]:
                if info["resolution"] > max_resolution:
                    info["valid"] = False
                    info["issues"].append(f"Resolution too high ({info['resolution']} > {max_resolution})")
                    stats["bad_resolution"] += 1
        
        # Categorize issues
        if not info["valid"]:
            for issue in info["issues"]:
                if "No C3'" in issue:
                    stats["no_c3"] += 1
                elif "No RNA" in issue:
                    stats["no_rna"] += 1
                elif "Error" in issue:
                    stats["error"] += 1
        
        if info["valid"]:
            valid_files.append((pdb_file, info))
            stats["valid"] += 1
        else:
            invalid_files.append((pdb_file, info))
    
    # Print summary
    print(f"\n{'='*60}")
    print(f"Validation Summary:")
    print(f"  Total files:        {stats['total']:6d}")
    print(f"  Valid:              {stats['valid']:6d} ({100*stats['valid']/stats['total']:.1f}%)")
    print(f"  Invalid:            {stats['total']-stats['valid']:6d} ({100*(stats['total']-stats['valid'])/stats['total']:.1f}%)")
    
    if invalid_files:
        print(f"\nInvalid file reasons:")
        print(f"  No C3' atoms:       {stats['no_c3']:6d}")
        print(f"  Too few residues:   {stats['too_small']:6d}")
        print(f"  Bad resolution:     {stats['bad_resolution']:6d}")
        print(f"  No RNA chains:      {stats['no_rna']:6d}")
        print(f"  Read errors:        {stats['error']:6d}")
    
    # Copy valid files if output directory specified
    if output_dir and valid_files:
        output_dir.mkdir(parents=True, exist_ok=True)
        validated_dir = output_dir / "validated"
        validated_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nCopying {len(valid_files)} valid files to {output_dir}")
        
        import shutil
        for pdb_file, info in valid_files:
            shutil.copy2(pdb_file, validated_dir / pdb_file.name)
    
    # Save validation report
    report_file = (output_dir or input_dir) / "validation_report.txt"
    with open(report_file, 'w') as f:
        f.write("RNA Structure Validation Report\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total files: {stats['total']}\n")
        f.write(f"Valid files: {stats['valid']}\n")
        f.write(f"Invalid files: {stats['total']-stats['valid']}\n\n")
        
        f.write("VALID STRUCTURES:\n")
        f.write("-" * 40 + "\n")
        for pdb_file, info in valid_files:
            f.write(f"{info['filename']:12} | Chains: {info['chain_count']} | ")
            f.write(f"Residues: {info['rna_residues']:4} | C3': {info['c3_atoms']:4}")
            if info['resolution']:
                f.write(f" | Res: {info['resolution']:.2f}Å")
            f.write("\n")
        
        if invalid_files:
            f.write("\n\nINVALID STRUCTURES:\n")
            f.write("-" * 40 + "\n")
            for pdb_file, info in invalid_files[:50]:  # First 50
                f.write(f"{info['filename']:12} | Issues: {', '.join(info['issues'])}\n")
    
    print(f"\nValidation report saved to: {report_file}")
    
    # Save lists of valid/invalid PDB IDs
    valid_ids_file = (output_dir or input_dir) / "valid_pdb_ids.txt"
    with open(valid_ids_file, 'w') as f:
        f.write("\n".join([pdb.stem for pdb, _ in valid_files]))
    print(f"Valid PDB IDs saved to: {valid_ids_file}")
    
    if invalid_files:
        invalid_ids_file = (output_dir or input_dir) / "invalid_pdb_ids.txt"
        with open(invalid_ids_file, 'w') as f:
            f.write("\n".join([pdb.stem for pdb, _ in invalid_files]))
        print(f"Invalid PDB IDs saved to: {invalid_ids_file}")
    
    return valid_files, invalid_files


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Validate RNA structures for folding objective function training"
    )
    
    parser.add_argument(
        "-i", "--input-dir",
        type=Path,
        help="Directory containing PDB files to validate"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default='rna_structures',
        help="Output directory for valid structures (optional)"
    )
    
    parser.add_argument(
        "-m", "--min-residues",
        type=int,
        default=10,
        help="Minimum number of RNA residues (default: 10)"
    )
    
    parser.add_argument(
        "-r", "--max-resolution",
        type=float,
        help="Maximum resolution in Angstroms (optional)"
    )
    
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Show detailed information for each file"
    )
    
    args = parser.parse_args()
    
    if not args.input_dir or not args.input_dir.is_dir():
        print("Error: Please provide a valid input directory containing PDB files.")
        parser.print_help()
        return 1
    
    valid_files, invalid_files = validate_dataset(
        args.input_dir,
        args.output,
        args.min_residues,
        args.max_resolution
    )
    
    if args.verbose and invalid_files:
        print(f"\nDetailed issues for invalid files:")
        for _, info in invalid_files[:20]:  # Show first 20
            print(f"\n{info['filename']}:")
            for issue in info['issues']:
                print(f"  - {issue}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())