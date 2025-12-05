"""
Command-line interface for RNA scoring library.
Orchestrates existing scripts for training, scoring, and plotting.
"""
import argparse
import sys
import os
import subprocess


def main():
    parser = argparse.ArgumentParser(description="RNA scoring pipeline")
    subparsers = parser.add_subparsers(dest='command', required=True)

    # Extract distances from structure IDs
    extract_parser = subparsers.add_parser('extract', help='Extract distances from structure IDs')
    extract_parser.add_argument('--list', type=str, help='Text file containing PDB IDs and optional chain IDs. Format: <PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]. Example: 1EHZ A, 1Y26 B C, 2OEU')
    extract_parser.add_argument('--folder', type=str, help='Directory with PDB/mmCIF files')
    extract_parser.add_argument('--pdb', type=str, help='Single PDB/mmCIF file')
    extract_parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help='Structure format')
    extract_parser.add_argument('--atom-mode', nargs='+', default=["C3'"], help='Atom selection mode(s)')
    extract_parser.add_argument('--chains', nargs='+', help='Specific chains to process')
    extract_parser.add_argument('--dist-mode', choices=['intra', 'inter'], default='intra')
    extract_parser.add_argument('--cutoff', type=float, default=20.0)
    extract_parser.add_argument('--seq-sep', type=int, default=4)
    extract_parser.add_argument('--bin-width', type=float, default=1.0)
    extract_parser.add_argument('--cores', type=int, help='Number of cores to use')
    extract_parser.add_argument('--out-dir', type=str, default='dist_data')
    extract_parser.add_argument('--save-detailed', action='store_true', help='Export detailed CSV log')
    extract_parser.add_argument('--all-models', action='store_true', help='Process all NMR models')
    extract_parser.add_argument('--method', choices=['histogram', 'kde'], default='histogram', help='Extraction method')

    # Train scoring function
    train_parser = subparsers.add_parser('train', help='Train scoring function from distances')
    train_parser.add_argument('--input-dir', type=str, required=True, help='Directory with histogram files')
    train_parser.add_argument('--output-dir', type=str, default='training_output')
    train_parser.add_argument('--max-score', type=float, default=10.0)
    train_parser.add_argument('--pseudocount', type=float, default=1e-6)
    train_parser.add_argument('--cutoff', type=float, default=20.0)
    train_parser.add_argument('--bin-width', type=float, default=1.0)
    train_parser.add_argument('--method', choices=['histogram', 'kde'], default='histogram', help='Training method: histogram (default) or kde')

    # Score structures
    score_parser = subparsers.add_parser('score', help='Score RNA structures')
    score_parser.add_argument('--list', type=str, help='Text file containing structure file paths (local files or PDB IDs). Example: /path/to/file1.pdb, /path/to/file2.cif, 1EHZ')
    score_parser.add_argument('--pdb', type=str, help='Single PDB/mmCIF file to score')
    score_parser.add_argument('--folder', type=str, help='Directory with PDB/mmCIF files')
    score_parser.add_argument('--tables', type=str, default='training_output', help='Directory with score tables')
    score_parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help='Input file format')
    score_parser.add_argument('--cutoff', type=float, default=20.0)
    score_parser.add_argument('--seq-sep', type=int, default=4)
    score_parser.add_argument('--detailed', action='store_true', help='Print detailed per-interaction scores')
    score_parser.add_argument('--output', type=str, help='Output file for scores (CSV)')

    # Plot scores
    plot_parser = subparsers.add_parser('plot', help='Plot score distributions')
    plot_parser.add_argument('--input-dir', type=str, default='training_output', help='Directory with score tables')
    plot_parser.add_argument('--output-dir', type=str, default='plots', help='Directory for output plots')
    plot_parser.add_argument('--combined', action='store_true', help='Generate combined plot with all profiles')

    # Full workflow (extract → train → score → plot)
    workflow_parser = subparsers.add_parser('workflow', help='Run complete pipeline: extract → train → score → plot')
    workflow_parser.add_argument('--train-list', type=str, help='Text file with PDB IDs for training. Format: <PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]. Example: 1EHZ A, 1Y26 B C, 2OEU')
    workflow_parser.add_argument('--train-folder', type=str, help='Directory with PDB/mmCIF files for training')
    workflow_parser.add_argument('--score-list', type=str, help='Text file with PDB IDs for scoring. Format: <PDB_ID> [CHAIN_ID1 CHAIN_ID2 ...]. Example: 1EHZ A, 1Y26 B C')
    workflow_parser.add_argument('--score-folder', type=str, help='Directory with PDB/mmCIF files for scoring')
    workflow_parser.add_argument('--atom-mode', nargs='+', default=["C3'"], help='Atom selection mode(s)')
    workflow_parser.add_argument('--chains', nargs='+', help='Specific chains to process')
    workflow_parser.add_argument('--dist-mode', choices=['intra', 'inter'], default='intra')
    workflow_parser.add_argument('--cutoff', type=float, default=20.0)
    workflow_parser.add_argument('--seq-sep', type=int, default=4)
    workflow_parser.add_argument('--bin-width', type=float, default=1.0)
    workflow_parser.add_argument('--method', choices=['histogram', 'kde'], default='histogram', help='Extraction/training method')
    workflow_parser.add_argument('--max-score', type=float, default=10.0, help='Maximum score for training')
    workflow_parser.add_argument('--pseudocount', type=float, default=1e-6)
    workflow_parser.add_argument('--output-dir', type=str, default='workflow_output', help='Base output directory')
    workflow_parser.add_argument('--no-plot', action='store_true', help='Skip plot generation')
    workflow_parser.add_argument('--cores', type=int, help='Number of cores for extraction')
    workflow_parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help='Structure format')
    workflow_parser.add_argument('--save-detailed', action='store_true', help='Save detailed extraction logs')

    # Access structures (download/organize)
    access_parser = subparsers.add_parser('access', help='Download RNA structures from PDB')
    access_parser.add_argument('-n', '--number', type=int, default=10,
                              help='Maximum number of structures to download (default: 10)')
    access_parser.add_argument('--all', action='store_true',
                              help='Download all available RNA structures (overrides -n)')
    access_parser.add_argument('--rna-only', action='store_true',
                              help='Download only pure RNA structures (no protein/DNA)')
    access_parser.add_argument('-f', '--formats', nargs='+', choices=['pdb', 'cif'], 
                              default=['pdb', 'cif'],
                              help='File formats to download (default: both)')
    access_parser.add_argument('-o', '--output', type=str, default='rna_structures',
                              help='Output directory (default: rna_structures)')
    access_parser.add_argument('-w', '--workers', type=int, default=5,
                              help='Number of parallel downloads (default: 5)')
    access_parser.add_argument('--info', action='store_true',
                              help='Show information about structures before downloading')
    access_parser.add_argument('--list-only', action='store_true',
                              help='Only list PDB IDs without downloading')
    access_parser.add_argument('--validate', action='store_true',
                              help='Run validation on downloaded structures')

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        print("""
    ====================================
         Welcome to rna-score CLI
    ====================================
Analyze, train, and score RNA structures with ease.
    """)
        sys.exit(0)
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(script_dir)
    
    
    try:
        if args.command == 'workflow':
            import shutil
            
            # Validate that training and scoring inputs are provided
            if not (args.train_list or args.train_folder):
                raise ValueError("Must provide --train-list or --train-folder for training data")
            if not (args.score_list or args.score_folder):
                raise ValueError("Must provide --score-list or --score-folder for scoring data")
            
            # Create base output directory
            base_dir = args.output_dir
            os.makedirs(base_dir, exist_ok=True)
            
            extract_dir = os.path.join(base_dir, 'extracted')
            training_dir = os.path.join(base_dir, 'training_output')
            scores_file = os.path.join(base_dir, 'scores.csv')
            plots_dir = os.path.join(base_dir, 'plots')
            
            print("=" * 60)
            print("RNA SCORE WORKFLOW")
            print("=" * 60)
            
            # Step 1: Extract distances from training data
            print("\n[Step 1/4] Extracting distances from training structures...")
            extract_cmd = [
                sys.executable, os.path.join(src_dir, 'extract_distances.py'),
                '--format', args.format,
                '--atom-mode'] + args.atom_mode + [
                '--dist-mode', args.dist_mode,
                '--cutoff', str(args.cutoff),
                '--seq-sep', str(args.seq_sep),
                '--bin-width', str(args.bin_width),
                '--method', args.method,
                '--out-dir', extract_dir
            ]
            if args.train_list:
                extract_cmd.extend(['--list', args.train_list])
            elif args.train_folder:
                extract_cmd.extend(['--folder', args.train_folder])
            if args.chains:
                extract_cmd.extend(['--chains'] + args.chains)
            if args.save_detailed:
                extract_cmd.append('--save-detailed')
            if args.cores:
                extract_cmd.extend(['--cores', str(args.cores)])
            
            subprocess.run(extract_cmd, check=True)
            print(f"✓ Distances extracted to {extract_dir}")
            
            # Step 2: Train scoring function
            print("\n[Step 2/4] Training scoring function...")
            train_cmd = [
                sys.executable, os.path.join(src_dir, 'train.py'),
                '--input-dir', extract_dir,
                '--output-dir', training_dir,
                '--max-score', str(args.max_score),
                '--pseudocount', str(args.pseudocount),
                '--cutoff', str(args.cutoff),
                '--bin-width', str(args.bin_width),
                '--method', args.method
            ]
            subprocess.run(train_cmd, check=True)
            print(f"✓ Training complete. Output in {training_dir}")
            
            # Step 3: Score test structures
            print("\n[Step 3/4] Scoring test structures...")
            score_cmd = [
                sys.executable, os.path.join(src_dir, 'score_structures.py'),
                '--format', args.format,
                '--tables', training_dir,
                '--cutoff', str(args.cutoff),
                '--seq-sep', str(args.seq_sep),
                '--output', scores_file
            ]
            
            # Score the structures from the scoring dataset
            if args.score_list:
                score_cmd.extend(['--list', args.score_list])
            elif args.score_folder:
                score_cmd.extend(['--folder', args.score_folder])
            
            subprocess.run(score_cmd, check=True)
            print(f"✓ Scoring complete. Results in {scores_file}")
            
            # Step 4: Generate plots
            if not args.no_plot:
                print("\n[Step 4/4] Generating plots...")
                plot_cmd = [
                    sys.executable, os.path.join(src_dir, 'plot_scores.py'),
                    '--input-dir', training_dir,
                    '--output-dir', plots_dir
                ]
                subprocess.run(plot_cmd, check=True)
                print(f"✓ Plots generated in {plots_dir}")
            else:
                print("\n[Step 4/4] Skipping plot generation")
            
            print("\n" + "=" * 60)
            print(f"✓ WORKFLOW COMPLETE")
            print("=" * 60)
            print(f"Output directory: {os.path.abspath(base_dir)}")
            print(f"  - Extracted distances: {extract_dir}")
            print(f"  - Training output: {training_dir}")
            print(f"  - Scores file: {scores_file}")
            if not args.no_plot:
                print(f"  - Plots: {plots_dir}")
            print("=" * 60)
        
        elif args.command == 'extract':
            cmd = [
                sys.executable, os.path.join(src_dir, 'extract_distances.py'),
                '--format', args.format,
                '--atom-mode'] + args.atom_mode + [
                '--dist-mode', args.dist_mode,
                '--cutoff', str(args.cutoff),
                '--seq-sep', str(args.seq_sep),
                '--bin-width', str(args.bin_width),
                '--method', args.method,
                '--out-dir', args.out_dir
            ]
            if args.list:
                cmd.extend(['--list', args.list])
            elif args.folder:
                cmd.extend(['--folder', args.folder])
            elif args.pdb:
                cmd.extend(['--pdb', args.pdb])
            else:
                raise ValueError("Must provide --list, --folder, or --pdb")
            if args.chains:
                cmd.extend(['--chains'] + args.chains)
            if args.all_models:
                cmd.append('--all-models')
            if args.save_detailed:
                cmd.append('--save-detailed')
            if args.cores:
                cmd.extend(['--cores', str(args.cores)])
            subprocess.run(cmd, check=True)
            print(f"✓ Distances extracted to {args.out_dir}")
        
        elif args.command == 'train':
            cmd = [
                sys.executable, os.path.join(src_dir, 'train.py'),
                '--input-dir', args.input_dir,
                '--output-dir', args.output_dir,
                '--max-score', str(args.max_score),
                '--pseudocount', str(args.pseudocount),
                '--cutoff', str(args.cutoff),
                '--bin-width', str(args.bin_width),
                '--method', args.method
            ]
            subprocess.run(cmd, check=True)
            print(f"✓ Training complete. Output in {args.output_dir}")
        
        elif args.command == 'score':
            cmd = [
                sys.executable, os.path.join(src_dir, 'score_structures.py'),
                '--format', args.format,
                '--tables', args.tables,
                '--cutoff', str(args.cutoff),
                '--seq-sep', str(args.seq_sep)
            ]
            if args.list:
                cmd.extend(['--list', args.list])
            elif args.pdb:
                cmd.extend(['--pdb', args.pdb])
            elif args.folder:
                cmd.extend(['--folder', args.folder])
            else:
                raise ValueError("Must provide --list, --pdb, or --folder")
            if args.detailed:
                cmd.append('--detailed')
            if args.output:
                cmd.extend(['--output', args.output])
            subprocess.run(cmd, check=True)
            print(f"✓ Scoring complete")
        
        elif args.command == 'plot':
            cmd = [
                sys.executable, os.path.join(src_dir, 'plot_scores.py'),
                '--input-dir', args.input_dir,
                '--output-dir', args.output_dir
            ]
            if args.combined:
                cmd.append('--combined')
            subprocess.run(cmd, check=True)
            print(f"✓ Plots generated in {args.output_dir}")
        
        elif args.command == 'access':
            cmd = [
                sys.executable, os.path.join(src_dir, 'access_rna_structures.py'),
                '-n', str(args.number),
                '-f'] + args.formats + [
                '-o', args.output,
                '-w', str(args.workers)
            ]
            if args.all:
                cmd.append('--all')
            if args.rna_only:
                cmd.append('--rna-only')
            if args.info:
                cmd.append('--info')
            if args.list_only:
                cmd.append('--list-only')
            
            subprocess.run(cmd, check=True)
            print(f"✓ Structures downloaded to {args.output}")
            
            if args.validate:
                validate_cmd = [
                    sys.executable, os.path.join(src_dir, 'utils', 'validate_pdb_files.py'),
                    '--input-dir', os.path.join(args.output, 'pdb')
                ]
                subprocess.run(validate_cmd, check=True)
                print(f"✓ Validation complete")
    
    except subprocess.CalledProcessError as e:
        print(f"Error running script: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()

