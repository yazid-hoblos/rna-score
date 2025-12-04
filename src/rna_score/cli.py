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
    extract_parser.add_argument('--list', type=str, help='File with structure IDs or filenames')
    extract_parser.add_argument('--folder', type=str, help='Directory with PDB/mmCIF files')
    extract_parser.add_argument('--pdb', type=str, help='Single PDB/mmCIF file')
    extract_parser.add_argument('--format', choices=['pdb', 'mmcif'], default='pdb', help='Structure format')
    extract_parser.add_argument('--atom-mode', nargs='+', default=['"C3\'"'], help='Atom selection mode(s)')
    extract_parser.add_argument('--chains', nargs='+', help='Specific chains to process')
    extract_parser.add_argument('--dist-mode', choices=['intra', 'inter'], default='intra')
    extract_parser.add_argument('--cutoff', type=float, default=20.0)
    extract_parser.add_argument('--seq-sep', type=int, default=4)
    extract_parser.add_argument('--bin-width', type=float, default=1.0)
    extract_parser.add_argument('--cores', type=int, help='Number of cores to use')
    extract_parser.add_argument('--out-dir', type=str, default='dist_data')
    extract_parser.add_argument('--save-detailed', action='store_true', help='Export detailed CSV log')
    extract_parser.add_argument('--all-models', action='store_true', help='Process all NMR models')

    # Train scoring function
    train_parser = subparsers.add_parser('train', help='Train scoring function from distances')
    train_parser.add_argument('--hist-dir', type=str, required=True, help='Directory with histogram files')
    train_parser.add_argument('--out-dir', type=str, default='training_output')
    train_parser.add_argument('--max-score', type=float, default=10.0)
    train_parser.add_argument('--pseudocount', type=float, default=1e-6)
    train_parser.add_argument('--cutoff', type=float, default=20.0)
    train_parser.add_argument('--bin-width', type=float, default=1.0)

    # Score structures
    score_parser = subparsers.add_parser('score', help='Score RNA structures')
    score_parser.add_argument('--list', type=str, help='File with list of structure file paths')
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
        if args.command == 'extract':
            cmd = [
                sys.executable, os.path.join(src_dir, 'extract_distances.py'),
                '--format', args.format,
                '--atom-mode'] + args.atom_mode + [
                '--dist-mode', args.dist_mode,
                '--cutoff', str(args.cutoff),
                '--seq-sep', str(args.seq_sep),
                '--bin-width', str(args.bin_width),
                '--method', 'histogram',
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
                '--input-dir', args.hist_dir,
                '--output-dir', args.out_dir,
                '--max-score', str(args.max_score),
                '--pseudocount', str(args.pseudocount),
                '--cutoff', str(args.cutoff),
                '--bin-width', str(args.bin_width)
            ]
            subprocess.run(cmd, check=True)
            print(f"✓ Training complete. Output in {args.out_dir}")
        
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

