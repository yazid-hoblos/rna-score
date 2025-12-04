#!/usr/bin/env python3
"""
Download RNA structures from PDB in PDB and mmCIF formats.
Allows setting maximum number of structures to download.
Allows filtering for RNA-only structures.
"""

import sys
import time
import argparse
import requests
from pathlib import Path
import json
from typing import List, Dict, Optional
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm


class RNAStructureDownloader:
    """Download RNA structures from PDB database."""
    
    def __init__(self, output_dir: str = "rna_structures", 
                 max_workers: int = 5,
                 delay: float = 0.5):
        """
        Initialize downloader.
        
        Args:
            output_dir: Directory to save structures
            max_workers: Number of parallel downloads
            delay: Delay between requests (seconds)
        """
        self.output_dir = Path(output_dir)
        self.pdb_dir = self.output_dir / "pdb"
        self.cif_dir = self.output_dir / "mmcif"
        self.max_workers = max_workers
        self.delay = delay
        
        # Create directories
        self.pdb_dir.mkdir(parents=True, exist_ok=True)
        self.cif_dir.mkdir(parents=True, exist_ok=True)
        
        # PDB API endpoints
        self.search_url = "https://search.rcsb.org/rcsbsearch/v2/query"
        self.download_base = "https://files.rcsb.org/download/"
        
    def search_rna_structures(self, max_results: Optional[int] = None) -> List[str]:
        """
        Search RNA structures in PDB using external query file.
        Args:
            max_results: Maximum number of results to return
        Returns:
            List of PDB IDs
        """
        query_path = Path(__file__).parent.parent / "queries" / "rna_structures.json"
        with open(query_path, "r") as f:
            query = json.load(f)
        # Update max_results in query
        if max_results is not None:
            query["request_options"]["paginate"]["rows"] = max_results
        print(f"Searching for RNA structures in PDB...")
        try:
            response = requests.post(
                self.search_url,
                json=query,
                headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()
            data = response.json()
            pdb_ids = [item["identifier"] for item in data.get("result_set", [])]
            print(f"Found {len(pdb_ids)} RNA structures")
            if max_results and len(pdb_ids) > max_results:
                pdb_ids = pdb_ids[:max_results]
                print(f"Limited to {max_results} structures as requested")
            return pdb_ids
        except requests.exceptions.RequestException as e:
            print(f"Error searching PDB: {e}")
            return []
    
    def get_rna_only_structures(self, max_results: Optional[int] = None) -> List[str]:
        """
        Get structures that contain ONLY RNA (no protein or DNA) using external query file.
        Args:
            max_results: Maximum number of results
        Returns:
            List of PDB IDs
        """
        query_path = Path(__file__).parent.parent / "queries" / "rna_only.json"
        with open(query_path, "r") as f:
            query = json.load(f)
        # Update max_results in query
        if max_results is not None:
            query["request_options"]["paginate"]["rows"] = max_results
        print(f"Searching for RNA-only structures...")
        try:
            response = requests.post(
                self.search_url,
                json=query,
                headers={"Content-Type": "application/json"}
            )
            response.raise_for_status()
            data = response.json()
            pdb_ids = [item["identifier"] for item in data.get("result_set", [])]
            print(f"Found {len(pdb_ids)} RNA-only structures")
            if max_results and len(pdb_ids) > max_results:
                pdb_ids = pdb_ids[:max_results]
                print(f"Limited to {max_results} structures")
            return pdb_ids
        except requests.exceptions.RequestException as e:
            print(f"Error searching PDB: {e}")
            return []
    
    def download_structure(self, pdb_id: str, format: str = "pdb") -> bool:
        """
        Download a single structure.
        
        Args:
            pdb_id: PDB identifier
            format: Format to download ('pdb' or 'cif')
            
        Returns:
            True if successful, False otherwise
        """
        pdb_id = pdb_id.lower()
        
        # Determine file extension and directory
        if format == "pdb":
            ext = ".pdb"
            output_dir = self.pdb_dir
        else:
            ext = ".cif"
            output_dir = self.cif_dir
        
        output_file = output_dir / f"{pdb_id}{ext}"
        
        # Skip if already downloaded
        if output_file.exists():
            return True
        
        # Download URL
        url = f"{self.download_base}{pdb_id}{ext}"
        
        try:
            response = requests.get(url, timeout=30)
            if response.status_code != 200:
                print(f"Download failed for {pdb_id}.{format}: HTTP {response.status_code} - {response.reason}")
                return False
            # Save file
            output_file.write_bytes(response.content)
            time.sleep(self.delay)  # otherwise we get blocked
            return True
        except requests.exceptions.RequestException as e:
            print(f"Download failed for {pdb_id}.{format}: {e}")
            return False
        except Exception as e:
            print(f"Unexpected error for {pdb_id}.{format}: {e}")
            return False
    
    def download_batch(self, pdb_ids: List[str], formats: List[str] = ["pdb", "cif"]):
        """
        Download multiple structures in parallel.
        
        Args:
            pdb_ids: List of PDB IDs
            formats: List of formats to download
        """
        # Create download tasks
        tasks = []
        for pdb_id in pdb_ids:
            for fmt in formats:
                tasks.append((pdb_id, fmt))
        
        successful = 0
        failed = 0
        failed_ids = []
        
        print(f"\nDownloading {len(pdb_ids)} structures in {formats} format(s)...")
        print(f"Output directory: {self.output_dir}")
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            futures = {
                executor.submit(self.download_structure, pdb_id, fmt): (pdb_id, fmt)
                for pdb_id, fmt in tasks
            }
            
            with tqdm(total=len(tasks), desc="Downloading") as pbar:
                for future in as_completed(futures):
                    pdb_id, fmt = futures[future]
                    try:
                        success = future.result()
                        if success:
                            successful += 1
                        else:
                            failed += 1
                            failed_ids.append(f"{pdb_id}.{fmt}")
                            print(f"Debug: Download failed for {pdb_id}.{fmt}. See above for details.")
                    except Exception as e:
                        failed += 1
                        failed_ids.append(f"{pdb_id}.{fmt}")
                        print(f"\nError downloading {pdb_id}.{fmt}: {e}")
                    
                    pbar.update(1)
                    pbar.set_postfix({"Success": successful, "Failed": failed})
        
        print(f"\n{'='*50}")
        print(f"Download Summary:")
        print(f"  Total attempted: {len(tasks)}")
        print(f"  Successful: {successful}")
        print(f"  Failed: {failed}")
        
        if failed_ids:
            print(f"\nFailed downloads:")
            for fid in failed_ids[:10]:  # Show first 10
                print(f"  - {fid}")
            if len(failed_ids) > 10:
                print(f"  ... and {len(failed_ids) - 10} more")
            
            # Save failed IDs to file
            failed_file = self.output_dir / "failed_downloads.txt"
            failed_file.write_text("\n".join(failed_ids))
            print(f"\nFailed IDs saved to: {failed_file}")
    
    def get_structure_info(self, pdb_ids: List[str]) -> Dict:
        """
        Get information about structures.
        
        Args:
            pdb_ids: List of PDB IDs
            
        Returns:
            Dictionary with structure information
        """
        info = []
        
        print("\nFetching structure information...")
        
        for pdb_id in tqdm(pdb_ids[:10], desc="Getting info"):  # Sample first 10
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
            
            try:
                response = requests.get(url)
                if response.status_code == 200:
                    data = response.json()
                    info.append({
                        "pdb_id": pdb_id,
                        "title": data.get("struct", {}).get("title", "N/A"),
                        "resolution": data.get("rcsb_entry_info", {}).get("resolution_combined", "N/A"),
                        "method": data.get("exptl", [{}])[0].get("method", "N/A") if data.get("exptl") else "N/A"
                    })
            except:
                continue
        
        return info


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Download RNA structures from PDB database",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Download 100 RNA structures in both formats
  python download_rna_structures.py -n 100
  
  # Download 50 RNA-only structures in PDB format only
  python download_rna_structures.py -n 50 --rna-only --formats pdb
  
  # Download all RNA structures (be careful!)
  python download_rna_structures.py --all
        """
    )
    
    parser.add_argument(
        "-n", "--number",
        type=int,
        default=10,
        help="Maximum number of structures to download (default: 10)"
    )
    
    parser.add_argument(
        "--all",
        action="store_true",
        help="Download all available RNA structures (overrides -n)"
    )
    
    parser.add_argument(
        "--rna-only",
        action="store_true",
        help="Download only pure RNA structures (no protein/DNA)"
    )
    
    parser.add_argument(
        "-f", "--formats",
        nargs="+",
        choices=["pdb", "cif"],
        default=["pdb", "cif"],
        help="File formats to download (default: both)"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="rna_structures",
        help="Output directory (default: rna_structures)"
    )
    
    parser.add_argument(
        "-w", "--workers",
        type=int,
        default=5,
        help="Number of parallel downloads (default: 5)"
    )
    
    parser.add_argument(
        "--info",
        action="store_true",
        help="Show information about structures before downloading"
    )
    
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="Only list PDB IDs without downloading"
    )
    
    args = parser.parse_args()
    
    # Initialize downloader
    downloader = RNAStructureDownloader(
        output_dir=args.output,
        max_workers=args.workers
    )
    
    # Determine max results
    max_results = None if args.all else args.number
    
    # Search for structures
    if args.rna_only:
        pdb_ids = downloader.get_rna_only_structures(max_results)
    else:
        pdb_ids = downloader.search_rna_structures(max_results)
    
    if not pdb_ids:
        print("No structures found!")
        return 1
    
    # List only mode
    if args.list_only:
        print("\nPDB IDs:")
        for i, pdb_id in enumerate(pdb_ids, 1):
            print(f"{i:4d}. {pdb_id}")
        
        # Save to file
        list_file = Path(args.output) / "pdb_ids.txt"
        list_file.parent.mkdir(parents=True, exist_ok=True)
        list_file.write_text("\n".join(pdb_ids))
        print(f"\nPDB IDs saved to: {list_file}")
        return 0
    
    # Show info if requested
    if args.info:
        info = downloader.get_structure_info(pdb_ids)
        if info:
            print("\nSample structures:")
            print(f"{'PDB ID':<8} {'Resolution':<12} {'Method':<15} Title")
            print("-" * 80)
            for item in info:
                title = item['title'][:40] + "..." if len(item['title']) > 40 else item['title']
                print(f"{item['pdb_id']:<8} {str(item['resolution']):<12} {item['method']:<15} {title}")
    
    # Download structures
    downloader.download_batch(pdb_ids, formats=args.formats)
    
    # Save PDB ID list
    list_file = downloader.output_dir / "downloaded_ids.txt"
    list_file.write_text("\n".join(pdb_ids))
    print(f"\nPDB IDs saved to: {list_file}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())