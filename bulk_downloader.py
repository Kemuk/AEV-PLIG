#!/usr/bin/env python3
"""
Simple script to bulk download PDB files and SDF files from a CSV containing PDB codes.
Uses RCSB APIs to get properly formatted SDF files instead of OpenBabel conversion.
Also creates a dataset.csv mapping file.
"""

import pandas as pd
import requests
import os
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

def get_ligands_from_entry(pdb_code):
    """Get ligand component IDs from a PDB entry using RCSB GraphQL API."""
    try:
        # Use GraphQL to get nonpolymer entities (ligands)
        query = """{ entry(entry_id: "%s") { nonpolymer_entities { pdbx_entity_nonpoly { comp_id name } } } }""" % pdb_code.lower()
        
        query_url = f"https://data.rcsb.org/graphql?query={query}"
        response = requests.get(query_url, timeout=30)
        
        if response.status_code == 200:
            data = response.json()
            ligands = []
            
            if 'data' in data and data['data'] and 'entry' in data['data']:
                entry_data = data['data']['entry']
                if entry_data and 'nonpolymer_entities' in entry_data:
                    nonpolymer_entities = entry_data['nonpolymer_entities']
                    
                    for entity in nonpolymer_entities:
                        if 'pdbx_entity_nonpoly' in entity:
                            comp_id = entity['pdbx_entity_nonpoly'].get('comp_id')
                            if comp_id:
                                # Filter out common solvents/ions (optional)
                                common_solvents = ['HOH', 'WAT', 'NA', 'CL', 'SO4', 'PO4', 'GOL', 'EDO']
                                if comp_id not in common_solvents:
                                    ligands.append(comp_id)
                    
            return ligands[:2]  # Limit to first 2 ligands
        else:
            # Debug: print status for troubleshooting
            print(f"  GraphQL query failed for {pdb_code}: HTTP {response.status_code}")
        return []
    except Exception as e:
        print(f"  Error getting ligands for {pdb_code}: {e}")
        return []

def download_sdf_from_modelserver(pdb_code, ligand_comp_ids, sdf_dir):
    """Download SDF files using RCSB ModelServer API with correct parameters."""
    try:
        if not ligand_comp_ids:
            return False
        
        for i, comp_id in enumerate(ligand_comp_ids):
            # Use ModelServer API with correct parameters
            url = f"https://models.rcsb.org/v1/{pdb_code.lower()}/ligand"
            
            # Request SDF format with correct parameters
            params = {
                'label_comp_id': comp_id,
                'encoding': 'sdf',
                'copy_all_categories': 'false'
            }
            
            response = requests.get(url, params=params, timeout=30)
            
            if response.status_code == 200 and response.text.strip():
                file_path = sdf_dir / f"{pdb_code}.sdf"
                with open(file_path, 'w') as f:
                    f.write(response.text)
                return True  # Successfully downloaded first ligand
            else:
                print(f"  Failed to download SDF for {pdb_code}/{comp_id}: HTTP {response.status_code}")
            
        return False
                
    except Exception:
        return False

def download_single_molecule(pdb_code, pdb_dir, sdf_dir):
    """Download PDB and SDF files for a single molecule using APIs."""
    pdb_success = False
    sdf_success = False
    
    # Download PDB file
    pdb_url = f"https://files.rcsb.org/download/{pdb_code.lower()}.pdb"
    try:
        response = requests.get(pdb_url, timeout=30)
        if response.status_code == 200:
            pdb_file_path = pdb_dir / f"{pdb_code}.pdb"
            with open(pdb_file_path, 'w') as f:
                f.write(response.text)
            pdb_success = True
    except Exception as e:
        print(f"  Error downloading PDB {pdb_code}: {e}")
    
    # Get ligands and download SDF using ModelServer API
    ligand_comp_ids = get_ligands_from_entry(pdb_code)
    if ligand_comp_ids:
        sdf_success = download_sdf_from_modelserver(pdb_code, ligand_comp_ids, sdf_dir)
    
    return {
        'pdb_code': pdb_code,
        'pdb_success': pdb_success,
        'sdf_success': sdf_success
    }

def download_from_csv(csv_file, pdb_column=None, output_dir="./structures", max_workers=5):
    """
    Download PDB and SDF files from CSV using RCSB APIs and create dataset mapping.
    
    Args:
        csv_file: Path to your CSV file
        pdb_column: Column name with PDB codes (auto-detects if None)
        output_dir: Where to save the files
        max_workers: Number of concurrent downloads
    """
    
    # Read CSV
    df = pd.read_csv(csv_file)
    print(f"CSV columns: {list(df.columns)}")
    
    # Auto-detect PDB column
    if pdb_column is None:
        pdb_columns = [col for col in df.columns if 'pdb' in col.lower()]
        if pdb_columns:
            pdb_column = pdb_columns[0]
        else:
            pdb_column = df.columns[0]  # Use first column
    
    print(f"Using column: {pdb_column}")
    
    # Get PDB codes
    pdb_codes = df[pdb_column].dropna().astype(str).str.strip().str.upper().tolist()
    print(f"Found {len(pdb_codes)} PDB codes")
    print(f"First few: {pdb_codes[:5]}")
    
    # Create output directories
    pdb_dir = Path(output_dir) / "pdb_files"
    sdf_dir = Path(output_dir) / "sdf_files"
    pdb_dir.mkdir(parents=True, exist_ok=True)
    sdf_dir.mkdir(parents=True, exist_ok=True)
    
    # Download files with threading using RCSB APIs
    results = []
    print(f"\nDownloading PDB files and SDF files using RCSB APIs with {max_workers} workers...")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(download_single_molecule, code, pdb_dir, sdf_dir): code 
                  for code in pdb_codes}
        
        for future in tqdm(futures, desc="Downloading molecules"):
            try:
                result = future.result(timeout=120)  # Increased timeout for API calls
                results.append(result)
            except Exception as e:
                pdb_code = futures[future]
                print(f"Failed to download {pdb_code}: {e}")
                results.append({
                    'pdb_code': pdb_code,
                    'pdb_success': False,
                    'sdf_success': False
                })
    
    # Create dataset CSV mapping
    dataset_rows = []
    successful_pdb = sum(1 for r in results if r['pdb_success'])
    successful_sdf = sum(1 for r in results if r['sdf_success'])
    
    for result in results:
        pdb_code = result['pdb_code']
        
        # Only include molecules where we have at least one file
        if result['pdb_success'] or result['sdf_success']:
            dataset_row = {
                'unique_id': pdb_code.lower(),
                'sdf_file': None,
                'pdb_file': None
            }
            
            if result['pdb_success']:
                dataset_row['pdb_file'] = f"pdb_files/{pdb_code}.pdb"
                
            if result['sdf_success']:
                dataset_row['sdf_file'] = f"sdf_files/{pdb_code}.sdf"
            
            dataset_rows.append(dataset_row)
    
    dataset_df = pd.DataFrame(dataset_rows)
    dataset_path = Path(output_dir) / "dataset.csv"
    dataset_df.to_csv(dataset_path, index=False)
    
    # Summary
    total = len(pdb_codes)
    both_success = sum(1 for r in results if r['pdb_success'] and r['sdf_success'])
    
    print(f"\nDownload complete!")
    print(f"Total molecules processed: {total}")
    print(f"PDB files downloaded: {successful_pdb} ({successful_pdb/total*100:.1f}%)")
    print(f"SDF files downloaded: {successful_sdf} ({successful_sdf/total*100:.1f}%)")
    print(f"Both PDB and SDF: {both_success} ({both_success/total*100:.1f}%)")
    print(f"Dataset CSV saved to: {dataset_path}")
    print(f"Dataset contains {len(dataset_df)} molecules with downloaded files")
    print()
    print("Notes:")
    print("- SDF files downloaded using RCSB ModelServer API for proper formatting")
    print("- Only structures with actual ligands will have SDF files")
    print("- Some structures may only have solvents/ions (no SDF generated)")
    
    # Show failed downloads
    failed = [r['pdb_code'] for r in results if not (r['pdb_success'] or r['sdf_success'])]
    if failed:
        print(f"\nCompletely failed downloads: {failed[:10]}{'...' if len(failed) > 10 else ''}")
    
    return dataset_df

# Usage
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Simple bulk downloader for PDB/SDF files using RCSB APIs")
    parser.add_argument("csv_file", help="Path to CSV file containing PDB codes")
    parser.add_argument("--pdb_column", help="Name of column containing PDB codes")
    parser.add_argument("--output_dir", default="./structures", help="Output directory")
    parser.add_argument("--max_workers", type=int, default=5, help="Number of concurrent downloads")
    
    args = parser.parse_args()
    
    dataset_df = download_from_csv(
        csv_file=args.csv_file,
        pdb_column=args.pdb_column,
        output_dir=args.output_dir,
        max_workers=args.max_workers
    )
    
    print("\nDataset CSV preview:")
    print(dataset_df.head())