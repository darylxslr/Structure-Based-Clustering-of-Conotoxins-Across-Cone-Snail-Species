#!/usr/bin/env python3
"""
Simple downloader to fetch PDB files from RCSB for a list of toxin names.
Saves files to the project's `pdb_files/` directory using the toxin name as filename.

Usage:
    python scripts/download_pdbs.py

Requirements:
    pip install requests

If you already know the PDB IDs for some names, add them to the `mapping` dict.
"""

import os
from pathlib import Path
import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = PROJECT_ROOT / "pdb_files"
PDB_DIR.mkdir(parents=True, exist_ok=True)

SEARCH_URLS = [
    "https://search.rcsb.org/rcsbsearch/v1/query",
    "https://search.rcsb.org/rcsbsearch/v2/query",
]


def search_pdb_ids(name):
    """Try multiple RCSB search endpoints to find PDB IDs for `name`.

    Returns a list of PDB identifiers (may be empty).
    """
    services = ["text", "full_text"]

    for service in services:
        body = {
            "query": {
                "type": "terminal",
                "service": service,
                "parameters": {"value": name},
            },
            "return_type": "entry",
        }

        for url in SEARCH_URLS:
            try:
                resp = requests.post(url, json=body, timeout=10)
            except requests.RequestException:
                continue

        if resp.status_code == 404:
            # endpoint not available, try next
            continue

        try:
            resp.raise_for_status()
        except requests.HTTPError:
            continue

        data = resp.json()
        # different versions may return under different keys
        ids = []
        if isinstance(data, dict):
            # v1 uses 'result_set'
            if "result_set" in data:
                ids = [r.get("identifier") for r in data.get("result_set", []) if r.get("identifier")]
            # v2 uses 'result_set' or 'data'
            elif "data" in data and isinstance(data["data"], list):
                ids = [r.get("identifier") for r in data.get("data", []) if r.get("identifier")]

        if ids:
            return ids

    return []


def download_pdb(pdb_id, out_path: Path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    resp = requests.get(url, timeout=10)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)


def main():
    names = ["MVIIA", "GI", "ImI", "GIIIA", "EVIA", "PVIIA"]

    # Optional: fill known mappings name -> pdb_id to avoid searching.
    mapping = {
        # e.g. "MVIIA": "1XYZ",
    }

    for name in names:
        pdb_id = mapping.get(name)
        if not pdb_id:
            # Try multiple search phrases to improve chances of finding peptide entries
            search_phrases = [
                name,
                f"conotoxin {name}",
                f"omega-conotoxin {name}",
                f"{name} conotoxin",
            ]
            ids = []
            for phrase in search_phrases:
                ids = search_pdb_ids(phrase)
                if ids:
                    pdb_id = ids[0]
                    print(f"Found {pdb_id} for {name} using query {phrase!r}")
                    break

            if not ids:
                print(f"No PDB entries found for {name!r} (tried {search_phrases})")
                continue

        out_path = PDB_DIR / f"{name}.pdb"
        try:
            download_pdb(pdb_id, out_path)
            print(f"Saved {out_path}")
        except Exception as e:
            print(f"Failed to download {pdb_id} for {name}: {e}")


if __name__ == "__main__":
    main()
