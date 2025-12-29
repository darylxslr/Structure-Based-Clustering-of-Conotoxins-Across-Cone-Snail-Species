#!/usr/bin/env python3
"""
Download conotoxin-related PDB files from RCSB until the pdb_files/ directory
contains 20 files in total. Skips files that already exist.

Usage:
    python scripts/download_more_pdbs.py

Requires `requests` (already in requirements.txt).
"""

import time
from pathlib import Path
import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]
PDB_DIR = PROJECT_ROOT / "pdb_files"
PDB_DIR.mkdir(parents=True, exist_ok=True)

SEARCH_URLS = [
    "https://search.rcsb.org/rcsbsearch/v1/query",
    "https://search.rcsb.org/rcsbsearch/v2/query",
]
DOWNLOAD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def search_conotoxins(size=100):
    # Try multiple RCSB search endpoints and return identifiers
    ids = []
    # try multiple query phrases to increase chance of matches
    phrases = [
        "conotoxin",
        "conopeptide",
        "alpha-conotoxin",
        "omega-conotoxin",
        "conotoxin peptide",
    ]

    for phrase in phrases:
        for url in SEARCH_URLS:
            services = ["text", "full_text"]
            body = None
            for service in services:
                body = {
                    "query": {
                        "type": "terminal",
                        "service": service,
                        "parameters": {"value": phrase}
                    },
                    "return_type": "entry",
                }

            try:
                resp = requests.post(url, json=body, timeout=15)
            except requests.RequestException:
                continue

            if resp.status_code >= 400:
                continue

            try:
                data = resp.json()
            except Exception:
                continue

            if isinstance(data, dict):
                if "data" in data and isinstance(data["data"], list):
                    ids = [r.get("identifier") for r in data.get("data", []) if r.get("identifier")]
                elif "result_set" in data:
                    ids = [r.get("identifier") for r in data.get("result_set", []) if r.get("identifier")]

            if ids:
                print(f"Query '{phrase}' on {url} returned {len(ids)} ids")
                break

        if ids:
            break

    # dedupe and limit
    seen = []
    collected = []
    for phrase in phrases:
        for url in SEARCH_URLS:
            for service in services:
                body = {
                    "query": {
                        "type": "terminal",
                        "service": service,
                        "parameters": {"value": phrase}
                    },
                    "return_type": "entry",
                }

                try:
                    resp = requests.post(url, json=body, timeout=15)
                except requests.RequestException:
                    continue

                if resp.status_code >= 400:
                    continue

                try:
                    data = resp.json()
                except Exception:
                    continue

                ids = []
                if isinstance(data, dict):
                    if "data" in data and isinstance(data["data"], list):
                        ids = [r.get("identifier") for r in data.get("data", []) if r.get("identifier")]
                    elif "result_set" in data:
                        ids = [r.get("identifier") for r in data.get("result_set", []) if r.get("identifier")]

                for i in ids:
                    if i and i not in collected:
                        collected.append(i)
                    if len(collected) >= size:
                        break

                if collected and len(collected) >= size:
                    break

            if collected and len(collected) >= size:
                break

        if collected and len(collected) >= size:
            break

    ids = collected
    return ids


def download_pdb(pdb_id, out_path: Path):
    url = DOWNLOAD_URL.format(pdb_id=pdb_id)
    resp = requests.get(url, timeout=15)
    resp.raise_for_status()
    out_path.write_bytes(resp.content)


def main(target_total=20):
    existing = {p.stem.upper() for p in PDB_DIR.glob("*.pdb")}
    if len(existing) >= target_total:
        print(f"Already have {len(existing)} PDB files (>= {target_total}). Nothing to do.")
        return

    needed = target_total - len(existing)
    print(f"Need {needed} more PDB files to reach {target_total} total.")

    ids = search_conotoxins(size=500)
    print(f"Found {len(ids)} search results for 'conotoxin' (using RCSB).")

    downloaded = 0
    for pdb_id in ids:
        pid = pdb_id.upper()
        if pid in existing:
            continue
        out_path = PDB_DIR / f"{pid}.pdb"
        try:
            download_pdb(pid, out_path)
            print(f"Downloaded {pid} -> {out_path.name}")
            downloaded += 1
            existing.add(pid)
        except Exception as e:
            print(f"Failed to download {pid}: {e}")

        if len(existing) >= target_total:
            break

        # be polite to the API
        time.sleep(0.2)

    print(f"Downloaded {downloaded} new files; total now {len(existing)}.")


if __name__ == "__main__":
    main()
