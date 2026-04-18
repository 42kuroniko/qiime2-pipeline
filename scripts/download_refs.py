#!/usr/bin/env python3

"""
Download GreenGenes2 reference database and classifier.
GreenGenes2 2024.09 (4mer proxy classifier)

Usage:
    python scripts/download_refs.py
"""

import os
import urllib.request
import sys

BASE_URL = "https://data.qiime2.org/2024.2/common"
ENV_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "envs")
os.makedirs(ENV_DIR, exist_ok=True)

FILES = {
    "gg_2024_09_12_4mer_seqs.qza": f"{BASE_URL}/gg_2024_09_12_4mer_seqs.qza",
    "gg_2024_09_12_4mer_tax.qza": f"{BASE_URL}/gg_2024_09_12_4mer_tax.qza",
    "gg_2024_09_12_4mer_proxy_tax_classifier.qza": (
        f"{BASE_URL}/gg_2024_09_12_4mer_proxy_tax_classifier.qza"
    ),
}


def download_file(url, dest):
    if os.path.exists(dest):
        size = os.path.getsize(dest) / (1024 * 1024)
        print(f"[SKIP] exists: {dest} ({size:.1f} MB)")
        return

    print(f"[DOWN] {url}")
    print(f"  -> {dest}")

    try:
        urllib.request.urlretrieve(url, dest)
        size = os.path.getsize(dest) / (1024 * 1024)
        print(f"[DONE] {dest} ({size:.1f} MB)")
    except Exception as e:
        os.remove(dest)
        print(f"[FAIL] {e}", file=sys.stderr)
        sys.exit(1)


def main():
    for fname, url in FILES.items():
        dest = os.path.join(ENV_DIR, fname)
        download_file(url, dest)

    print("\nDone. GreenGenes2 files ready.")


if __name__ == "__main__":
    main()
