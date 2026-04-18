#!/usr/bin/env python3
# Author: Rong Wu

"""
Generate manifest.tsv for QIIME2 pipeline.

Usage:
    python scripts/make_manifest.py /path/to/fastq_folder

Output: data/manifest.tsv
    sample-id, forward-absolute-filepath, reverse-absolute-filepath

Supports:
    sample_NAME_L001_R1_001.fastq.gz
    sample_NAME_S1_L001_R1_001.fastq.gz
    NAME_R1.fastq.gz / NAME_R2.fastq.gz
    *_R1_*.fastq.gz / *_R2_*.fastq.gz
"""

import os
import sys
import re
import csv

try:
    from pathlib import Path
except ImportError:
    from pathlib2 import Path  # Python 2 fallback


def extract_sample_name(fwd_path):
    """Extract sample name from forward fastq path."""
    fname = os.path.basename(fwd_path)

    patterns = [
        r'^(.+?)_L\d+',
        r'^(.+?)_S\d+',
        r'^(.+?)_R1',
        r'^(.+?)_r1',
    ]
    for pat in patterns:
        m = re.match(pat, fname, re.IGNORECASE)
        if m:
            return m.group(1)

    name = re.sub(r'[_ -][Rr][12].*\.fastq.*$', '', fname)
    return name


def find_fastq_pairs(folder):
    """Scan folder and find all R1/R2 pairs."""
    folder = os.path.abspath(folder)
    fwd_files = sorted(Path(folder).glob("**/*_R1_*.fastq.gz"))
    if not fwd_files:
        fwd_files = sorted(Path(folder).glob("**/*_r1_*.fastq.gz"))
    if not fwd_files:
        fwd_files = sorted(Path(folder).glob("**/*_R1.fastq.gz"))
    if not fwd_files:
        fwd_files = sorted(Path(folder).glob("**/*_R1.fq.gz"))
    if not fwd_files:
        fwd_files = sorted(Path(folder).glob("**/*_1.fastq.gz"))
    if not fwd_files:
        fwd_files = sorted(Path(folder).glob("**/*.fastq.gz"))

    pairs = []
    for fwd in fwd_files:
        fname = os.path.basename(str(fwd))

        rev_fname = re.sub(r'_R1_', '_R2_', fname)
        rev_fname = re.sub(r'_r1_', '_r2_', rev_fname)
        rev_fname = re.sub(r'_R1\.', '_R2.', rev_fname)
        rev_fname = re.sub(r'_r1\.', '_r2.', rev_fname)
        rev_fname = re.sub(r'_1\.', '_2.', rev_fname)

        rev_path = os.path.join(os.path.dirname(str(fwd)), rev_fname)

        if not os.path.exists(rev_path):
            alt_rev = str(fwd).replace('_R1_', '_R2_').replace('_r1_', '_r2_')
            if os.path.exists(alt_rev):
                rev_path = alt_rev

        if not os.path.exists(rev_path):
            print(f"[WARN] reverse not found for {fwd}", file=sys.stderr)
            continue

        sample_name = extract_sample_name(str(fwd))
        pairs.append((sample_name, str(fwd), rev_path))

    return pairs


def main():
    if len(sys.argv) < 2:
        print("Usage: python make_manifest.py /path/to/fastq_folder")
        sys.exit(1)

    folder = sys.argv[1]
    out_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "manifest.tsv")

    print(f"Scanning: {folder}")
    pairs = find_fastq_pairs(folder)

    if not pairs:
        print("[ERROR] No fastq pairs found!", file=sys.stderr)
        print("Supported: *_R1_*.fastq.gz / *_R2_*.fastq.gz", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(pairs)} samples")

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
        for sample, fwd, rev in pairs:
            writer.writerow([sample, fwd, rev])

    print(f"Written: {out_path}")


if __name__ == "__main__":
    main()
