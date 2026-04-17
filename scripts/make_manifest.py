#!/usr/bin/env python3
"""
生成 manifest.tsv

用法:
    python scripts/make_manifest.py /path/to/fastq_folder

生成 data/manifest.tsv，格式:
    sample-id, forward-absolute-filepath, reverse-absolute-filepath

支持以下文件命名格式:
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
    """从 forward fastq 路径提取样本名"""
    fname = os.path.basename(fwd_path)

    # 去掉常见后缀模式，提取样本名
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

    # fallback: 去掉 _R1/_R2 部分
    name = re.sub(r'[_ -][Rr][12].*\.fastq.*$', '', fname)
    return name


def find_fastq_pairs(folder):
    """扫描文件夹，找出所有 R1/R2 配对"""
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

        # 构建可能的 reverse 文件名
        rev_fname = re.sub(r'_R1_', '_R2_', fname)
        rev_fname = re.sub(r'_r1_', '_r2_', rev_fname)
        rev_fname = re.sub(r'_R1\.', '_R2.', rev_fname)
        rev_fname = re.sub(r'_r1\.', '_r2.', rev_fname)
        rev_fname = re.sub(r'_1\.', '_2.', rev_fname)

        rev_path = os.path.join(os.path.dirname(str(fwd)), rev_fname)

        if not os.path.exists(rev_path):
            # 尝试其他命名模式
            alt_rev = str(fwd).replace('_R1_', '_R2_').replace('_r1_', '_r2_')
            if os.path.exists(alt_rev):
                rev_path = alt_rev

        if not os.path.exists(rev_path):
            print(f"[WARN] 未找到 reverse 文件 for {fwd}", file=sys.stderr)
            continue

        sample_name = extract_sample_name(str(fwd))
        pairs.append((sample_name, str(fwd), rev_path))

    return pairs


def main():
    if len(sys.argv) < 2:
        print("用法: python make_manifest.py /path/to/fastq_folder")
        sys.exit(1)

    folder = sys.argv[1]
    out_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "manifest.tsv")

    print(f"扫描文件夹: {folder}")
    pairs = find_fastq_pairs(folder)

    if not pairs:
        print("[ERROR] 未找到任何 fastq 文件对！", file=sys.stderr)
        print("支持的文件名格式: *_R1_*.fastq.gz / *_R2_*.fastq.gz", file=sys.stderr)
        sys.exit(1)

    print(f"找到 {len(pairs)} 个样本")

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["sample-id", "forward-absolute-filepath", "reverse-absolute-filepath"])
        for sample, fwd, rev in pairs:
            writer.writerow([sample, fwd, rev])

    print(f"已生成: {out_path}")


if __name__ == "__main__":
    main()
