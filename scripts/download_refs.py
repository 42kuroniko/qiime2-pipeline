#!/usr/bin/env python3
"""
下载 GreenGenes2 参考数据库和分类器
GreenGenes2 2024.09 版本 (4mer proxy classifier)

用法:
    python scripts/download_refs.py

下载内容:
    - envs/gg_2024_09_12_4mer_seqs.qza
    - envs/gg_2024_09_12_4mer_tax.qza
    - envs/gg_2024_09_12_4mer_proxy_tax_classifier.qza
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
        print(f"[SKIP] 已存在: {dest} ({size:.1f} MB)")
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

    print("\n所有参考文件下载完成！")
    print("现在可以开始运行 Snakemake 流程了：")


if __name__ == "__main__":
    main()
