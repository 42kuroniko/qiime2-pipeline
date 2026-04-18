#!/usr/bin/env bash
# Author: Rong Wu
# Setup conda env + download GreenGenes2 references.
# Usage: bash scripts/setup_env.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_DIR="$PROJECT_DIR/envs"

echo "=========================================="
echo "qiime2-pipeline setup"
echo "=========================================="
echo ""

if ! command -v conda &> /dev/null; then
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
    else
        echo "[ERROR] conda or mamba not found"
        echo "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
else
    CONDA_CMD="conda"
fi
echo "[OK] conda/mamba: $CONDA_CMD"

ENV_FILE="$ENV_DIR/qiime2.yaml"
if conda env list | grep -q "^qiime2 "; then
    echo "[SKIP] conda env 'qiime2' already exists"
else
    echo ""
    echo "[INFO] Creating conda env (5-10 min)..."
    $CONDA_CMD env create -f "$ENV_FILE"
fi

echo ""
echo "[INFO] Downloading GreenGenes2 reference data..."
source "$CONDA_CMD" activate qiime2 2>/dev/null || conda activate qiime2

python "$SCRIPT_DIR/download_refs.py"

echo ""
echo "=========================================="
echo "Setup complete."
echo "=========================================="
echo ""
echo "Activate: conda activate qiime2"
echo "Test:     qiime --help"
echo "Run:      snakemake -np  # preview"
echo "          snakemake --cores 4"
