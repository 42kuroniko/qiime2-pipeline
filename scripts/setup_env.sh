#!/usr/bin/env bash
# =============================================================================
# 初始化脚本: 创建 conda 环境 + 下载 GreenGenes2 参考数据库
# 用法: bash scripts/setup_env.sh
# =============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
ENV_DIR="$PROJECT_DIR/envs"

echo "=========================================="
echo "qiime2-pipeline 环境部署"
echo "=========================================="
echo ""

# 1. 检查 conda
if ! command -v conda &> /dev/null; then
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
    else
        echo "[ERROR] 未找到 conda 或 mamba"
        echo "请先安装 Miniconda: https://docs.conda.io/en/latest/miniconda.html"
        exit 1
    fi
else
    CONDA_CMD="conda"
fi
echo "[OK] conda/mamba: $CONDA_CMD"

# 2. 创建环境
ENV_FILE="$ENV_DIR/qiime2.yaml"
if conda env list | grep -q "^qiime2 "; then
    echo "[SKIP] conda 环境 'qiime2' 已存在"
else
    echo ""
    echo "[INFO] 创建 conda 环境 (约需 5-10 分钟)..."
    $CONDA_CMD env create -f "$ENV_FILE"
fi

# 3. 下载 GreenGenes2 参考数据
echo ""
echo "[INFO] 下载 GreenGenes2 参考数据库..."
source "$CONDA_CMD" activate qiime2 2>/dev/null || conda activate qiime2

python "$SCRIPT_DIR/download_refs.py"

# 4. 验证
echo ""
echo "=========================================="
echo "环境部署完成！"
echo "=========================================="
echo ""
echo "激活环境:"
echo "  conda activate qiime2"
echo ""
echo "快速测试:"
echo "  qiime --help"
echo "  qiime info"
echo ""
echo "查看 GreenGenes2 分类器:"
echo "  ls -lh $PROJECT_DIR/envs/*.qza"
echo ""
echo "开始分析:"
echo "  snakemake -np  # 预览"
echo "  snakemake --cores 4"
echo ""
