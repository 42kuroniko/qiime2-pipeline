# qiime2-pipeline

Snakemake workflow for 16S amplicon bioinformatics analysis using QIIME2 + GreenGenes2.

## 流程概览

```
fastq files
    │
    ▼
qiime2 import (manifest)
    │
    ▼
DADA2 denoise (QC, trim, deblur)
    │
    ▼
ASV table + Representative sequences
    │
    ├──► Taxonomy (GreenGenes2 Naive Bayes)
    │         │
    │         ▼
    │    taxonomy.tsv
    │
    ├──► Phylogenetic tree (MAFFT + FastTree)
    │         │
    │         ▼
    │    rooted-tree.qza
    │
    ├──► Rarefaction + Alpha/Beta diversity
    │         │
    │         ▼
    │    unweighted-unifrac, weighted-unifrac,
    │    bray-curtis, jaccard
    │    + group significance tests
    │
    ▼
Export for R (phyloseq)
    │
    ├──► feature-table.biom
    ├──► taxonomy.tsv
    ├──► tree.nwk
    ├──► sam.tsv (sample metadata)
    └──► ASV-sequences.fasta
```

## 目录结构

```
qiime2-pipeline/
├── Snakefile              # 主流程定义
├── config.yaml            # 配置文件
├── README.md
├── .gitignore
├── data/
│   ├── manifest.tsv       # 样本列表 (需填写)
│   └── metadata.tsv       # 样本元数据 (需填写)
├── envs/
│   ├── qiime2.yaml        # conda 环境
│   └── gg_2024_09_12_4mer_proxy_tax_classifier.qza   # 下载得到
├── scripts/
│   ├── setup_env.sh       # 一键环境部署
│   ├── download_refs.py    # 下载 GreenGenes2 参考库
│   └── make_manifest.py   # 从 fastq 文件夹自动生成 manifest
├── results/
│   ├── tables/            # ASV 表
│   ├── rep_seqs/          # 代表性序列
│   ├── taxonomy/         # 物种注释
│   ├── tree/              # 系统发育树
│   ├── diversity/         # 多样性分析结果
│   └── export/            # R 导出
└── logs/                  # 运行日志
```

## 快速开始

### 1. 部署环境

```bash
bash scripts/setup_env.sh
```

这会：
- 创建 `qiime2` conda 环境
- 下载 GreenGenes2 2024 参考数据库和分类器

手动激活环境：
```bash
conda activate qiime2
```

### 2. 准备数据文件

**2a. 放置你的 fastq 文件**

把你的 `.fastq.gz` 文件放到 `data/` 或任意目录。

**2b. 生成 manifest.tsv**

```bash
python scripts/make_manifest.py /path/to/your/fastq/folder
```

这会自动扫描文件夹，生成 `data/manifest.tsv`。
检查并确认路径为绝对路径，路径格式正确。

**2c. 填写 metadata.tsv**

复制模板并编辑：
```bash
cp data/metadata_template.tsv data/metadata.tsv
# 用文本编辑器填入你的实验分组信息
```

格式示例 (tab 分隔)：
```
sample-id  treatment  group  day
ctrl_01    control    A     0
ctrl_02    control    A     7
treat_01   treated    B     0
treat_02   treated    B     7
```

### 3. 配置参数

编辑 `config.yaml`：

```yaml
# manifest 和 metadata 路径
manifest: "data/manifest.tsv"
metadata: "data/metadata.tsv"

# DADA2 参数（根据你的测序区域调整）
dada2:
  trim_left_f: 20
  trim_left_r: 20
  trunc_len_f: 0
  trunc_len_r: 0
  max_ee: 2.0
  trunc_q: 2

# 稀有化深度（设为最低样本序列数的 90%）
diversity:
  sampling_depth: 1000
```

> **DADA2 参数参考：**
> - **V4 区** (250bp PE): `trunc_len_f=0, trunc_len_r=0, trim_left_f=17, trim_left_r=21`
> - **V1-V2 区**: `trunc_len_f=240, trunc_len_r=160, trim_left_f=17, trim_left_r=21`
> - 如果 DADA2 报 `trim overhang` 错误，减少 `trunc_len_f/r`

### 4. 运行流程

```bash
# 预览 (不实际执行)
snakemake -np

# 实际运行 (使用 4 核)
snakemake --cores 4

# 后台运行
nohup snakemake --cores 8 > logs/snakemake.log 2>&1 &

# 从特定步骤恢复
snakemake --cores 4 --resume
```

### 5. R 中加载结果

```r
library(phyloseq)
library(biom)

# 加载数据
biom_file <- "results/export/feature-table.biom"
tax_file <- "results/export/taxonomy.tsv"
tree_file <- "results/export/tree.nwk"
sam_file <- "results/export/sam.tsv"

# 读取
ps <- phyloseq(
  otu_table(biom_file, parseTaxonomy = FALSE),
  tax_table(as.matrix(read.table(tax_file, sep="\t", header=TRUE, row.names=1))),
  phy_tree(read_tree(tree_file)),
  sample_data(read.table(sam_file, sep="\t", header=TRUE, row.names=1))
)

# 基本分析
plot_richness(ps, x="treatment", measures=c("Shannon", "observed"))
ord <- ordinate(ps, "PCoA", "unifrac")
plot_ordination(ps, ord, color="treatment")
```

## 常用参数

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `--cores N` | 使用 N 核并行 | 1 |
| `-np` | 预览模式（不执行） | - |
| `-F` | 强制重新运行所有规则 | - |
| `--unlock` | 解锁被中断的运行 | - |
| `-k` | 某个 job 失败时继续运行其他 | - |
| `--latency-wait 60` | 等待文件系统同步 | - |

## 输出文件说明

| 文件 | 说明 |
|------|------|
| `merged-table.qza` | 合并后的 ASV 表 |
| `rarefied-table.qza` | 稀有化后的 ASV 表（用于多样性） |
| `taxonomy.qza/tsv` | 物种注释 |
| `rooted-tree.qza` | 有根系统发育树 |
| `unweighted-unifrac_distance-matrix.qza` | Beta 多样性距离矩阵 |
| `core-metrics-results/` | Alpha/beta 多样性全套结果 |
| `feature-table.biom` | R phyloseq 可直接导入的格式 |

## 依赖

- **conda / mamba** (推荐用 mamba 加速环境创建)
- **qiime2 >= 2024.2**
- **snakemake >= 7.0**
- **GreenGenes2 2024.09 分类器**

## 参考

- QIIME2 文档: https://docs.qiime2.org/
- GreenGenes2: https://docs.qiime2.org/2024.2/tutorials/using-qiime2-files/
- Snakemake: https://snakemake.readthedocs.io/
- phyloseq (R): https://joey711.github.io/phyloseq/
