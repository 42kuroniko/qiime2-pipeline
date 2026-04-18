# qiime2-pipeline

Dual-end 16S amplicon analysis using QIIME2 + Snakemake + DADA2.

## Study

**PRJNA1309757** - Urothelial Carcinoma and Benign Bladder Tissue Microbiome

72 paired samples from formalin-fixed TURBT bladder specimens. Illumina MiniSeq, 16S rRNA gene amplicon. Comparison of urothelial carcinoma (UCC) vs benign bladder tissue, including BCG-treated samples.

## Pipeline

```
fastq -> import -> DADA2 denoise -> ASV table -> taxonomy
       -> phylogeny (MAFFT + FastTree) -> diversity (alpha/beta)
       -> BIOM export
```

## Setup

```bash
bash scripts/setup_env.sh
mamba env create -f envs/qiime2.yaml
python scripts/download_refs.py
snakemake -n
snakemake --use-conda -j 4
```

## Output

```
results/
  tables/          ASV tables
  taxonomy/        GreenGenes2 taxonomy
  tree/            Phylogenetic tree
  diversity/       Alpha/beta diversity
  export/          BIOM + TSV for R/phyloseq
```

## Citation

MICROBIAL PROFILING OF UROTHELIAL CARCINOMA AND BENIGN BLADDER TISSUE FROM FORMALIN FIXED SPECIMENS. Rush University Medical Center. NCBI BioProject PRJNA1309757, 2025.
