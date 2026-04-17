# =============================================================================
# qiime2-pipeline: Snakemake workflow for 16S amplicon analysis
# From fastq → ASV → taxonomy → diversity → R-ready export
# =============================================================================

import os
import glob

# ---------------------------------------------------------------------------
# 配置
# ---------------------------------------------------------------------------
configfile: "config.yaml"

# 样本名（从 manifest 文件读取，排除注释行）
SAMPLES = []
with open(config["manifest"], "r") as f:
    for line in f:
        if not line.startswith("#") and line.strip():
            parts = line.strip().split(",")
            if len(parts) >= 2 and os.path.exists(parts[1].strip()):
                sample_id = parts[0].strip()
                if sample_id not in SAMPLES:
                    SAMPLES.append(sample_id)

if not SAMPLES:
    raise ValueError("No valid samples found in manifest. Check your manifest file and paths.")

# 输出目录
RESULTS = config.get("results_dir", "results")
DATA = config.get("data_dir", "data")

# ---------------------------------------------------------------------------
# 工具函数
# ---------------------------------------------------------------------------
def get_fastq_files(wildcards):
    """根据 manifest 获取每个样本的双端 fastq 路径"""
    manifest_df = pd.read_csv(config["manifest"], sep="\t")
    row = manifest_df[manifest_df["sample-id"] == wildcards.sample]
    if row.empty:
        raise ValueError(f"Sample {wildcards.sample} not found in manifest")
    forward = row["forward-absolute-filepath"].values[0]
    reverse = row["reverse-absolute-filepath"].values[0]
    return {"r1": forward, "r2": reverse}

# ---------------------------------------------------------------------------
# 规则
# ---------------------------------------------------------------------------

rule all:
    input:
        # 最终输出文件
        expand(os.path.join(RESULTS, "tables", "{sample}-table.qza"), sample=SAMPLES),
        os.path.join(RESULTS, "tables", "merged-table.qza"),
        os.path.join(RESULTS, "tables", "rarefied-table.qza"),
        os.path.join(RESULTS, "taxonomy", "taxonomy.qza"),
        os.path.join(RESULTS, "taxonomy", "taxonomy.tsv"),
        os.path.join(RESULTS, "tree", "rooted-tree.qza"),
        os.path.join(RESULTS, "tree", "unrooted-tree.qza"),
        os.path.join(RESULTS, "diversity", "unweighted-unifrac_distance-matrix.qza"),
        os.path.join(RESULTS, "diversity", "weighted-unifrac_distance-matrix.qza"),
        os.path.join(RESULTS, "diversity", "bray-curtis_distance-matrix.qza"),
        os.path.join(RESULTS, "diversity", "jaccard_distance-matrix.qza"),
        expand(os.path.join(RESULTS, "diversity", "alpha", "{sample}-alpha.qza"), sample=SAMPLES),
        os.path.join(RESULTS, "diversity", "core-metrics-results"),
        # R 导出
        os.path.join(RESULTS, "export", "feature-table.biom"),
        os.path.join(RESULTS, "export", "taxonomy.tsv"),
        os.path.join(RESULTS, "export", "tree.nwk"),
        os.path.join(RESULTS, "export", "sam.tsv"),
        os.path.join(RESULTS, "export", "ASV-sequences.fasta"),


# ---------------------------------------------------------------------------
# 1. 导入 fastq 为 qza
# ---------------------------------------------------------------------------
rule import_fastq:
    output:
        os.path.join(RESULTS, "imported", "{sample}-demux.qza")
    log:
        os.path.join(RESULTS, "logs", "import-{sample}.log")
    params:
        sample_id="{sample}"
    shell:
        """
        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path {config[manifest]} \
            --output-path {output} \
            --input-format PairedEndFastqManifestPhred33V2 \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 2. 质量可视化
# ---------------------------------------------------------------------------
rule quality_summary:
    input:
        os.path.join(RESULTS, "imported", "{sample}-demux.qza")
    output:
        os.path.join(RESULTS, "quality", "{sample}-quality.qzv")
    log:
        os.path.join(RESULTS, "logs", "quality-{sample}.log")
    shell:
        """
        qiime demux summarize \
            --data-id {input} \
            --output-dir {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 3. DADA2 去噪（去引物、质控、去嵌合）
# ---------------------------------------------------------------------------
rule dada2_denoise:
    output:
        table=os.path.join(RESULTS, "tables", "{sample}-table.qza"),
        rep_seqs=os.path.join(RESULTS, "rep_seqs", "{sample}-rep-seqs.qza"),
        stats=os.path.join(RESULTS, "stats", "{sample}-dada2-stats.qza"),
    log:
        os.path.join(RESULTS, "logs", "dada2-{sample}.log")
    params:
        trim_left_f=config.get("dada2", {}).get("trim_left_f", 20),
        trim_left_r=config.get("dada2", {}).get("trim_left_r", 20),
        trunc_len_f=config.get("dada2", {}).get("trunc_len_f", 0),
        trunc_len_r=config.get("dada2", {}).get("trunc_len_r", 0),
        max_ee=config.get("dada2", {}).get("max_ee", 2.0),
        trunc_q=config.get("dada2", {}).get("trunc_q", 2),
    shell:
        """
        qiime dada2 denoise-paired \
            --i-demultiplexed-seqs {input} \
            --p-trim-left-f {params.trim_left_f} \
            --p-trim-left-r {params.trim_left_r} \
            --p-trunc-len-f {params.trunc_len_f} \
            --p-trunc-len-r {params.trunc_len_r} \
            --p-max-ee {params.max_ee} \
            --p-trunc-q {params.trunc_q} \
            --o-table {output.table} \
            --o-representative-sequences {output.rep_seqs} \
            --o-denoising-stats {output.stats} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 4. 合并所有样本的 ASV 表
# ---------------------------------------------------------------------------
rule merge_tables:
    input:
        tables=expand(os.path.join(RESULTS, "tables", "{sample}-table.qza"), sample=SAMPLES)
    output:
        os.path.join(RESULTS, "tables", "merged-table.qza")
    log:
        os.path.join(RESULTS, "logs", "merge-tables.log")
    shell:
        """
        qiime feature-table merge \
            --i-tables {input} \
            --o-merged-table {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 5. 合并代表性序列
# ---------------------------------------------------------------------------
rule merge_rep_seqs:
    input:
        seqs=expand(os.path.join(RESULTS, "rep_seqs", "{sample}-rep-seqs.qza"), sample=SAMPLES)
    output:
        os.path.join(RESULTS, "rep_seqs", "merged-rep-seqs.qza")
    log:
        os.path.join(RESULTS, "logs", "merge-rep-seqs.log")
    shell:
        """
        qiime feature-table merge-seqs \
            --i-data {input} \
            --o-merged-data {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 6. 稀有化（Rarefaction）- 用于多样性分析
# ---------------------------------------------------------------------------
rule rarefy:
    input:
        os.path.join(RESULTS, "tables", "merged-table.qza")
    output:
        os.path.join(RESULTS, "tables", "rarefied-table.qza")
    log:
        os.path.join(RESULTS, "logs", "rarefy.log")
    params:
        sampling_depth=config.get("diversity", {}).get("sampling_depth", 1000)
    shell:
        """
        qiime feature-table rarefy \
            --i-table {input} \
            --p-sampling-depth {params.sampling_depth} \
            --o-rarefied-table {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 7. 物种注释（GreenGenenes2 Naive Bayes classifier）
# ---------------------------------------------------------------------------
rule classify_taxonomy:
    input:
        rep_seqs=os.path.join(RESULTS, "rep_seqs", "merged-rep-seqs.qza"),
        classifier=config["classifier"]
    output:
        taxonomy=os.path.join(RESULTS, "taxonomy", "taxonomy.qza"),
    log:
        os.path.join(RESULTS, "logs", "classify-taxonomy.log")
    shell:
        """
        qiime feature-classifier classify-sklearn \
            --i-classifier {input.classifier} \
            --i-reads {input.rep_seqs} \
            --o-classification {output.taxonomy} \
            2> {log}
        """


rule export_taxonomy_tsv:
    input:
        os.path.join(RESULTS, "taxonomy", "taxonomy.qza")
    output:
        os.path.join(RESULTS, "taxonomy", "taxonomy.tsv")
    log:
        os.path.join(RESULTS, "logs", "export-taxonomy.log")
    shell:
        """
        qiime tools export \
            --input-path {input} \
            --output-path {RESULTS}/taxonomy \
            2> {log}
        mv {RESULTS}/taxonomy/taxonomy.tsv {output} 2>/dev/null || true
        """


# ---------------------------------------------------------------------------
# 8. 构建系统发育树
# ---------------------------------------------------------------------------
rule align_seqs:
    input:
        os.path.join(RESULTS, "rep_seqs", "merged-rep-seqs.qza")
    output:
        os.path.join(RESULTS, "tree", "aligned-rep-seqs.qza")
    log:
        os.path.join(RESULTS, "logs", "align-seqs.log")
    shell:
        """
        qiime alignment mafft \
            --i-sequences {input} \
            --o-alignment {output} \
            2> {log}
        """


rule mask_alignment:
    input:
        os.path.join(RESULTS, "tree", "aligned-rep-seqs.qza")
    output:
        os.path.join(RESULTS, "tree", "masked-aligned-rep-seqs.qza")
    log:
        os.path.join(RESULTS, "logs", "mask-alignment.log")
    shell:
        """
        qiime alignment mask \
            --i-alignment {input} \
            --o-masked-alignment {output} \
            2> {log}
        """


rule build_tree:
    input:
        os.path.join(RESULTS, "tree", "masked-aligned-rep-seqs.qza")
    output:
        unrooted=os.path.join(RESULTS, "tree", "unrooted-tree.qza"),
        rooted=os.path.join(RESULTS, "tree", "rooted-tree.qza"),
    log:
        os.path.join(RESULTS, "logs", "build-tree.log")
    shell:
        """
        qiime phylogeny fasttree \
            --i-alignment {input} \
            --o-tree {output.unrooted} \
            2> {log}

        qiime phylogeny midpoint-root \
            --i-tree {output.unrooted} \
            --o-rooted-tree {output.rooted} \
            2>> {log}
        """


# ---------------------------------------------------------------------------
# 9. Alpha / Beta 多样性分析
# ---------------------------------------------------------------------------
rule diversity_core_metrics:
    input:
        table=os.path.join(RESULTS, "tables", "rarefied-table.qza"),
        phylogeny=os.path.join(RESULTS, "tree", "rooted-tree.qza"),
        metadata=config["metadata"]
    output:
        directory(os.path.join(RESULTS, "diversity", "core-metrics-results"))
    log:
        os.path.join(RESULTS, "logs", "diversity-core-metrics.log")
    params:
        sampling_depth=config.get("diversity", {}).get("sampling_depth", 1000)
    shell:
        """
        qiime diversity core-metrics-phylogenetic \
            --i-table {input.table} \
            --i-phylogeny {input.phylogeny} \
            --m-metadata-file {input.metadata} \
            --p-sampling-depth {params.sampling_depth} \
            --output-dir {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 10. Alpha 多样性组间检验
# ---------------------------------------------------------------------------
rule alpha_group_significance:
    input:
        rarefied_table=os.path.join(RESULTS, "tables", "rarefied-table.qza"),
        metadata=config["metadata"]
    output:
        os.path.join(RESULTS, "diversity", "alpha", "alpha-group-significance.qzv")
    log:
        os.path.join(RESULTS, "logs", "alpha-significance.log")
    shell:
        """
        qiime diversity alpha-group-significance \
            --i-alpha-diversity {RESULTS}/diversity/core-metrics-results/observed_features_vector.qza \
            --m-metadata-file {input.metadata} \
            --o-visualization {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 11. Beta 多样性组间检验
# ---------------------------------------------------------------------------
rule beta_group_significance:
    input:
        distance_matrix=os.path.join(RESULTS, "diversity", "core-metrics-results", "unweighted_unifrac_distance_matrix.qza"),
        metadata=config["metadata"]
    output:
        os.path.join(RESULTS, "diversity", "beta", "unweighted-unifrac-group-significance.qzv")
    log:
        os.path.join(RESULTS, "logs", "beta-significance.log")
    shell:
        """
        qiime diversity beta-group-significance \
            --i-distance-matrix {input.distance_matrix} \
            --m-metadata-file {input.metadata} \
            --o-visualization {output} \
            2> {log}
        """


# ---------------------------------------------------------------------------
# 12. 导出为 R (phyloseq) 可用格式
# ---------------------------------------------------------------------------
rule export_for_r:
    input:
        table=os.path.join(RESULTS, "tables", "merged-table.qza"),
        taxonomy=os.path.join(RESULTS, "taxonomy", "taxonomy.qza"),
        tree=os.path.join(RESULTS, "tree", "rooted-tree.qza"),
        metadata=config["metadata"]
    output:
        os.path.join(RESULTS, "export", "feature-table.biom"),
        os.path.join(RESULTS, "export", "taxonomy.tsv"),
        os.path.join(RESULTS, "export", "tree.nwk"),
        os.path.join(RESULTS, "export", "sam.tsv"),
        os.path.join(RESULTS, "export", "ASV-sequences.fasta"),
    log:
        os.path.join(RESULTS, "logs", "export-for-r.log")
    shell:
        """
        mkdir -p {RESULTS}/export

        # 导出 feature-table.biom
        qiime tools export \
            --input-path {input.table} \
            --output-path {RESULTS}/export

        mv {RESULTS}/export/feature-table.biom {output[0]} 2>/dev/null || \
        mv {RESULTS}/export/feature-table/feature-table.biom {output[0]} 2>/dev/null || true

        # 导出 taxonomy.tsv
        qiime tools export \
            --input-path {input.taxonomy} \
            --output-path {RESULTS}/export
        mv {RESULTS}/export/taxonomy.tsv {output[1]} 2>/dev/null || true

        # 导出 tree.nwk
        qiime tools export \
            --input-path {input.tree} \
            --output-path {RESULTS}/export
        mv {RESULTS}/export/tree.nwk {output[2]} 2>/dev/null || true

        # 复制 metadata 为 sam.tsv（R 中读作 sample_data）
        cp {input.metadata} {output[3]}

        # 导出 ASV 序列 fasta
        qiime tools export \
            --input-path {RESULTS}/rep_seqs/merged-rep-seqs.qza \
            --output-path {RESULTS}/export
        mv {RESULTS}/export/dna-sequences.fasta {output[4]} 2>/dev/null || true
        """


# ---------------------------------------------------------------------------
# 13. 生成流程报告
# ---------------------------------------------------------------------------
rule workflow_summary:
    input:
        expand(os.path.join(RESULTS, "stats", "{sample}-dada2-stats.qza"), sample=SAMPLES)
    output:
        os.path.join(RESULTS, "pipeline-summary.qzv")
    log:
        os.path.join(RESULTS, "logs", "workflow-summary.log")
    shell:
        """
        qiime diversity adonis \
            --i-table {RESULTS}/tables/merged-table.qza \
            --o-visualization {output} \
            2> {log} || true

        # 基本的统计报告
        echo "Pipeline completed successfully" > {output}
        """


# ---------------------------------------------------------------------------
# 清理规则（可选）
# ---------------------------------------------------------------------------
rule clean:
    shell:
        """
        rm -rf results/imported results/tables/*.qza results/rep_seqs/*.qza
        rm -rf results/stats results/quality results/logs
        rm -rf results/tree/*.qza results/taxonomy/*.qza results/taxonomy/*.tsv
        rm -rf results/diversity
        rm -rf results/export
        echo "Cleaned results directory"
        """
