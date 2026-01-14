"""
Cell Ranger Alignment and Counting
===================================

Rules for running Cell Ranger count on 10x Genomics data.
"""

rule cellranger_count:
    """
    Run Cell Ranger count to align reads and generate count matrices
    """
    input:
        fastqs = lambda wildcards: get_sample_fastqs(wildcards)
    output:
        h5 = "results/cellranger/{sample}/outs/filtered_feature_bc_matrix.h5",
        bam = "results/cellranger/{sample}/outs/possorted_genome_bam.bam",
        summary = "results/cellranger/{sample}/outs/web_summary.html",
        metrics = "results/cellranger/{sample}/outs/metrics_summary.csv"
    params:
        reference = config["cellranger_reference"],
        chemistry = config["cellranger"]["chemistry"],
        expect_cells = config["cellranger"]["expect_cells"],
        sample_id = "{sample}",
        outdir = "results/cellranger"
    threads: config["resources"]["cellranger"]["threads"]
    resources:
        mem_gb = config["resources"]["cellranger"]["memory_gb"]
    log:
        "logs/cellranger/{sample}_count.log"
    shell:
        """
        cellranger count \
            --id={params.sample_id} \
            --transcriptome={params.reference} \
            --fastqs={input.fastqs} \
            --sample={params.sample_id} \
            --chemistry={params.chemistry} \
            --expect-cells={params.expect_cells} \
            --localcores={threads} \
            --localmem={resources.mem_gb} \
            2>&1 | tee {log}
        
        # Move output to results directory
        mv {params.sample_id} {params.outdir}/
        """


rule aggregate_cellranger_metrics:
    """
    Aggregate Cell Ranger metrics across all samples
    """
    input:
        expand("results/cellranger/{sample}/outs/metrics_summary.csv",
               sample=SAMPLES)
    output:
        "results/cellranger/aggregated_metrics.csv"
    log:
        "logs/cellranger/aggregate_metrics.log"
    script:
        "../scripts/aggregate_cellranger_metrics.py"
