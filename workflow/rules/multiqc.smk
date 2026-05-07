rule multiqc:
    input:
        fastqc=expand(
            "results/qc/fastqc/raw/{sample}_{r}_fastqc.zip",
            sample=SAMPLES,
            r=[1, 2],
        ),
        fastp=expand("results/qc/fastp/{sample}.json", sample=SAMPLES),
        flagstat=expand("results/qc/flagstat/{sample}.flagstat", sample=SAMPLES),
        markdup=expand("results/qc/dedup/{sample}.markdup.stats", sample=SAMPLES),
    output:
        report="results/multiqc/multiqc_report.html",
    log:
        "logs/multiqc.log",
    conda:
        "../envs/multiqc.yaml"
    resources:
        mem_mb=8000,
    params:
        outdir=lambda wc, output: str(Path(output.report).parent),
        scan_dirs="results/qc logs/fastp",
    shell:
        "multiqc --force --outdir {params.outdir} -n multiqc_report.html {params.scan_dirs} > {log} 2>&1"
