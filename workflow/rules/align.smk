rule bwa_align_dedup:
    input:
        reads=get_trimmed_reads,
        idx=multiext(
            REF_FASTA,
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    output:
        bam="results/dedup/{sample}.bam",
        bai="results/dedup/{sample}.bam.bai",
    log:
        "logs/bwa_dedup/{sample}.log",
    benchmark:
        "benchmarks/bwa_dedup/{sample}.tsv"
    threads: 16
    resources:
        mem_mb=24000,
    params:
        rg=get_read_group,
        idx_prefix=lambda wc, input: input.idx[0].rsplit(".0123", 1)[0],
    shell:
        "(bwa-mem2 mem -t {threads} -R '{params.rg}' {params.idx_prefix} {input.reads} | "
        "samtools fixmate -@ {threads} -m -u - - | "
        "samtools sort -@ {threads} -u -m 1G - | "
        "samtools markdup -@ {threads} -r --write-index - {output.bam}) > {log} 2>&1"
