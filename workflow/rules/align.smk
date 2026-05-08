rule bwa_mem2_mem:
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
        bam=temp("results/mapped/{sample}.sorted.bam"),
    log:
        "logs/bwa_mem2/{sample}.log",
    conda:
        "../envs/bwa-mem2.yaml"
    threads: 4
    resources:
        mem_mb=8000,
        cpus_per_task=4,
    params:
        rg=get_read_group,
        idx_prefix=lambda wc, input: input.idx[0].rsplit(".0123", 1)[0],
    shell:
        "(bwa-mem2 mem -t {threads} -R '{params.rg}' {params.idx_prefix} {input.reads} | "
        "samtools sort -@ {threads} -o {output.bam} -) > {log} 2>&1"


rule mark_duplicates:
    input:
        "results/mapped/{sample}.sorted.bam",
    output:
        bam="results/dedup/{sample}.bam",
    log:
        "logs/markdup/{sample}.log",
    conda:
        "../envs/bwa-mem2.yaml"
    group:
        "post_align"
    threads: 2
    resources:
        mem_mb=4000,
        cpus_per_task=2,
    shell:
        # samtools markdup needs name-sorted -> fixmate -> coord-sorted -> markdup
        "(samtools collate -@ {threads} -O -u {input} | "
        "samtools fixmate -@ {threads} -m -u - - | "
        "samtools sort -@ {threads} -u - | "
        "samtools markdup -@ {threads} - {output.bam}) > {log} 2>&1"


rule index_bam:
    input:
        "results/dedup/{sample}.bam",
    output:
        "results/dedup/{sample}.bam.bai",
    log:
        "logs/index_bam/{sample}.log",
    conda:
        "../envs/bwa-mem2.yaml"
    group:
        "post_align"
    resources:
        mem_mb=1000,
    shell:
        "samtools index {input} > {log} 2>&1"
