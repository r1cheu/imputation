rule bwa_mem2_index:
    input:
        REF_FASTA,
    output:
        multiext(
            REF_FASTA,
            ".0123",
            ".amb",
            ".ann",
            ".bwt.2bit.64",
            ".pac",
        ),
    log:
        "logs/reference/bwa_mem2_index.log",
    conda:
        "../envs/bwa-mem2.yaml"
    cache: True
    threads: 8
    resources:
        mem_mb=32000,
        cpus_per_task=8,
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"


rule samtools_faidx:
    input:
        REF_FASTA,
    output:
        f"{REF_FASTA}.fai",
    log:
        "logs/reference/samtools_faidx.log",
    conda:
        "../envs/bwa-mem2.yaml"
    cache: True
    resources:
        mem_mb=2000,
    shell:
        "samtools faidx {input} > {log} 2>&1"


