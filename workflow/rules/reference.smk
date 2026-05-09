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
    cache: True
    threads: 1
    resources:
        mem_mb=12000,
    shell:
        "bwa-mem2 index {input} > {log} 2>&1"


rule samtools_faidx:
    input:
        REF_FASTA,
    output:
        f"{REF_FASTA}.fai",
    log:
        "logs/reference/samtools_faidx.log",
    cache: True
    threads: 1
    resources:
        mem_mb=500,
    shell:
        "samtools faidx {input} > {log} 2>&1"


