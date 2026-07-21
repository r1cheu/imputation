rule bwa_mem2_index:
    conda:
        "../envs/bwa-mem2.yml"
    input:
        config["reference"]["fasta"],
    output:
        multiext(
            config["reference"]["fasta"],
            idx0123=".0123",
            amb=".amb",
            ann=".ann",
            bwt=".bwt.2bit.64",
            pac=".pac",
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
    conda:
        "../envs/samtools.yml"
    input:
        fasta=config["reference"]["fasta"],
    output:
        fai=config["reference"]["fasta"] + ".fai",
    log:
        "logs/reference/samtools_faidx.log",
    cache: True
    threads: 1
    resources:
        mem_mb=500,
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"
