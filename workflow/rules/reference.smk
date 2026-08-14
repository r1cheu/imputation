rule minibwa_index:
    input:
        config["reference"]["fasta"],
    output:
        multiext(config["reference"]["fasta"], l2b=".l2b", mbw=".mbw"),
    log:
        "logs/reference/minibwa_index.log",
    cache: True
    conda:
        "../envs/minibwa.yml"
    threads: 8
    resources:
        mem_mb=12000,
    shell:
        "minibwa index -t {threads} {input} > {log} 2>&1"


rule samtools_faidx:
    input:
        fasta=config["reference"]["fasta"],
    output:
        fai=config["reference"]["fasta"] + ".fai",
    log:
        "logs/reference/samtools_faidx.log",
    cache: True
    conda:
        "../envs/samtools.yml"
    threads: 1
    resources:
        mem_mb=500,
    shell:
        "samtools faidx {input.fasta} > {log} 2>&1"
