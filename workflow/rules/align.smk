rule minibwa_align_dedup:
    input:
        r1="results/trimmed/{sample}.1.fq.gz",
        r2="results/trimmed/{sample}.2.fq.gz",
        idx=multiext(config["reference"]["fasta"], l2b=".l2b", mbw=".mbw"),
    output:
        multiext("results/dedup/{sample}", bam=".bam", bai=".bam.bai"),
    log:
        "logs/minibwa_dedup/{sample}.log",
    benchmark:
        "benchmarks/minibwa_dedup/{sample}.tsv"
    conda:
        "../envs/align.yml"
    threads: 16
    resources:
        mem_mb=24000,
    params:
        rg=get_read_group,
        idx_prefix=config["reference"]["fasta"],
    shell:
        "(minibwa map -t {threads} -R '{params.rg}' {params.idx_prefix} {input.r1} {input.r2} | "
        "samtools fixmate -@ {threads} -m -u - - | "
        "samtools sort -@ {threads} -u -m 1G - | "
        "samtools markdup -@ {threads} -r --write-index - {output.bam}##idx##{output.bai}) > {log} 2>&1"
