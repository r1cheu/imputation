rule bwa_align_dedup:
    conda:
        "../envs/align.yml"
    input:
        reads=multiext("results/trimmed/{sample}", r1=".1.fq.gz", r2=".2.fq.gz"),
        idx=multiext(
            config["reference"]["fasta"],
            idx0123=".0123",
            amb=".amb",
            ann=".ann",
            bwt=".bwt.2bit.64",
            pac=".pac",
        ),
    output:
        multiext("results/dedup/{sample}", bam=".bam", bai=".bam.bai"),
    log:
        "logs/bwa_dedup/{sample}.log",
    benchmark:
        "benchmarks/bwa_dedup/{sample}.tsv"
    threads: 16
    resources:
        mem_mb=24000,
    params:
        rg=get_read_group,
        idx_prefix=config["reference"]["fasta"],
    shell:
        "(bwa-mem2 mem -t {threads} -R '{params.rg}' {params.idx_prefix} {input.reads.r1} {input.reads.r2} | "
        "samtools fixmate -@ {threads} -m -u - - | "
        "samtools sort -@ {threads} -u -m 1G - | "
        "samtools markdup -@ {threads} -r --write-index - {output.bam}##idx##{output.bai}) > {log} 2>&1"
