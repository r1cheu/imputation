rule fastp_trim:
    conda:
        "../envs/fastp.yml"
    input:
        unpack(get_fastq),
    output:
        multiext("results/trimmed/{sample}", r1=".1.fq.gz", r2=".2.fq.gz"),
    log:
        "logs/fastp/{sample}.log",
    benchmark:
        "benchmarks/fastp/{sample}.tsv"
    threads: 8
    resources:
        mem_mb=4000,
    params:
        extra=config["fastp"]["extra"],
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-w {threads} {params.extra} > {log} 2>&1"
