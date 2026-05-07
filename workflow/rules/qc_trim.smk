rule fastqc_raw:
    input:
        unpack(get_fastq),
    output:
        r1_html="results/qc/fastqc/raw/{sample}_1_fastqc.html",
        r1_zip="results/qc/fastqc/raw/{sample}_1_fastqc.zip",
        r2_html="results/qc/fastqc/raw/{sample}_2_fastqc.html",
        r2_zip="results/qc/fastqc/raw/{sample}_2_fastqc.zip",
    log:
        "logs/fastqc/raw/{sample}.log",
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    resources:
        mem_mb=2000,
        cpus_per_task=2,
        runtime=60,
    params:
        outdir=lambda wc, output: str(Path(output.r1_html).parent),
    shell:
        # rename to {sample}_{1,2} so fastqc output filenames are stable
        """
        mkdir -p {params.outdir}
        tmp=$(mktemp -d)
        ln -sf $(readlink -f {input.r1}) $tmp/{wildcards.sample}_1.fq.gz
        ln -sf $(readlink -f {input.r2}) $tmp/{wildcards.sample}_2.fq.gz
        fastqc -t {threads} -o {params.outdir} $tmp/{wildcards.sample}_1.fq.gz $tmp/{wildcards.sample}_2.fq.gz > {log} 2>&1
        rm -rf $tmp
        """


rule fastp_trim:
    input:
        unpack(get_fastq),
    output:
        r1=temp("results/trimmed/{sample}.1.fq.gz"),
        r2=temp("results/trimmed/{sample}.2.fq.gz"),
        html="results/qc/fastp/{sample}.html",
        json="results/qc/fastp/{sample}.json",
    log:
        "logs/fastp/{sample}.log",
    conda:
        "../envs/fastp.yaml"
    threads: 2
    resources:
        mem_mb=2000,
        cpus_per_task=2,
        runtime=120,
    params:
        extra=config["fastp"]["extra"],
    shell:
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "-h {output.html} -j {output.json} -w {threads} {params.extra} > {log} 2>&1"
