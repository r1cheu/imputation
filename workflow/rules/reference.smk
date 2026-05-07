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
    resources:
        mem_mb=32000,
        runtime=120,
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
    shell:
        "samtools faidx {input} > {log} 2>&1"


rule panel_index:
    input:
        PANEL_VCF,
    output:
        f"{PANEL_VCF}.tbi",
    log:
        "logs/reference/panel_index.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "tabix -p vcf {input} > {log} 2>&1"


rule split_panel:
    input:
        vcf=PANEL_VCF,
        tbi=f"{PANEL_VCF}.tbi",
    output:
        vcf="results/panel/{chrom}.vcf.gz",
        tbi="results/panel/{chrom}.vcf.gz.tbi",
    log:
        "logs/panel/split_{chrom}.log",
    conda:
        "../envs/bcftools.yaml"
    threads: 2
    resources:
        mem_mb=4000,
        cpus_per_task=2,
    shell:
        "(bcftools view -r {wildcards.chrom} -Oz -o {output.vcf} --threads {threads} {input.vcf} && "
        "tabix -p vcf {output.vcf}) > {log} 2>&1"
