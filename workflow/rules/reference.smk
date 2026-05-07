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


rule panel_convert:
    input:
        PANEL_TEMPLATE,
    output:
        bcf="results/panel/{chrom}.bcf",
        csi="results/panel/{chrom}.bcf.csi",
    log:
        "logs/panel/convert_{chrom}.log",
    conda:
        "../envs/bcftools.yaml"
    cache: True
    threads: 2
    resources:
        mem_mb=2000,
        cpus_per_task=2,
    shell:
        # convert user-supplied per-chrom .vcf.gz to .bcf once; downstream consumes BCF only
        "(bcftools view -Ob -o {output.bcf} --threads {threads} {input} && "
        "bcftools index -f {output.bcf} --threads {threads}) > {log} 2>&1"


rule panel_sites:
    input:
        bcf="results/panel/{chrom}.bcf",
        csi="results/panel/{chrom}.bcf.csi",
    output:
        bcf="results/panel/{chrom}.sites.bcf",
        csi="results/panel/{chrom}.sites.bcf.csi",
    log:
        "logs/panel/sites_{chrom}.log",
    conda:
        "../envs/bcftools.yaml"
    cache: True
    threads: 2
    resources:
        mem_mb=2000,
        cpus_per_task=2,
    shell:
        # sites-only BCF for GLIMPSE2_chunk (no GT — much faster to load)
        "(bcftools view -G -Ob -o {output.bcf} --threads {threads} {input.bcf} && "
        "bcftools index -f {output.bcf} --threads {threads}) > {log} 2>&1"
