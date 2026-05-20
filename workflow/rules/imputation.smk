rule compute_gl:
    input:
        bam="results/dedup/{sample}.bam",
        bai="results/dedup/{sample}.bam.bai",
        ref=REF_FASTA,
        fai=f"{REF_FASTA}.fai",
        panel=PANEL_FULL,
        panel_csi=f"{PANEL_FULL}.csi",
        sites_tsv=PANEL_SITES_TSV,
        sites_tbi=f"{PANEL_SITES_TSV}.tbi",
    output:
        bcf=temp("results/gl/{sample}/{chrom}.bcf"),
        csi=temp("results/gl/{sample}/{chrom}.bcf.csi"),
    log:
        "logs/compute_gl/{sample}_{chrom}.log",
    threads: 2
    resources:
        mem_mb=4000,
    shell:
        "(bcftools mpileup -f {input.ref} -I -E -a 'FORMAT/DP' "
        "-T {input.panel} -r {wildcards.chrom} {input.bam} -Ou | "
        "bcftools call -Aim -C alleles -T {input.sites_tsv} -Ob -o {output.bcf} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


rule concat_gl:
    input:
        bcfs=expand("results/gl/{{sample}}/{chrom}.bcf", chrom=CHROMS),
        csis=expand("results/gl/{{sample}}/{chrom}.bcf.csi", chrom=CHROMS),
    output:
        bcf="results/gl/{sample}.bcf",
        csi="results/gl/{sample}.bcf.csi",
    log:
        "logs/concat_gl/{sample}.log",
    threads: 2
    resources:
        mem_mb=2000,
    shell:
        "(bcftools concat --threads {threads} -Ob -o {output.bcf} {input.bcfs} && "
        "bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1"


rule merge_gl:
    input:
        bcfs=lambda wc: [f"results/gl/{s}.bcf" for s in SAMPLES],
        csis=lambda wc: [f"results/gl/{s}.bcf.csi" for s in SAMPLES],
    output:
        bcf="results/gl_merged/{chrom}.bcf",
        csi="results/gl_merged/{chrom}.bcf.csi",
    log:
        "logs/merge_gl/{chrom}.log",
    threads: 4
    resources:
        mem_mb=4000,
    shell:
        "(bcftools merge -m none -r {wildcards.chrom} --threads {threads} "
        "-Ob -o {output.bcf} {input.bcfs} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


checkpoint glimpse_chunk:
    input:
        bcf=PANEL_FULL,
        csi=f"{PANEL_FULL}.csi",
        gmap=get_map,
    output:
        "results/chunks/{chrom}.txt",
    log:
        "logs/glimpse_chunk/{chrom}.log",
    threads: 1
    resources:
        mem_mb=2000,
    params:
        window_size=config["glimpse_chunk"]["window_size"],
        buffer_size=config["glimpse_chunk"]["buffer_size"],
        extra=config["glimpse_chunk"]["extra"],
    shell:
        "GLIMPSE_chunk --input {input.bcf} --map {input.gmap} --region {wildcards.chrom} "
        "--window-size {params.window_size} --buffer-size {params.buffer_size} "
        "--thread {threads} --output {output} {params.extra} > {log} 2>&1"


rule glimpse_phase:
    input:
        gl="results/gl_merged/{chrom}.bcf",
        gl_csi="results/gl_merged/{chrom}.bcf.csi",
        ref=PANEL_FULL,
        ref_csi=f"{PANEL_FULL}.csi",
        gmap=get_map,
    output:
        bcf=temp("results/phased/{chrom}/chunk_{idx}.bcf"),
        csi=temp("results/phased/{chrom}/chunk_{idx}.bcf.csi"),
    log:
        "logs/glimpse_phase/{chrom}_chunk_{idx}.log",
    threads: 16
    resources:
        mem_mb=16000,
    params:
        input_region=lambda wc: get_chunk_region(wc, "input_region"),
        output_region=lambda wc: get_chunk_region(wc, "output_region"),
    shell:
        "(GLIMPSE_phase --input {input.gl} --reference {input.ref} --map {input.gmap} "
        "--input-region {params.input_region} --output-region {params.output_region} "
        "--thread {threads} --output {output.bcf} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


rule glimpse_ligate:
    input:
        bcfs=phased_chunks,
        csis=phased_chunks_idx,
    output:
        bcf="results/imputed/{chrom}.bcf",
        csi="results/imputed/{chrom}.bcf.csi",
    log:
        "logs/glimpse_ligate/{chrom}.log",
    threads: 4
    resources:
        mem_mb=8000,
    params:
        listfile="results/phased/{chrom}/ligate.list",
    shell:
        """
        printf '%s\\n' {input.bcfs} > {params.listfile}
        (GLIMPSE_ligate --input {params.listfile} --output {output.bcf} --thread {threads} && \
         bcftools index -f {output.bcf} --threads {threads}) > {log} 2>&1
        """


rule concat_imputed:
    input:
        bcfs=expand("results/imputed/{chrom}.bcf", chrom=CHROMS),
        csis=expand("results/imputed/{chrom}.bcf.csi", chrom=CHROMS),
    output:
        bcf="results/imputed/all.bcf",
        csi="results/imputed/all.bcf.csi",
    log:
        "logs/concat_imputed/all.log",
    threads: 4
    resources:
        mem_mb=8000,
    shell:
        "(bcftools concat --threads {threads} -Ob -o {output.bcf} {input.bcfs} && "
        "bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1"
