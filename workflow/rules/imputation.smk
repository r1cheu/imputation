rule compute_gl:
    conda:
        "../envs/bcftools.yml"
    input:
        multiext("results/dedup/{sample}", bam=".bam", bai=".bam.bai"),
        multiext(PANEL_PREFIX, panel_vcf=".vcf.gz", panel_csi=".vcf.gz.csi"),
        multiext(SITES_PREFIX, sites_tsv=".tsv.gz", sites_tbi=".tsv.gz.tbi"),
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fasta"] + ".fai",
    output:
        multiext("results/gl/{sample}/{chrom}", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/compute_gl/{sample}_{chrom}.log",
    threads: 1
    resources:
        mem_mb=4000,
    shell:
        "(bcftools mpileup -f {input.ref} -I -E -a 'FORMAT/DP' "
        "-T {input.panel_vcf} -r {wildcards.chrom} {input.bam} -Ou | "
        "bcftools call -Aim -C alleles -T {input.sites_tsv} -Ob -o {output.bcf} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


rule concat_gl:
    conda:
        "../envs/bcftools.yml"
    input:
        bcfs=expand("results/gl/{{sample}}/{chrom}.bcf", chrom=CHROMS),
        csis=expand("results/gl/{{sample}}/{chrom}.bcf.csi", chrom=CHROMS),
    output:
        multiext("results/gl/{sample}", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/concat_gl/{sample}.log",
    threads: 1
    resources:
        mem_mb=2000,
    shell:
        "(bcftools concat --threads {threads} -Ob -o {output.bcf} {input.bcfs} && "
        "bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1"


rule merge_gl_batch:
    conda:
        "../envs/bcftools.yml"
    input:
        bcfs=merge_gl_batch_bcfs,
        csis=merge_gl_batch_csis,
    output:
        temp(
            multiext("results/gl_merged_batches/{chrom}/batch_{batch}", bcf=".bcf", csi=".bcf.csi")
        ),
    log:
        "logs/merge_gl_batch/{chrom}_batch_{batch}.log",
    threads: 2
    resources:
        mem_mb=4000,
    params:
        listfile="results/gl_merged_batches/{chrom}/batch_{batch}.list",
    shell:
        """
        printf '%s\\n' {input.bcfs} > {params.listfile}
        (bcftools merge -m none -r {wildcards.chrom} --threads {threads} \
            -Ob -o {output.bcf} -l {params.listfile} && \
         bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1
        """


rule merge_gl:
    conda:
        "../envs/bcftools.yml"
    input:
        bcfs=expand("results/gl_merged_batches/{{chrom}}/batch_{batch}.bcf", batch=MERGE_GL_BATCHES),
        csis=expand("results/gl_merged_batches/{{chrom}}/batch_{batch}.bcf.csi", batch=MERGE_GL_BATCHES),
    output:
        multiext("results/gl_merged/{chrom}", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/merge_gl/{chrom}.log",
    threads: 8
    resources:
        mem_mb=4000,
    params:
        listfile="results/gl_merged/{chrom}.list",
    shell:
        """
        printf '%s\\n' {input.bcfs} > {params.listfile}
        (bcftools merge -m none -r {wildcards.chrom} --threads {threads} \
            -Ob -o {output.bcf} -l {params.listfile} && \
         bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1
        """


checkpoint glimpse_chunk:
    input:
        multiext(PANEL_PREFIX, panel_vcf=".vcf.gz", panel_csi=".vcf.gz.csi"),
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
        "GLIMPSE_chunk --input {input.panel_vcf} --region {wildcards.chrom} "
        "--window-size {params.window_size} --buffer-size {params.buffer_size} "
        "--thread {threads} --output {output} {params.extra} > {log} 2>&1"


rule glimpse_phase:
    input:
        multiext("results/gl_merged/{chrom}", gl_bcf=".bcf", gl_csi=".bcf.csi"),
        multiext(PANEL_PREFIX, ref_vcf=".vcf.gz", ref_csi=".vcf.gz.csi"),
        gmap=get_map,
    output:
        temp(
            multiext("results/phased/{chrom}/chunk_{idx}", bcf=".bcf", csi=".bcf.csi")
        ),
    log:
        "logs/glimpse_phase/{chrom}_chunk_{idx}.log",
    threads: 24
    resources:
        mem_mb=16000,
    params:
        input_region=lambda wc: get_chunk_region(wc, "input_region"),
        output_region=lambda wc: get_chunk_region(wc, "output_region"),
    shell:
        "(GLIMPSE_phase --input {input.gl_bcf} --reference {input.ref_vcf} --map {input.gmap} "
        "--input-region {params.input_region} --output-region {params.output_region} "
        "--thread {threads} --output {output.bcf} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


rule glimpse_ligate:
    input:
        bcfs=phased_chunks,
        csis=phased_chunks_idx,
    output:
        multiext("results/imputed/{chrom}", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/glimpse_ligate/{chrom}.log",
    threads: 8
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
    conda:
        "../envs/bcftools.yml"
    input:
        bcfs=expand("results/imputed/{chrom}.bcf", chrom=CHROMS),
        csis=expand("results/imputed/{chrom}.bcf.csi", chrom=CHROMS),
    output:
        multiext("results/imputed/all", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/concat_imputed/all.log",
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        "(bcftools concat --threads {threads} -Ob -o {output.bcf} {input.bcfs} && "
        "bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1"
