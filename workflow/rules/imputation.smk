localrules:
    bam_list,


rule glimpse2_split_reference:
    input:
        panel=PANEL,
        panel_csi=PANEL + ".csi",
        gmap=get_map,
        fai=config["reference"]["fasta"] + ".fai",
    output:
        bin="results/split_ref/{chrom}.bin",
    log:
        "logs/glimpse2_split_reference/{chrom}.log",
    conda:
        "../envs/glimpse2.yml"
    threads: 4
    resources:
        mem_mb=8000,
    params:
        prefix="results/split_ref/{chrom}.tmp",
        region=chrom_region,
    shell:
        """
        (GLIMPSE2_split_reference --reference {input.panel} --map {input.gmap} \
            --input-region {params.region} --output-region {params.region} \
            --threads {threads} --output {params.prefix} \
            && mv {params.prefix}_*.bin {output.bin}) >{log} 2>&1
        """


rule bam_list:
    input:
        bams=batch_bams,
    output:
        "results/bam_lists/batch_{batch}.txt",
    params:
        samples=batch_samples,
    run:
        with open(output[0], "w") as f:
            for bam, sample in zip(input.bams, params.samples):
                f.write(f"{bam} {sample}\n")


rule glimpse2_phase:
    input:
        bams=batch_bams,
        bais=batch_bais,
        bam_list="results/bam_lists/batch_{batch}.txt",
        ref="results/split_ref/{chrom}.bin",
    output:
        temp(
            multiext(
                "results/phased/{chrom}/batch_{batch}", bcf=".bcf", csi=".bcf.csi"
            )
        ),
    log:
        "logs/glimpse2_phase/{chrom}_batch_{batch}.log",
    benchmark:
        "benchmarks/glimpse2_phase/{chrom}_batch_{batch}.tsv"
    conda:
        "../envs/glimpse2.yml"
    threads: 24
    resources:
        mem_mb=32000,
    shell:
        "(GLIMPSE2_phase --bam-list {input.bam_list} --reference {input.ref} "
        "--threads {threads} --output {output.bcf} && "
        "bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1"


rule merge_imputed:
    input:
        bcfs=phased_batches,
        csis=phased_batches_idx,
    output:
        multiext("results/imputed/{chrom}", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/merge_imputed/{chrom}.log",
    conda:
        "../envs/bcftools.yml"
    threads: 8
    resources:
        mem_mb=8000,
    params:
        listfile="results/phased/{chrom}/merge.list",
    shell:
        """
        printf '%s\\n' {input.bcfs} >{params.listfile}
        (bcftools merge --force-single -m none --threads {threads} \
            -Ob -o {output.bcf} -l {params.listfile} \
            && bcftools index -f --threads {threads} {output.bcf}) >{log} 2>&1
        """


rule concat_imputed:
    input:
        bcfs=expand("results/imputed/{chrom}.bcf", chrom=CHROMS),
        csis=expand("results/imputed/{chrom}.bcf.csi", chrom=CHROMS),
    output:
        multiext("results/imputed/all", bcf=".bcf", csi=".bcf.csi"),
    log:
        "logs/concat_imputed/all.log",
    conda:
        "../envs/bcftools.yml"
    threads: 8
    resources:
        mem_mb=8000,
    shell:
        "(bcftools concat --threads {threads} -Ob -o {output.bcf} {input.bcfs} && "
        "bcftools index -f --threads {threads} {output.bcf}) > {log} 2>&1"
