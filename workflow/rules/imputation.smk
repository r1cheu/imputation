checkpoint glimpse2_chunk:
    input:
        bcf=PANEL_SITES,
        csi=f"{PANEL_SITES}.csi",
        gmap=get_map,
    output:
        "results/chunks/{chrom}.txt",
    log:
        "logs/glimpse2_chunk/{chrom}.log",
    threads: 1
    resources:
        mem_mb=2000,
    params:
        chrom=lambda wc: wc.chrom,
        window_mb=config["glimpse2_chunk"]["window_mb"],
        buffer_mb=config["glimpse2_chunk"]["buffer_mb"],
        extra=config["glimpse2_chunk"]["extra"],
    shell:
        "GLIMPSE2_chunk --input {input.bcf} --map {input.gmap} --region {params.chrom} "
        "--window-mb {params.window_mb} --buffer-mb {params.buffer_mb} --sequential "
        "--threads {threads} --output {output} {params.extra} > {log} 2>&1"


rule glimpse2_split_reference:
    input:
        bcf=PANEL_FULL,
        csi=f"{PANEL_FULL}.csi",
        gmap=get_map,
        chunks="results/chunks/{chrom}.txt",
    output:
        bin="results/refbin/{chrom}/chunk_{idx}.bin",
    log:
        "logs/glimpse2_split_reference/{chrom}_chunk_{idx}.log",
    cache: True
    threads: 2
    resources:
        mem_mb=8000,
    params:
        prefix=lambda wc: f"results/refbin/{wc.chrom}/chunk_{wc.idx}",
        input_region=lambda wc: get_chunk_region(wc, "input_region"),
        output_region=lambda wc: get_chunk_region(wc, "output_region"),
    shell:
        # GLIMPSE2_split_reference appends "_<chr>_<start>_<end>.bin" to --output;
        # rename the single produced file to the canonical {prefix}.bin.
        "GLIMPSE2_split_reference --reference {input.bcf} --map {input.gmap} "
        "--input-region {params.input_region} --output-region {params.output_region} "
        "--threads {threads} --output {params.prefix} > {log} 2>&1 && "
        "mv {params.prefix}_*.bin {output.bin}"


rule make_bam_list:
    input:
        bams=all_sample_bams(),
        bais=all_sample_bais(),
    output:
        "results/bam_list.txt",
    log:
        "logs/make_bam_list.log",
    resources:
        mem_mb=1000,
    shell:
        "(for f in {input.bams}; do readlink -f $f; done > {output}) 2> {log}"


rule glimpse2_phase:
    input:
        bam_list="results/bam_list.txt",
        bams=all_sample_bams(),
        bais=all_sample_bais(),
        ref_bin="results/refbin/{chrom}/chunk_{idx}.bin",
    output:
        bcf="results/phased/{chrom}/chunk_{idx}.bcf",
        csi="results/phased/{chrom}/chunk_{idx}.bcf.csi",
    log:
        "logs/glimpse2_phase/{chrom}_chunk_{idx}.log",
    threads: 8
    resources:
        mem_mb=32000,
    shell:
        "(GLIMPSE2_phase --bam-list {input.bam_list} --reference {input.ref_bin} "
        "--threads {threads} --output {output.bcf} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


rule glimpse2_ligate:
    input:
        bcfs=phased_chunks,
        csis=phased_chunks_idx,
    output:
        bcf="results/imputed/{chrom}.bcf",
        csi="results/imputed/{chrom}.bcf.csi",
    log:
        "logs/glimpse2_ligate/{chrom}.log",
    threads: 4
    resources:
        mem_mb=8000,
    params:
        listfile="results/phased/{chrom}/ligate.list",
    shell:
        """
        printf '%s\n' {input.bcfs} > {params.listfile}
        (GLIMPSE2_ligate --input {params.listfile} --output {output.bcf} --threads {threads} && \
         bcftools index -f {output.bcf} --threads {threads}) > {log} 2>&1
        """
