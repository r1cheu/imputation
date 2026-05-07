checkpoint glimpse2_chunk:
    input:
        vcf="results/panel/{chrom}.vcf.gz",
        tbi="results/panel/{chrom}.vcf.gz.tbi",
        gmap=get_map,
    output:
        "results/chunks/{chrom}.txt",
    log:
        "logs/glimpse2_chunk/{chrom}.log",
    conda:
        "../envs/glimpse2.yaml"
    threads: 2
    resources:
        mem_mb=4000,
        cpus_per_task=2,
    params:
        window_mb=config["glimpse2_chunk"]["window_mb"],
        buffer_mb=config["glimpse2_chunk"]["buffer_mb"],
        extra=config["glimpse2_chunk"]["extra"],
    shell:
        "GLIMPSE2_chunk --input {input.vcf} --map {input.gmap} --region {wildcards.chrom} "
        "--window-mb {params.window_mb} --buffer-mb {params.buffer_mb} --sequential "
        "--threads {threads} --output {output} {params.extra} > {log} 2>&1"


rule glimpse2_split_reference:
    input:
        vcf="results/panel/{chrom}.vcf.gz",
        tbi="results/panel/{chrom}.vcf.gz.tbi",
        gmap=get_map,
        chunks="results/chunks/{chrom}.txt",
    output:
        bin="results/refbin/{chrom}/chunk_{idx}.bin",
    log:
        "logs/glimpse2_split_reference/{chrom}_chunk_{idx}.log",
    conda:
        "../envs/glimpse2.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        cpus_per_task=2,
        runtime=120,
    params:
        prefix=lambda wc: f"results/refbin/{wc.chrom}/chunk_{wc.idx}",
    shell:
        # read the matching chunk row by id (column 1) and pass --input-region / --output-region
        """
        line=$(awk -v i={wildcards.idx} '$1==i' {input.chunks})
        if [ -z "$line" ]; then echo "chunk {wildcards.idx} not found in {input.chunks}" >&2; exit 1; fi
        in_region=$(echo "$line" | awk '{{print $3}}')
        out_region=$(echo "$line" | awk '{{print $4}}')
        GLIMPSE2_split_reference --reference {input.vcf} --map {input.gmap} \
            --input-region $in_region --output-region $out_region \
            --threads {threads} --output {params.prefix} > {log} 2>&1
        """


rule make_bam_list:
    input:
        bams=all_sample_bams(),
        bais=all_sample_bais(),
    output:
        "results/bam_list.txt",
    log:
        "logs/make_bam_list.log",
    conda:
        "../envs/bwa-mem2.yaml"
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
    conda:
        "../envs/glimpse2.yaml"
    threads: config["glimpse2_phase"]["threads"]
    resources:
        mem_mb=16000,
        cpus_per_task=config["glimpse2_phase"]["threads"],
        runtime=720,
    shell:
        "(GLIMPSE2_phase --bam-list {input.bam_list} --reference {input.ref_bin} "
        "--threads {threads} --output {output.bcf} && "
        "bcftools index -f {output.bcf}) > {log} 2>&1"


rule glimpse2_ligate:
    input:
        bcfs=phased_chunks,
        csis=phased_chunks_idx,
    output:
        vcf="results/imputed/{chrom}.vcf.gz",
        tbi="results/imputed/{chrom}.vcf.gz.tbi",
    log:
        "logs/glimpse2_ligate/{chrom}.log",
    conda:
        "../envs/glimpse2.yaml"
    threads: 2
    resources:
        mem_mb=8000,
        cpus_per_task=2,
        runtime=240,
    params:
        listfile="results/phased/{chrom}/ligate.list",
    shell:
        """
        printf '%s\n' {input.bcfs} > {params.listfile}
        (GLIMPSE2_ligate --input {params.listfile} --output {output.vcf} --threads {threads} && \
         tabix -p vcf {output.vcf}) > {log} 2>&1
        """
