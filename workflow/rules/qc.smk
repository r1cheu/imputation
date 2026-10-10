localrules:
    filter_samples,


rule bam_qc:
    input:
        multiext("results/dedup/{sample}", bam=".bam", bai=".bam.bai"),
    output:
        stats="results/qc/{sample}.stats",
    log:
        "logs/bam_qc/{sample}.log",
    conda:
        "../envs/samtools.yml"
    threads: 2
    resources:
        mem_mb=2000,
    shell:
        "samtools stats -@ {threads} {input.bam} > {output.stats} 2> {log}"


checkpoint filter_samples:
    input:
        stats=expand("results/qc/{sample}.stats", sample=SAMPLES),
        fai=config["reference"]["fasta"] + ".fai",
    output:
        qc="results/qc/sample_qc.tsv",
        passed="results/qc/samples_pass.txt",
    log:
        "logs/filter_samples/all.log",
    params:
        min_depth=MIN_DEPTH,
        max_unmapped_frac=MAX_UNMAPPED_FRAC,
    run:
        size = genome_size(input.fai)
        rows = []
        for sample, path in zip(SAMPLES, input.stats):
            stats = read_bam_stats(path)
            depth = stats["mapped_bases"] / size
            unmapped_frac = (
                stats["unmapped"] / stats["total"] if stats["total"] else 1.0
            )
            rows.append(
                {
                    "sample": sample,
                    "depth": round(depth, 4),
                    "unmapped_frac": round(unmapped_frac, 4),
                    "pass": depth >= params.min_depth
                    and unmapped_frac <= params.max_unmapped_frac,
                }
            )
        qc = pd.DataFrame(rows)
        passed = qc.loc[qc["pass"], "sample"]
        with open(log[0], "w") as f:
            f.write(
                f"{len(passed)}/{len(qc)} samples passed "
                f"(depth >= {params.min_depth}, unmapped_frac <= {params.max_unmapped_frac})\n"
            )
        if passed.empty:
            raise ValueError("no sample passed the BAM filter")
        qc.to_csv(output.qc, sep="\t", index=False)
        passed.to_csv(output.passed, index=False, header=False)
