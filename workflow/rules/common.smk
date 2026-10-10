from pathlib import Path

import pandas as pd
from snakemake.utils import validate

samples = (
    pd.read_csv(config["sample_sheet"], sep="\t", dtype=str)
    .set_index("sample", drop=False)
    .sort_index()
)
validate(samples, schema="../schemas/samples.schema.yaml")
validate(config, schema="../schemas/config.schema.yaml")

SAMPLES = samples.index.tolist()
CHROMS = config["chromosomes"]
PHASE_BATCH_SIZE = int(config.get("phase", {}).get("batch_size", 200))
PANEL = config["panel"]["full_template"]
MIN_DEPTH = float(config.get("bam_filter", {}).get("min_depth", 1.0))
MAX_UNMAPPED_FRAC = float(config.get("bam_filter", {}).get("max_unmapped_frac", 0.8))


wildcard_constraints:
    sample="|".join(SAMPLES),
    chrom="|".join(CHROMS),
    batch=r"\d+",


def get_fastq(wc):
    row = samples.loc[wc.sample]
    return {"r1": row["fq1"], "r2": row["fq2"]}


def get_read_group(wc):
    platform = samples.loc[wc.sample, "platform"]
    return rf"@RG\tID:{wc.sample}\tSM:{wc.sample}\tLB:{wc.sample}\tPL:{platform}"


def get_map(wc):
    return config["genetic_map"]["template"].format(chrom=wc.chrom)


STATS_KEYS = {
    "raw total sequences:": "total",
    "reads unmapped:": "unmapped",
    "bases mapped (cigar):": "mapped_bases",
}


def read_bam_stats(path):
    """Parse the SN section of a `samtools stats` report."""
    stats = {}
    with open(path) as f:
        for line in f:
            if not line.startswith("SN\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[1] in STATS_KEYS:
                stats[STATS_KEYS[fields[1]]] = int(fields[2])
    return stats


def genome_size(fai):
    return int(pd.read_csv(fai, sep="\t", header=None, usecols=[1])[1].sum())


def chrom_region(wc, input):
    """Whole-chromosome region (chrX:1-LEN) as required by GLIMPSE2_split_reference."""
    lengths = pd.read_csv(
        input.fai, sep="\t", header=None, usecols=[0, 1], index_col=0
    )[1]
    return f"{wc.chrom}:1-{lengths[wc.chrom]}"


def passed_samples():
    path = checkpoints.filter_samples.get().output.passed
    with open(path) as f:
        return [line.strip() for line in f if line.strip()]


def phase_batches():
    n = len(passed_samples())
    return [f"{i:04d}" for i in range((n + PHASE_BATCH_SIZE - 1) // PHASE_BATCH_SIZE)]


def batch_samples(wc):
    start = int(wc.batch) * PHASE_BATCH_SIZE
    return passed_samples()[start : start + PHASE_BATCH_SIZE]


def batch_bams(wc):
    return [f"results/dedup/{sample}.bam" for sample in batch_samples(wc)]


def batch_bais(wc):
    return [f"results/dedup/{sample}.bam.bai" for sample in batch_samples(wc)]


def phased_batches(wc):
    return [f"results/phased/{wc.chrom}/batch_{batch}.bcf" for batch in phase_batches()]


def phased_batches_idx(wc):
    return [
        f"results/phased/{wc.chrom}/batch_{batch}.bcf.csi" for batch in phase_batches()
    ]
