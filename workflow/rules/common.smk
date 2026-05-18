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

REF_FASTA = config["reference"]["fasta"]
REF_PREFIX = REF_FASTA  # bwa-mem2 index prefix == fasta path
PANEL_FULL = config["panel"]["full_template"]
PANEL_SITES = config["panel"]["sites_template"]
MAP_TEMPLATE = config["genetic_map"]["template"]


def get_fastq(wc):
    row = samples.loc[wc.sample]
    return {"r1": row["fq1"], "r2": row["fq2"]}


def get_trimmed_reads(wc):
    return [
        f"results/trimmed/{wc.sample}.1.fq.gz",
        f"results/trimmed/{wc.sample}.2.fq.gz",
    ]


def get_read_group(wc):
    platform = samples.loc[wc.sample, "platform"]
    return rf"@RG\tID:{wc.sample}\tSM:{wc.sample}\tLB:{wc.sample}\tPL:{platform}"


def get_map(wc):
    return MAP_TEMPLATE.format(chrom=wc.chrom)


CHUNK_COLS = ["idx", "chrom", "input_region", "output_region"]


def read_chunks(chrom):
    path = checkpoints.glimpse_chunk.get(chrom=chrom).output[0]
    df = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        comment="#",
        usecols=[0, 1, 2, 3],
        names=CHUNK_COLS,
        index_col=False,
        dtype=str,
    )
    return df


def get_chunk_region(wc, kind):
    df = read_chunks(wc.chrom)
    row = df.loc[df["idx"] == wc.idx]
    if row.empty:
        raise ValueError(f"chunk {wc.idx} not in chunks file for {wc.chrom}")
    return row.iloc[0][kind]


def phased_chunks(wc):
    return [
        f"results/phased/{wc.chrom}/chunk_{idx}.bcf"
        for idx in read_chunks(wc.chrom)["idx"]
    ]


def phased_chunks_idx(wc):
    return [
        f"results/phased/{wc.chrom}/chunk_{idx}.bcf.csi"
        for idx in read_chunks(wc.chrom)["idx"]
    ]
