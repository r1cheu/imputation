# Snakemake workflow: rice imputation

[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.20-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for genotype imputation of low-coverage rice sequencing data using **GLIMPSE1**.

## Overview

```
fastp ──▶ bwa-mem2 + markdup ──▶ bcftools mpileup + call ──▶ bcftools concat + merge
                                                                       │
panel VCF ──▶ GLIMPSE_chunk ──▶ GLIMPSE_phase ◀────────────────────────┘
                                       │
                                GLIMPSE_ligate ──▶ bcftools concat ──▶ imputed VCF
```

## Requirements

- Linux + [pixi](https://pixi.sh) (manages snakemake itself)
- conda / mamba (auto-managed by snakemake `--sdm conda`)
- A SLURM cluster (optional; local execution also works)

User-supplied resources (paths in `config/config.yaml`):

| Item                                          | Where to put it                            |
| --------------------------------------------- | ------------------------------------------ |
| Reference FASTA                               | `config.reference.fasta`                   |
| Per-chrom full-GT panel VCF.GZ + .csi         | matching `config.panel.full_template`      |
| Per-chrom sites-only TSV.GZ + .tbi            | matching `config.panel.sites_tsv_template` |
| Per-chromosome genetic maps (GLIMPSE1 format) | matching `config.genetic_map.template`     |
| Per-sample paired-end FASTQ                   | listed in `config/samples.tsv`             |

> **GLIMPSE1** binaries (`GLIMPSE_chunk`, `GLIMPSE_phase`, `GLIMPSE_ligate`) must be on `PATH`.
> Follow the [GLIMPSE1 installation guide](https://odelaneau.github.io/GLIMPSE/glimpse1/installation.html) to compile before running the workflow.

## Usage

```bash
pixi install
pixi run envs
pixi run dry
pixi run local
pixi run run
```

See `config/README.md` for the configuration reference.

### Between-workflow caching

Reference- and panel-only rules are marked `cache: True` so their outputs can be
shared across runs/projects:

- `bwa_mem2_index`, `samtools_faidx`

`pixi.toml` defaults `SNAKEMAKE_OUTPUT_CACHE=.snakemake-cache`. To share the cache
across projects (recommended for a fixed reference + panel), set it to a shared
path before invoking pixi tasks:

```bash
export SNAKEMAKE_OUTPUT_CACHE=/shared/snakemake-cache
pixi run run
```

## TODO

- Replace bwa-mem2 with minibwa for alignment

## Layout

```
config/         user-facing config (config.yaml + samples.tsv)
slurm/          snakemake-executor-plugin-slurm profile
workflow/
  Snakefile     main entry; rule all = per-chrom imputed VCF
  rules/        common, reference, trim, align, imputation
  envs/         per-rule conda environments
  schemas/      config and sample-sheet JSON schemas
pixi.toml       pixi project manifest
```

## References

- Rubinacci S, Ribeiro D, Hofmeister R, Delaneau O. _Efficient phasing and imputation of low-coverage sequencing data using large reference panels._ Nature Genetics 53.1 (2021): 120-126. (GLIMPSE1)
- Vasimuddin Md. et al. _Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems._ IPDPS 2019. (bwa-mem2)
- Chen S. _Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp._ iMeta 2023.
- Köster J. et al. _Sustainable data analysis with Snakemake._ F1000Research 2021.
