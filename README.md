# Snakemake workflow: rice imputation

[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.20-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for genotype imputation of low-coverage rice sequencing data using **GLIMPSE2**.

## Overview

```
fastp ──▶ minibwa + markdup ──▶ samtools stats ──▶ sample filter (depth / unmapped)
                                                          │
panel VCF ──▶ GLIMPSE2_split_reference (whole chrom) ──▶ GLIMPSE2_phase (--bam-list, per batch)
                                                          │
                                  bcftools merge (batches) ──▶ bcftools concat ──▶ raw.bcf
                                                          │
                          bcftools +impute-info ──▶ INFO filter ──▶ final.bcf ──▶ plink2 ──▶ final.bed
```

Each chromosome is imputed as a single region (no chunking, so no ligation).
Samples with mapped depth below `bam_filter.min_depth` or a fraction of unmapped
reads above `bam_filter.max_unmapped_frac` are dropped before imputation; the
per-sample metrics are written to `results/qc/sample_qc.tsv`.

Batches are imputed separately, so the INFO score GLIMPSE2 writes only reflects
one batch. It is recomputed across all samples with `bcftools +impute-info`, and
sites below `site_filter.min_info` are dropped. Outputs in `results/imputed/`:

- `raw.bcf` — all imputed sites, unfiltered
- `final.bcf` — sites passing the INFO filter (INFO recomputed)
- `final.{bed,bim,fam}` — `final.bcf` converted with plink2

## Requirements

- Linux + [pixi](https://pixi.sh) (manages snakemake itself)
- conda / mamba (auto-managed by snakemake `--sdm conda`)
- A SLURM cluster (optional; local execution also works)

User-supplied resources (paths in `config/config.yaml`):

| Item                                       | Where to put it                        |
| ------------------------------------------ | -------------------------------------- |
| Reference FASTA                            | `config.reference.fasta`               |
| Per-chrom phased panel VCF.GZ/BCF + .csi   | matching `config.panel.full_template`  |
| Per-chromosome genetic maps (`pos chr cM`) | matching `config.genetic_map.template` |
| Per-sample paired-end FASTQ                | listed in `config/samples.tsv`         |

## Usage

```bash
pixi install
pixi run dry
pixi run local
pixi run run
```

See `config/README.md` for the configuration reference.

### Between-workflow caching

Reference- and panel-only rules are marked `cache: True` so their outputs can be
shared across runs/projects:

- `minibwa_index`, `samtools_faidx`

`pixi.toml` defaults `SNAKEMAKE_OUTPUT_CACHE=.snakemake-cache`. To share the cache
across projects (recommended for a fixed reference + panel), set it to a shared
path before invoking pixi tasks:

```bash
export SNAKEMAKE_OUTPUT_CACHE=/shared/snakemake-cache
pixi run run
```

## Layout

```
config/         user-facing config (config.yaml + samples.tsv)
slurm/          snakemake-executor-plugin-slurm profile
workflow/
  Snakefile     main entry; rule all = INFO-filtered BCF + PLINK bed
  rules/        common, reference, trim, align, qc, imputation
  envs/         per-rule conda environments
  schemas/      config and sample-sheet JSON schemas
pixi.toml       pixi project manifest
```

## References

- Rubinacci S, Hofmeister R, Sousa da Mota B, Delaneau O. _Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes._ Nature Genetics 55 (2023): 1088-1090. (GLIMPSE2)
- Li H. _Minibwa: successor of bwa-mem for short-read alignment._ https://github.com/lh3/minibwa
- Chen S. _Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication using fastp._ iMeta 2023.
- Köster J. et al. _Sustainable data analysis with Snakemake._ F1000Research 2021.
