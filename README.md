# Snakemake workflow: rice imputation

[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A58.20-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow for genotype imputation of low-coverage rice sequencing data using **GLIMPSE2**.

## Overview

```
FASTQ ──▶ fastqc / fastp trim ──▶ bwa-mem2 mem ──▶ samtools sort + markdup ──▶ BAM
                                                                                │
phased panel VCF ──▶ split per-chrom ──▶ GLIMPSE2_chunk ──▶ split_reference     │
                                                                │               │
                                                                └──▶ GLIMPSE2_phase (--bam-list)
                                                                              │
                                                                              ▼
                                                              GLIMPSE2_ligate ──▶ per-chrom VCF
```

## Requirements

- Linux + [pixi](https://pixi.sh) (manages snakemake itself)
- conda / mamba (auto-managed by snakemake `--sdm conda`)
- A SLURM cluster (optional; local execution also works)

User-supplied resources (paths in `config/config.yaml`):

| Item | Where to put it |
|---|---|
| Reference FASTA | `reference.fasta` |
| Phased reference panel VCF (whole genome, bgzipped) | `panel.vcf` |
| Per-chromosome genetic maps in GLIMPSE2 format | matching `genetic_map.template` |
| Per-sample paired-end FASTQ | listed in `config/samples.tsv` |

## Usage

```bash
pixi install                       # install snakemake + slurm plugin into .pixi/
pixi run envs                      # pre-create per-rule conda envs (optional)
pixi run dry                       # dry-run
pixi run local                     # local execution with 8 cores
pixi run run                       # SLURM submission via slurm/config.yaml
```

See `config/README.md` for the configuration reference.

## Layout

```
config/         user-facing config (config.yaml + samples.tsv)
slurm/          snakemake-executor-plugin-slurm profile
workflow/
  Snakefile     main entry; rule all = per-chrom imputed VCF + multiqc report
  rules/        common, reference, qc_trim, align, imputation, multiqc
  envs/         per-rule conda environments
  schemas/      config and sample-sheet JSON schemas
pixi.toml       pixi project manifest
```

## References

- Rubinacci S. et al. *Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes.* Nat. Genet. 2023. (GLIMPSE2)
- Vasimuddin Md. et al. *Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.* IPDPS 2019. (bwa-mem2)
- Köster J. et al. *Sustainable data analysis with Snakemake.* F1000Research 2021.
