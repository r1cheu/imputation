# Rice imputation pipeline

A bash/SLURM pipeline for genotype imputation of low-coverage rice sequencing data.

The old Snakemake workflow is still kept under `workflow/` for comparison and
fallback, but the recommended HPC entry point is now the bash runner in
`scripts/`.

## Overview

```
FASTQ -> fastp trim -> bwa-mem2 mem -> samtools fixmate/sort/markdup -> BAM
                                                                           |
panel BCF + maps -> GLIMPSE_chunk                                          |
                                                                           v
BAM -> genotype likelihoods -> merged per-chrom GL -> GLIMPSE_phase -> GLIMPSE_ligate -> all.bcf
```

## Requirements

- Linux with bash and a SLURM cluster
- `sbatch`
- `fastp`, `bwa-mem2`, `samtools`, `bcftools`
- `GLIMPSE_chunk`, `GLIMPSE_phase`, `GLIMPSE_ligate`

`scripts/install_bins.sh` can install the native binaries into `bin/`:

```bash
bash scripts/install_bins.sh
export PATH="$PWD/bin:$PATH"
```

GLIMPSE is expected to already be available on the cluster `PATH`.

User-supplied resources (paths in `config/pipeline.env`):

| Item | Where to put it |
|---|---|
| Reference FASTA | `REFERENCE_FASTA` in `config/pipeline.env` |
| Per-chrom full-GT panel BCF + .csi | matching `PANEL_FULL_TEMPLATE` |
| Per-chrom sites TSV + .tbi | matching `PANEL_SITES_TSV_TEMPLATE` |
| Per-chromosome genetic maps | matching `MAP_TEMPLATE` |
| Per-sample paired-end FASTQ | listed in `config/samples.tsv` |

## Usage

```bash
bash scripts/submit_pipeline.sh config/pipeline.env
```

The submit script creates the SLURM dependency chain and prints the submitted
job ids. Dynamic GLIMPSE phase jobs are submitted by a continuation job after
all chunk files are available.

The bash runner skips stages whose final outputs already exist, so rerunning the
same command resumes from the missing or failed outputs.

## Configuration

- Edit `config/pipeline.env` for paths, chromosomes, GLIMPSE/fastp options, and
  SLURM resource settings.
- Edit `config/samples.tsv` for sample ids, platform strings, and FASTQ paths.
- See `config/README.md` for the field reference.

## Snakemake fallback

The original Snakemake workflow remains available:

```bash
pixi run run
```

Use this only when you want to compare behavior with the previous workflow.

## Layout

```
config/         user-facing config (config.yaml + samples.tsv)
scripts/        bash/SLURM runner and native binary installer
slurm/          legacy snakemake-executor-plugin-slurm profile
workflow/
  Snakefile     legacy Snakemake entry
  rules/        legacy Snakemake rules
  schemas/      config and sample-sheet JSON schemas
pixi.toml       legacy Snakemake/pixi project manifest
```

## References

- Rubinacci S. et al. *Imputation of low-coverage sequencing data from 150,119 UK Biobank genomes.* Nat. Genet. 2023. (GLIMPSE2)
- Vasimuddin Md. et al. *Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems.* IPDPS 2019. (bwa-mem2)
- Köster J. et al. *Sustainable data analysis with Snakemake.* F1000Research 2021.
