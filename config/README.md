# Configuration

## Bash runner: `pipeline.env`

| Key | Meaning |
|---|---|
| `SAMPLE_SHEET` | TSV listing samples (see below) |
| `REFERENCE_FASTA` | Reference genome FASTA. `bwa-mem2 index` and `samtools faidx` are run by the pipeline. |
| `PANEL_FULL_TEMPLATE` | Per-chrom full-GT panel BCF template. Use `{chrom}` for the chromosome placeholder. |
| `PANEL_SITES_TSV_TEMPLATE` | Per-chrom sites TSV template, bgzipped with `.tbi` alongside. Use `{chrom}` for the chromosome placeholder. |
| `MAP_TEMPLATE` | Per-chromosome genetic map template. Use `{chrom}` for the chromosome placeholder. |
| `CHROMS` | Bash array of chromosome names. Must match reference, panel, sites, and map files. |
| `FASTP_EXTRA` | Extra flags passed verbatim to `fastp`. |
| `GLIMPSE_CHUNK_WINDOW_SIZE` | `GLIMPSE_chunk --window-size`. |
| `GLIMPSE_CHUNK_BUFFER_SIZE` | `GLIMPSE_chunk --buffer-size`. |
| `GLIMPSE_CHUNK_EXTRA` | Extra flags passed verbatim to `GLIMPSE_chunk`. |
| `SLURM_ACCOUNT` | Account passed to `sbatch --account`. |
| `SLURM_PARTITION` | Partition passed to `sbatch --partition`. |
| `*_CPUS`, `*_MEM`, `*_TIME` | Per-stage SLURM resources. |
| `*_ARRAY_LIMIT` | Per-stage maximum simultaneous SLURM array tasks. |

Run with:

```bash
bash scripts/submit_pipeline.sh config/pipeline.env
```

## `samples.tsv`

Tab-separated, one sample per row.

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | unique sample id (used as RG ID/SM) |
| `platform` | yes | sequencing platform string (RG PL), e.g. `ILLUMINA` |
| `fq1` | yes | path to read 1 fastq.gz |
| `fq2` | yes | path to read 2 fastq.gz |

Rows are sorted by `sample` before submission, matching the legacy Snakemake
workflow behavior.

## Legacy Snakemake config

`config/config.yaml` is still used by the old Snakemake workflow under
`workflow/`. It is not read by the bash/SLURM runner.
