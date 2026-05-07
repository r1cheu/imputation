# Configuration

## `config.yaml`

| Key | Meaning |
|---|---|
| `sample_sheet` | TSV listing samples (see below) |
| `reference.fasta` | Reference genome FASTA (e.g. IRGSP-1.0). Indices are built by the workflow. |
| `panel.vcf` | Whole-genome phased reference panel VCF, bgzipped. Will be split per chromosome. |
| `chromosomes` | List of chromosome names. Must match both reference and panel. |
| `genetic_map.template` | Per-chromosome genetic map path template, e.g. `resources/maps/{chrom}.gmap`. Format expected by GLIMPSE2: `pos chr cM`. |
| `glimpse2_chunk.window_mb` | GLIMPSE2_chunk `--window-mb` (default 4.0) |
| `glimpse2_chunk.buffer_mb` | GLIMPSE2_chunk `--buffer-mb` (default 0.5) |
| `glimpse2_chunk.extra` | Extra flags passed verbatim to GLIMPSE2_chunk |
| `fastp.threads` | Threads per fastp job |
| `bwa_mem2.threads` | Threads per bwa-mem2 mem job |
| `glimpse2_phase.threads` | Threads per GLIMPSE2_phase job |

## `samples.tsv`

Tab-separated, one sample per row.

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | unique sample id (used as RG ID/SM) |
| `platform` | yes | sequencing platform string (RG PL), e.g. `ILLUMINA` |
| `fq1` | yes | path to read 1 fastq.gz |
| `fq2` | yes | path to read 2 fastq.gz |
