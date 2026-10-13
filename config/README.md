# Configuration

## `config.yaml`

| Key | Meaning |
|---|---|
| `sample_sheet` | TSV listing samples (see below) |
| `reference.fasta` | Reference genome FASTA (e.g. IRGSP-1.0). Indices are built by the workflow. |
| `panel.full_template` | Per-chromosome phased reference panel template (VCF.GZ or BCF, with `.csi` alongside), e.g. `resources/ref_panel/{chrom}.vcf.gz`. |
| `chromosomes` | List of chromosome names. Must match both reference and panel. |
| `genetic_map.template` | Per-chromosome genetic map path template, e.g. `resources/map/{chrom}.map`. Format expected by GLIMPSE2: `pos chr cM`. |
| `bam_filter.min_depth` | Drop samples whose mapped depth (mapped bases after duplicate removal / reference length) is below this (default 1.0). |
| `bam_filter.max_unmapped_frac` | Drop samples whose fraction of unmapped reads is above this (default 0.8). |
| `phase.batch_size` | Maximum number of BAMs imputed in one GLIMPSE2_phase job (default 200). |
| `site_filter.min_info` | Keep sites whose IMPUTE INFO score, recomputed across all samples with `bcftools +impute-info`, is at least this (default 0.95). |
| `fastp.extra` | Extra flags passed verbatim to fastp (e.g. quality cutoffs) |

Threads/memory per rule are hard-coded in `workflow/rules/*.smk` (tuned for ~600 cores / ~5000 samples). Override by editing those `resources:` blocks.

## `samples.tsv`

Tab-separated, one sample per row.

| Column | Required | Meaning |
|---|---|---|
| `sample` | yes | unique sample id (used as RG ID/SM) |
| `platform` | yes | sequencing platform string (RG PL), e.g. `ILLUMINA` |
| `fq1` | yes | path to read 1 fastq.gz |
| `fq2` | yes | path to read 2 fastq.gz |
