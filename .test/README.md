# Test data

End-to-end smoke test using GLIMPSE2 tutorial data (NA12878 chr22 1× WGS).
**Not run in CI** — downloads ~600 MB and takes ~15 min wall time.

## Setup

```bash
pixi install -e test
bash .test/setup_data.sh
```

`setup_data.sh` downloads:
- chr22 reference (UCSC GRCh38)
- 1000G chr22 phased panel (filtered to drop NA12878 family)
- GLIMPSE2 chr22 genetic map
- NA12878 1× BAM from GLIMPSE GitHub, converted to PE FASTQ via `samtools fastq`

## Run

```bash
pixi run snakemake --directory .test --sdm conda --cores 8
```

Expected output: `.test/results/imputed/chr22.vcf.gz` (~10 MB, ~930k SNPs).

## Files in git

Only `config/` is tracked. Data, results, logs, snakemake state are gitignored.
