#!/usr/bin/env bash
# Download and prepare test data for end-to-end pipeline test.
# Source: GLIMPSE2 tutorial https://odelaneau.github.io/GLIMPSE/docs/tutorials/getting_started/
# Run from repo root: bash .test/setup_data.sh
set -euo pipefail

cd "$(dirname "$0")"
mkdir -p resources/{reference,panel,maps} data

PIXI_RUN="pixi run -e test"

echo "[1/6] Reference chr22 FASTA"
$PIXI_RUN wget -q -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz -O resources/reference/chr22.fa.gz
gunzip -kf resources/reference/chr22.fa.gz

echo "[2/6] Genetic map"
$PIXI_RUN wget -q -c https://raw.githubusercontent.com/odelaneau/GLIMPSE/master/maps/genetic_maps.b38/chr22.b38.gmap.gz -O resources/maps/chr22.b38.gmap.gz
gunzip -c resources/maps/chr22.b38.gmap.gz | awk 'NR==1{print; next} {$2="chr22"; print}' OFS='\t' > resources/maps/chr22.gmap

echo "[3/6] NA12878 1x BAM"
$PIXI_RUN wget -q -c https://raw.githubusercontent.com/odelaneau/GLIMPSE/master/tutorial/NA12878_1x_bam/NA12878.bam -O data/NA12878.bam
$PIXI_RUN wget -q -c https://raw.githubusercontent.com/odelaneau/GLIMPSE/master/tutorial/NA12878_1x_bam/NA12878.bam.bai -O data/NA12878.bam.bai

echo "[4/6] BAM -> paired FASTQ"
$PIXI_RUN bash -c '
  samtools sort -n -@ 4 -o /tmp/NA12878.namesorted.bam data/NA12878.bam
  samtools fastq -@ 4 -1 data/NA12878_1.fq.gz -2 data/NA12878_2.fq.gz -0 /dev/null -s /dev/null -n /tmp/NA12878.namesorted.bam
  rm -f /tmp/NA12878.namesorted.bam
'

echo "[5/6] 1000G chr22 phased panel (~500 MB)"
PANEL_BASE=http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased
PANEL_FILE=CCDG_14151_B01_GRM_WGS_2020-08-05_chr22.filtered.shapeit2-duohmm-phased.vcf.gz
$PIXI_RUN wget -q -c "$PANEL_BASE/$PANEL_FILE" -O resources/panel/raw.vcf.gz
$PIXI_RUN wget -q -c "$PANEL_BASE/$PANEL_FILE.tbi" -O resources/panel/raw.vcf.gz.tbi

echo "[6/6] Filter panel + emit full BCF (with .csi) and sites VCF (with .csi)"
$PIXI_RUN bash -c '
  bcftools norm -m -any resources/panel/raw.vcf.gz -Ou --threads 4 |
  bcftools view -m 2 -M 2 -v snps -s ^NA12878,NA12891,NA12892 --threads 4 -Ob -o resources/panel/1000GP.chr22.noNA12878.bcf
  bcftools index -f resources/panel/1000GP.chr22.noNA12878.bcf
  bcftools view -G resources/panel/1000GP.chr22.noNA12878.bcf -Oz -o resources/panel/1000GP.chr22.noNA12878.sites.vcf.gz
  bcftools index -f resources/panel/1000GP.chr22.noNA12878.sites.vcf.gz
  rm -f resources/panel/raw.vcf.gz resources/panel/raw.vcf.gz.tbi
  rm -f resources/reference/chr22.fa.gz resources/maps/chr22.b38.gmap.gz
'

echo "Done. Run: pixi run snakemake --directory .test --sdm conda --cores 8"
