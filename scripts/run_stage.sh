#!/usr/bin/env bash
set -euo pipefail

if (( $# != 2 )); then
    printf 'usage: %s <stage> <config.env>\n' "$0" >&2
    exit 2
fi

STAGE=$1
CONFIG_PATH=$2
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

if [[ "$CONFIG_PATH" != /* ]]; then
    CONFIG_PATH="$ROOT/$CONFIG_PATH"
fi

cd "$ROOT"
source "$CONFIG_PATH"

export PATH="$ROOT/bin:$PATH"

awk -F '\t' '
    NR == 1 {
        for (i = 1; i <= NF; i++) column[$i] = i
        exit !(("sample" in column) && ("platform" in column) && ("fq1" in column) && ("fq2" in column))
    }
' "$SAMPLE_SHEET"

mapfile -t SAMPLE_ROWS < <(
    awk -F '\t' '
        NR == 1 {
            for (i = 1; i <= NF; i++) column[$i] = i
            next
        }
        NF > 1 {
            print $column["sample"] "\t" $column["platform"] "\t" $column["fq1"] "\t" $column["fq2"]
        }
    ' "$SAMPLE_SHEET" | sort -k1,1
)

SAMPLE_COUNT=${#SAMPLE_ROWS[@]}
CHROM_COUNT=${#CHROMS[@]}

case "$STAGE" in
    reference)
        mkdir -p logs/reference

        BWA_OUTPUTS=(
            "$REFERENCE_FASTA.0123"
            "$REFERENCE_FASTA.amb"
            "$REFERENCE_FASTA.ann"
            "$REFERENCE_FASTA.bwt.2bit.64"
            "$REFERENCE_FASTA.pac"
        )
        BWA_INDEX_DONE=1
        for BWA_OUTPUT in "${BWA_OUTPUTS[@]}"; do
            if [[ ! -s "$BWA_OUTPUT" ]]; then
                BWA_INDEX_DONE=0
            fi
        done

        if (( BWA_INDEX_DONE == 0 )); then
            rm -f "${BWA_OUTPUTS[@]}"
            bwa-mem2 index "$REFERENCE_FASTA" > logs/reference/bwa_mem2_index.log 2>&1
        fi

        if [[ ! -s "$REFERENCE_FASTA.fai" ]]; then
            rm -f "$REFERENCE_FASTA.fai"
            samtools faidx "$REFERENCE_FASTA" > logs/reference/samtools_faidx.log 2>&1
        fi
        ;;

    trim)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for trim}
        SAMPLE_INDEX=$((TASK_ID - 1))
        IFS=$'\t' read -r SAMPLE PLATFORM FQ1 FQ2 <<< "${SAMPLE_ROWS[$SAMPLE_INDEX]}"

        mkdir -p results/trimmed logs/fastp
        TRIM_R1="results/trimmed/$SAMPLE.1.fq.gz"
        TRIM_R2="results/trimmed/$SAMPLE.2.fq.gz"
        LOG="logs/fastp/$SAMPLE.log"

        if [[ -s "$TRIM_R1" && -s "$TRIM_R2" ]]; then
            exit 0
        fi

        rm -f "$TRIM_R1" "$TRIM_R2"
        fastp -i "$FQ1" -I "$FQ2" -o "$TRIM_R1" -O "$TRIM_R2" -w "$TRIM_CPUS" $FASTP_EXTRA > "$LOG" 2>&1
        ;;

    align)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for align}
        SAMPLE_INDEX=$((TASK_ID - 1))
        IFS=$'\t' read -r SAMPLE PLATFORM FQ1 FQ2 <<< "${SAMPLE_ROWS[$SAMPLE_INDEX]}"

        mkdir -p results/dedup logs/bwa_dedup
        TRIM_R1="results/trimmed/$SAMPLE.1.fq.gz"
        TRIM_R2="results/trimmed/$SAMPLE.2.fq.gz"
        BAM="results/dedup/$SAMPLE.bam"
        BAI="results/dedup/$SAMPLE.bam.bai"
        LOG="logs/bwa_dedup/$SAMPLE.log"
        READ_GROUP="@RG\tID:$SAMPLE\tSM:$SAMPLE\tLB:$SAMPLE\tPL:$PLATFORM"

        if [[ -s "$BAM" && -s "$BAI" ]]; then
            exit 0
        fi

        rm -f "$BAM" "$BAI"
        {
            bwa-mem2 mem -t "$ALIGN_CPUS" -R "$READ_GROUP" "$REFERENCE_FASTA" "$TRIM_R1" "$TRIM_R2" |
                samtools fixmate -@ "$ALIGN_CPUS" -m -u - - |
                samtools sort -@ "$ALIGN_CPUS" -u -m 1G - |
                samtools markdup -@ "$ALIGN_CPUS" -r --write-index - "$BAM"
        } > "$LOG" 2>&1
        ;;

    compute_gl)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for compute_gl}
        TASK_INDEX=$((TASK_ID - 1))
        SAMPLE_INDEX=$((TASK_INDEX / CHROM_COUNT))
        CHROM_INDEX=$((TASK_INDEX % CHROM_COUNT))
        IFS=$'\t' read -r SAMPLE PLATFORM FQ1 FQ2 <<< "${SAMPLE_ROWS[$SAMPLE_INDEX]}"
        CHROM=${CHROMS[$CHROM_INDEX]}

        mkdir -p "results/gl/$SAMPLE" logs/compute_gl
        BAM="results/dedup/$SAMPLE.bam"
        PANEL_FULL=${PANEL_FULL_TEMPLATE//\{chrom\}/$CHROM}
        PANEL_SITES_TSV=${PANEL_SITES_TSV_TEMPLATE//\{chrom\}/$CHROM}
        OUT_BCF="results/gl/$SAMPLE/$CHROM.bcf"
        OUT_CSI="results/gl/$SAMPLE/$CHROM.bcf.csi"
        LOG="logs/compute_gl/${SAMPLE}_${CHROM}.log"

        if [[ -s "$OUT_BCF" && -s "$OUT_CSI" ]]; then
            exit 0
        fi

        rm -f "$OUT_BCF" "$OUT_CSI"
        {
            bcftools mpileup -f "$REFERENCE_FASTA" -I -E -a 'FORMAT/DP' -T "$PANEL_FULL" -r "$CHROM" "$BAM" -Ou |
                bcftools call -Aim -C alleles -T "$PANEL_SITES_TSV" -Ob -o "$OUT_BCF"
            bcftools index -f "$OUT_BCF"
        } > "$LOG" 2>&1
        ;;

    concat_gl)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for concat_gl}
        SAMPLE_INDEX=$((TASK_ID - 1))
        IFS=$'\t' read -r SAMPLE PLATFORM FQ1 FQ2 <<< "${SAMPLE_ROWS[$SAMPLE_INDEX]}"

        mkdir -p results/gl logs/concat_gl
        OUT_BCF="results/gl/$SAMPLE.bcf"
        OUT_CSI="results/gl/$SAMPLE.bcf.csi"
        LOG="logs/concat_gl/$SAMPLE.log"

        if [[ -s "$OUT_BCF" && -s "$OUT_CSI" ]]; then
            exit 0
        fi

        GL_INPUTS=()
        for CHROM in "${CHROMS[@]}"; do
            GL_INPUTS+=("results/gl/$SAMPLE/$CHROM.bcf")
        done

        rm -f "$OUT_BCF" "$OUT_CSI"
        {
            bcftools concat --threads "$CONCAT_GL_CPUS" -Ob -o "$OUT_BCF" "${GL_INPUTS[@]}"
            bcftools index -f --threads "$CONCAT_GL_CPUS" "$OUT_BCF"
        } > "$LOG" 2>&1
        ;;

    merge_gl)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for merge_gl}
        CHROM_INDEX=$((TASK_ID - 1))
        CHROM=${CHROMS[$CHROM_INDEX]}

        mkdir -p results/gl_merged logs/merge_gl
        OUT_BCF="results/gl_merged/$CHROM.bcf"
        OUT_CSI="results/gl_merged/$CHROM.bcf.csi"
        LOG="logs/merge_gl/$CHROM.log"

        if [[ -s "$OUT_BCF" && -s "$OUT_CSI" ]]; then
            exit 0
        fi

        SAMPLE_GL_INPUTS=()
        for SAMPLE_ROW in "${SAMPLE_ROWS[@]}"; do
            IFS=$'\t' read -r SAMPLE PLATFORM FQ1 FQ2 <<< "$SAMPLE_ROW"
            SAMPLE_GL_INPUTS+=("results/gl/$SAMPLE.bcf")
        done

        rm -f "$OUT_BCF" "$OUT_CSI"
        {
            bcftools merge -m none -r "$CHROM" --threads "$MERGE_GL_CPUS" -Ob -o "$OUT_BCF" "${SAMPLE_GL_INPUTS[@]}"
            bcftools index -f "$OUT_BCF"
        } > "$LOG" 2>&1
        ;;

    chunk)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for chunk}
        CHROM_INDEX=$((TASK_ID - 1))
        CHROM=${CHROMS[$CHROM_INDEX]}

        mkdir -p results/chunks logs/glimpse_chunk
        PANEL_FULL=${PANEL_FULL_TEMPLATE//\{chrom\}/$CHROM}
        MAP=${MAP_TEMPLATE//\{chrom\}/$CHROM}
        OUT_TXT="results/chunks/$CHROM.txt"
        LOG="logs/glimpse_chunk/$CHROM.log"

        if [[ -s "$OUT_TXT" ]]; then
            exit 0
        fi

        rm -f "$OUT_TXT"
        GLIMPSE_chunk --input "$PANEL_FULL" --map "$MAP" --region "$CHROM" \
            --window-size "$GLIMPSE_CHUNK_WINDOW_SIZE" --buffer-size "$GLIMPSE_CHUNK_BUFFER_SIZE" \
            --thread "$CHUNK_CPUS" --output "$OUT_TXT" $GLIMPSE_CHUNK_EXTRA > "$LOG" 2>&1
        ;;

    phase)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for phase}
        PHASE_TASK=$(sed -n "${TASK_ID}p" results/chunks/phase_tasks.tsv)
        IFS=$'\t' read -r CHROM CHUNK_INDEX INPUT_REGION OUTPUT_REGION <<< "$PHASE_TASK"

        mkdir -p "results/phased/$CHROM" logs/glimpse_phase
        PANEL_FULL=${PANEL_FULL_TEMPLATE//\{chrom\}/$CHROM}
        MAP=${MAP_TEMPLATE//\{chrom\}/$CHROM}
        GL_BCF="results/gl_merged/$CHROM.bcf"
        OUT_BCF="results/phased/$CHROM/chunk_$CHUNK_INDEX.bcf"
        OUT_CSI="results/phased/$CHROM/chunk_$CHUNK_INDEX.bcf.csi"
        LOG="logs/glimpse_phase/${CHROM}_chunk_${CHUNK_INDEX}.log"

        if [[ -s "$OUT_BCF" && -s "$OUT_CSI" ]]; then
            exit 0
        fi

        rm -f "$OUT_BCF" "$OUT_CSI"
        {
            GLIMPSE_phase --input "$GL_BCF" --reference "$PANEL_FULL" --map "$MAP" \
                --input-region "$INPUT_REGION" --output-region "$OUTPUT_REGION" \
                --thread "$PHASE_CPUS" --output "$OUT_BCF"
            bcftools index -f "$OUT_BCF"
        } > "$LOG" 2>&1
        ;;

    ligate)
        TASK_ID=${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required for ligate}
        CHROM_INDEX=$((TASK_ID - 1))
        CHROM=${CHROMS[$CHROM_INDEX]}

        mkdir -p "results/phased/$CHROM" results/imputed logs/glimpse_ligate
        LIST_FILE="results/phased/$CHROM/ligate.list"
        OUT_BCF="results/imputed/$CHROM.bcf"
        OUT_CSI="results/imputed/$CHROM.bcf.csi"
        LOG="logs/glimpse_ligate/$CHROM.log"

        if [[ -s "$OUT_BCF" && -s "$OUT_CSI" ]]; then
            exit 0
        fi

        awk -v chrom="$CHROM" '
            !/^#/ && NF >= 4 {
                print "results/phased/" chrom "/chunk_" $1 ".bcf"
            }
        ' "results/chunks/$CHROM.txt" > "$LIST_FILE"

        rm -f "$OUT_BCF" "$OUT_CSI"
        {
            GLIMPSE_ligate --input "$LIST_FILE" --output "$OUT_BCF" --thread "$LIGATE_CPUS"
            bcftools index -f "$OUT_BCF" --threads "$LIGATE_CPUS"
        } > "$LOG" 2>&1
        ;;

    concat_imputed)
        mkdir -p results/imputed logs/concat_imputed
        OUT_BCF="results/imputed/all.bcf"
        OUT_CSI="results/imputed/all.bcf.csi"
        LOG="logs/concat_imputed/all.log"

        if [[ -s "$OUT_BCF" && -s "$OUT_CSI" ]]; then
            exit 0
        fi

        IMPUTED_INPUTS=()
        for CHROM in "${CHROMS[@]}"; do
            IMPUTED_INPUTS+=("results/imputed/$CHROM.bcf")
        done

        rm -f "$OUT_BCF" "$OUT_CSI"
        {
            bcftools concat --threads "$CONCAT_IMPUTED_CPUS" -Ob -o "$OUT_BCF" "${IMPUTED_INPUTS[@]}"
            bcftools index -f --threads "$CONCAT_IMPUTED_CPUS" "$OUT_BCF"
        } > "$LOG" 2>&1
        ;;

    *)
        printf 'unknown stage: %s\n' "$STAGE" >&2
        exit 2
        ;;
esac
