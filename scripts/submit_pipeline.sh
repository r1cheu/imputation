#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
CONTINUE_AFTER_CHUNK=0

if [[ "${1:-}" == "--continue-after-chunk" ]]; then
    CONTINUE_AFTER_CHUNK=1
    CONFIG_PATH=${2:-config/pipeline.env}
else
    CONFIG_PATH=${1:-config/pipeline.env}
fi

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

mkdir -p "$SLURM_OUTPUT_DIR"

if (( CONTINUE_AFTER_CHUNK == 1 )); then
    mkdir -p results/chunks
    PHASE_TASKS_TMP="results/chunks/phase_tasks.tsv.tmp"
    PHASE_TASKS="results/chunks/phase_tasks.tsv"
    : > "$PHASE_TASKS_TMP"

    for CHROM in "${CHROMS[@]}"; do
        awk -v chrom="$CHROM" '
            !/^#/ && NF >= 4 {
                print chrom "\t" $1 "\t" $3 "\t" $4
            }
        ' "results/chunks/$CHROM.txt" >> "$PHASE_TASKS_TMP"
    done

    mv "$PHASE_TASKS_TMP" "$PHASE_TASKS"
    PHASE_COUNT=$(wc -l < "$PHASE_TASKS")

    if (( PHASE_COUNT == 0 )); then
        printf 'no GLIMPSE phase tasks found in %s\n' "$PHASE_TASKS" >&2
        exit 1
    fi

    PHASE_JOB=$(sbatch --parsable \
        --job-name=imp_phase \
        --account="$SLURM_ACCOUNT" \
        --partition="$SLURM_PARTITION" \
        --cpus-per-task="$PHASE_CPUS" \
        --mem="$PHASE_MEM" \
        --time="$PHASE_TIME" \
        --array="1-${PHASE_COUNT}%${PHASE_ARRAY_LIMIT}" \
        --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
        --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
        --chdir="$ROOT" \
        "$ROOT/scripts/run_stage.sh" phase "$CONFIG_PATH")

    LIGATE_JOB=$(sbatch --parsable \
        --job-name=imp_ligate \
        --account="$SLURM_ACCOUNT" \
        --partition="$SLURM_PARTITION" \
        --cpus-per-task="$LIGATE_CPUS" \
        --mem="$LIGATE_MEM" \
        --time="$LIGATE_TIME" \
        --array="1-${CHROM_COUNT}%${LIGATE_ARRAY_LIMIT}" \
        --dependency="afterok:$PHASE_JOB" \
        --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
        --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
        --chdir="$ROOT" \
        "$ROOT/scripts/run_stage.sh" ligate "$CONFIG_PATH")

    CONCAT_IMPUTED_JOB=$(sbatch --parsable \
        --job-name=imp_concat_all \
        --account="$SLURM_ACCOUNT" \
        --partition="$SLURM_PARTITION" \
        --cpus-per-task="$CONCAT_IMPUTED_CPUS" \
        --mem="$CONCAT_IMPUTED_MEM" \
        --time="$CONCAT_IMPUTED_TIME" \
        --dependency="afterok:$LIGATE_JOB" \
        --output="$SLURM_OUTPUT_DIR/%x_%j.out" \
        --error="$SLURM_OUTPUT_DIR/%x_%j.err" \
        --chdir="$ROOT" \
        "$ROOT/scripts/run_stage.sh" concat_imputed "$CONFIG_PATH")

    printf 'phase: %s\n' "$PHASE_JOB"
    printf 'ligate: %s\n' "$LIGATE_JOB"
    printf 'concat_imputed: %s\n' "$CONCAT_IMPUTED_JOB"
    exit 0
fi

for TOOL in sbatch fastp bwa-mem2 samtools bcftools GLIMPSE_chunk GLIMPSE_phase GLIMPSE_ligate; do
    command -v "$TOOL" >/dev/null || {
        printf 'missing required tool on PATH: %s\n' "$TOOL" >&2
        exit 1
    }
done

if (( SAMPLE_COUNT == 0 )); then
    printf 'no samples found in %s\n' "$SAMPLE_SHEET" >&2
    exit 1
fi

if (( CHROM_COUNT == 0 )); then
    printf 'no chromosomes configured in %s\n' "$CONFIG_PATH" >&2
    exit 1
fi

COMPUTE_GL_COUNT=$((SAMPLE_COUNT * CHROM_COUNT))

REFERENCE_JOB=$(sbatch --parsable \
    --job-name=imp_reference \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$REFERENCE_CPUS" \
    --mem="$REFERENCE_MEM" \
    --time="$REFERENCE_TIME" \
    --output="$SLURM_OUTPUT_DIR/%x_%j.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%j.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" reference "$CONFIG_PATH")

TRIM_JOB=$(sbatch --parsable \
    --job-name=imp_trim \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$TRIM_CPUS" \
    --mem="$TRIM_MEM" \
    --time="$TRIM_TIME" \
    --array="1-${SAMPLE_COUNT}%${TRIM_ARRAY_LIMIT}" \
    --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" trim "$CONFIG_PATH")

ALIGN_JOB=$(sbatch --parsable \
    --job-name=imp_align \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$ALIGN_CPUS" \
    --mem="$ALIGN_MEM" \
    --time="$ALIGN_TIME" \
    --array="1-${SAMPLE_COUNT}%${ALIGN_ARRAY_LIMIT}" \
    --dependency="afterok:$REFERENCE_JOB:$TRIM_JOB" \
    --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" align "$CONFIG_PATH")

COMPUTE_GL_JOB=$(sbatch --parsable \
    --job-name=imp_compute_gl \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$COMPUTE_GL_CPUS" \
    --mem="$COMPUTE_GL_MEM" \
    --time="$COMPUTE_GL_TIME" \
    --array="1-${COMPUTE_GL_COUNT}%${COMPUTE_GL_ARRAY_LIMIT}" \
    --dependency="afterok:$ALIGN_JOB" \
    --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" compute_gl "$CONFIG_PATH")

CONCAT_GL_JOB=$(sbatch --parsable \
    --job-name=imp_concat_gl \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$CONCAT_GL_CPUS" \
    --mem="$CONCAT_GL_MEM" \
    --time="$CONCAT_GL_TIME" \
    --array="1-${SAMPLE_COUNT}%${CONCAT_GL_ARRAY_LIMIT}" \
    --dependency="afterok:$COMPUTE_GL_JOB" \
    --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" concat_gl "$CONFIG_PATH")

MERGE_GL_JOB=$(sbatch --parsable \
    --job-name=imp_merge_gl \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$MERGE_GL_CPUS" \
    --mem="$MERGE_GL_MEM" \
    --time="$MERGE_GL_TIME" \
    --array="1-${CHROM_COUNT}%${MERGE_GL_ARRAY_LIMIT}" \
    --dependency="afterok:$CONCAT_GL_JOB" \
    --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" merge_gl "$CONFIG_PATH")

CHUNK_JOB=$(sbatch --parsable \
    --job-name=imp_chunk \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$CHUNK_CPUS" \
    --mem="$CHUNK_MEM" \
    --time="$CHUNK_TIME" \
    --array="1-${CHROM_COUNT}%${CHUNK_ARRAY_LIMIT}" \
    --output="$SLURM_OUTPUT_DIR/%x_%A_%a.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%A_%a.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/run_stage.sh" chunk "$CONFIG_PATH")

CONTINUE_JOB=$(sbatch --parsable \
    --job-name=imp_continue \
    --account="$SLURM_ACCOUNT" \
    --partition="$SLURM_PARTITION" \
    --cpus-per-task="$SUBMIT_CPUS" \
    --mem="$SUBMIT_MEM" \
    --time="$SUBMIT_TIME" \
    --dependency="afterok:$MERGE_GL_JOB:$CHUNK_JOB" \
    --output="$SLURM_OUTPUT_DIR/%x_%j.out" \
    --error="$SLURM_OUTPUT_DIR/%x_%j.err" \
    --chdir="$ROOT" \
    "$ROOT/scripts/submit_pipeline.sh" --continue-after-chunk "$CONFIG_PATH")

printf 'reference: %s\n' "$REFERENCE_JOB"
printf 'trim: %s\n' "$TRIM_JOB"
printf 'align: %s\n' "$ALIGN_JOB"
printf 'compute_gl: %s\n' "$COMPUTE_GL_JOB"
printf 'concat_gl: %s\n' "$CONCAT_GL_JOB"
printf 'merge_gl: %s\n' "$MERGE_GL_JOB"
printf 'chunk: %s\n' "$CHUNK_JOB"
printf 'continue: %s\n' "$CONTINUE_JOB"
