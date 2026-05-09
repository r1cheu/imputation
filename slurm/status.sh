#!/usr/bin/env bash
# cluster-generic status checker. arg: SLURM job id. stdout: running|success|failed.
set -euo pipefail
jobid=$1

# Try sacct first (post-completion), fall back to squeue (in queue / running).
state=$(sacct -j "$jobid" --format=State --noheader --parsable2 2>/dev/null | head -1 | awk '{print $1}')
if [[ -z "$state" ]]; then
    state=$(squeue -j "$jobid" --noheader --format=%T 2>/dev/null | head -1)
fi

case "$state" in
    COMPLETED)                                              echo success ;;
    FAILED|CANCELLED|TIMEOUT|NODE_FAIL|OUT_OF_MEMORY|BOOT_FAIL|DEADLINE|PREEMPTED|REVOKED|SPECIAL_EXIT)
                                                            echo failed ;;
    PENDING|RUNNING|REQUEUED|RESIZING|SUSPENDED|CONFIGURING|COMPLETING|"")
                                                            echo running ;;
    *)                                                      echo running ;;
esac
