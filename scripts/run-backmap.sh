#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Restart-aware entrypoint for backmapping simulations.
Detects existing restart/sentinel files and resumes from the correct phase.

Options:
  -np N         Number of MPI processes (default: 1)
  -in FILE      Master input script name (default: in.backmap)
  -h, --help    Show this help message

The script checks for restart.backmap / restart.backmap2 and sentinel files
to determine whether to start fresh or resume.  Per-phase scripts (if
present) define the layout:

  * \${INPUT}.phase3 exists  → three segments (CG equil → backmap → AT production)
  * \${INPUT}.phase2 exists  → two segments (either CG equil → backmap, or
                               backmap → AT production, depending on settings)
  * only \${INPUT}.phase1   → backmapping only (writes *_hybrid.data)

Sentinels: phase_1.done, phase_2.done, phase_3.done as applicable.
EOF
    exit "${1:-0}"
}

NP=1
INPUT="in.backmap"

while [[ $# -gt 0 ]]; do
    case "$1" in
        -np)
            NP="$2"
            shift 2
            ;;
        -in)
            INPUT="$2"
            shift 2
            ;;
        -h|--help)
            usage 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage 1
            ;;
    esac
done

run_lmp() {
    local script="$1"
    if [[ "$NP" -gt 1 ]]; then
        mpirun --allow-run-as-root -np "$NP" lmp -in "$script"
    else
        lmp -in "$script"
    fi
}

HAS_RESTART=false
if [[ -f restart.backmap ]] || [[ -f restart.backmap2 ]]; then
    HAS_RESTART=true
fi

LAYOUT="master"
if [[ -f "${INPUT}.phase3" ]]; then
    LAYOUT="three"
elif [[ -f "${INPUT}.phase2" ]]; then
    LAYOUT="two"
elif [[ -f "${INPUT}.phase1" ]]; then
    LAYOUT="one"
fi

all_complete() {
    case "$LAYOUT" in
        three)
            [[ -f phase_3.done ]]
            ;;
        two)
            [[ -f phase_2.done ]]
            ;;
        one)
            [[ -f phase_1.done ]]
            ;;
        *)
            return 1
            ;;
    esac
}

if all_complete; then
    echo "All phases complete — nothing to do."
    exit 0
fi

if [[ "$HAS_RESTART" == "false" ]]; then
    echo "No restart file found — starting fresh."
    run_lmp "$INPUT"
    exit $?
fi

# Restart mode: resume the correct per-phase script
case "$LAYOUT" in
    three)
        if [[ -f phase_2.done ]]; then
            echo "Resuming phase 3 (AT production) from restart file."
            run_lmp "${INPUT}.phase3"
        elif [[ -f phase_1.done ]]; then
            echo "Resuming phase 2 (backmapping) from restart file."
            run_lmp "${INPUT}.phase2"
            echo "Continuing to phase 3..."
            run_lmp "${INPUT}.phase3"
        else
            echo "Restart file found but phase 1 not complete — restarting phase 1."
            run_lmp "${INPUT}.phase1"
            echo "Continuing to phase 2..."
            run_lmp "${INPUT}.phase2"
            echo "Continuing to phase 3..."
            run_lmp "${INPUT}.phase3"
        fi
        ;;
    two)
        if [[ -f phase_1.done ]]; then
            echo "Resuming phase 2 from restart file."
            run_lmp "${INPUT}.phase2"
        else
            echo "Restart file found but phase 1 not complete — restarting phase 1."
            run_lmp "${INPUT}.phase1"
            echo "Continuing to phase 2..."
            run_lmp "${INPUT}.phase2"
        fi
        ;;
    one)
        echo "Resuming phase 1 (backmapping only) from restart file."
        run_lmp "${INPUT}.phase1"
        ;;
    *)
        echo "Resuming with master input ${INPUT}."
        run_lmp "$INPUT"
        ;;
esac
