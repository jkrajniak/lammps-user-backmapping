#!/usr/bin/env bash
# Run a LAMMPS input that follows the resume protocol until it completes, continuing from
# the newest valid restart file after a crash, a kill, or a machine that went away.
#
#   LMP=/path/to/lmp NP=8 resumable_lmp.sh <tag> -in <input> -log <log> [-var name value ...]
#
# Run it in the job directory; call it again with the same arguments after an interruption.
# The input must accept `-var seg N`, `-var resume 0|1` and `-var resume_file F`, write its
# restart files as rst.s<seg>.a / rst.s<seg>.b, and name every per-segment output
# <name>.s<seg>[.<ext>] (see examples/common/in.at_reference_protocol).
#
# State in the job directory (all small):
#   .<tag>.seg      number of the last started segment
#   .<tag>.done     the run completed (the script then only re-merges and exits 0)
#   .<tag>.fail     the run failed repeatedly from the same checkpoint (exit 1)
#   rst.s*          restart files; only the two newest valid ones are kept
#   <log>.s<seg>    one log per segment (the -log name with .s<seg> inserted before the extension)
#
# Environment: LMP (required), NP (default 8), MAX_SAME_CKPT_FAILURES (default 2),
#              MPIRUN (default mpirun), MERGE (default: merge_segments.py next to this script).
set -u
LMP="${LMP:?set LMP to a LAMMPS binary with the packages the input needs}"
NP="${NP:-8}"
MPIRUN="${MPIRUN:-mpirun}"
MAX_FAIL="${MAX_SAME_CKPT_FAILURES:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"
MERGE="${MERGE:-python3 $HERE/merge_segments.py}"

tag="${1:?usage: resumable_lmp.sh <tag> -in <input> -log <log> [-var ...]}"
shift
args=("$@")
log_base=""
for ((i = 0; i < ${#args[@]}; i++)); do
  [ "${args[i]}" = "-log" ] && log_base="${args[i + 1]}"
done
[ -n "$log_base" ] || { echo "resumable_lmp.sh: -log <name> is required" >&2; exit 2; }

seg_file=".$tag.seg"; done_file=".$tag.done"; fail_file=".$tag.fail"

# step of a restart file, or empty when it is missing or unreadable (a half-written file)
restart_step() {
  [ -s "$1" ] || return 1
  printf 'read_restart %s\nprint "RSTSTEP $(step)"\n' "$1" > ".probe.$tag.in"
  "$LMP" -in ".probe.$tag.in" -log none -screen /dev/stdout 2> /dev/null | sed -n 's/^RSTSTEP //p' | head -1
  rm -f ".probe.$tag.in"
}

# newest valid restart file in $PWD; prints "<step> <file>" and deletes all but the two newest
newest_restart() {
  local f s list=""
  for f in rst.s*.a rst.s*.b; do
    [ -e "$f" ] || continue
    s="$(restart_step "$f")"
    [ -n "$s" ] && list="$list$s $f"$'\n'
  done
  [ -n "$list" ] || return 1
  list="$(printf '%s' "$list" | sort -k1,1nr)"
  printf '%s\n' "$list" | sed -n '3,$p' | while read -r _ old; do rm -f "$old"; done
  printf '%s\n' "$list" | sed -n '1p'
}

if [ -e "$done_file" ]; then
  $MERGE . || true
  exit 0
fi
[ -e "$fail_file" ] && { echo "resumable_lmp.sh: $tag previously failed ($fail_file); remove it to retry" >&2; exit 1; }

last_failed_step=""; same_fail=0
while true; do
  seg=$(( $(cat "$seg_file" 2> /dev/null || echo 0) + 1 ))
  echo "$seg" > "$seg_file"
  resume_args=()
  ckpt_step=""
  if best="$(newest_restart)"; then
    ckpt_step="${best%% *}"
    resume_args=(-var resume 1 -var resume_file "${best#* }")
    echo "resumable_lmp.sh: $tag segment $seg resumes from ${best#* } (step $ckpt_step) at $(date -u +%FT%TZ)"
  else
    echo "resumable_lmp.sh: $tag segment $seg starts fresh at $(date -u +%FT%TZ)"
  fi
  # per-segment log name: insert .s<seg> before the extension
  ext="${log_base##*.}"; stem="${log_base%.*}"
  [ "$ext" = "$log_base" ] && seg_log="$log_base.s$seg" || seg_log="$stem.s$seg.$ext"
  new_args=()
  skip=0
  for a in "${args[@]}"; do
    if [ "$skip" = 1 ]; then new_args+=("$seg_log"); skip=0; continue; fi
    new_args+=("$a")
    [ "$a" = "-log" ] && skip=1
  done
  "$MPIRUN" -np "$NP" "$LMP" "${new_args[@]}" -screen none -var seg "$seg" "${resume_args[@]}"
  rc=$?
  if [ "$rc" = 0 ]; then
    touch "$done_file"
    $MERGE . || echo "resumable_lmp.sh: merge failed (segments are kept)" >&2
    echo "resumable_lmp.sh: $tag finished in $seg segment(s)"
    exit 0
  fi
  echo "resumable_lmp.sh: $tag segment $seg exited with $rc" >&2
  if [ "${ckpt_step:-none}" = "${last_failed_step:-none2}" ]; then same_fail=$((same_fail + 1)); else same_fail=1; fi
  last_failed_step="${ckpt_step:-none}"
  newer="$(newest_restart)" && [ "${newer%% *}" != "${ckpt_step:-x}" ] && same_fail=0 && last_failed_step=""
  if [ "$same_fail" -ge "$MAX_FAIL" ]; then
    echo "resumable_lmp.sh: $tag failed $same_fail times from the same checkpoint; giving up" >&2
    echo "segment $seg rc $rc checkpoint ${ckpt_step:-none}" > "$fail_file"
    exit 1
  fi
done
