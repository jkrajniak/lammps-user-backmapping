#!/usr/bin/env bash
# Queue supervisor for spot or preemptible machines. Safe to start again at any time (after a reboot,
# a lost machine, or by hand): only one instance runs per queue directory, finished jobs are skipped,
# an unfinished job continues from its newest checkpoint (the job command is resumable).
#   queue_runner.sh <queue-dir>
# <queue-dir> holds
#   queue.txt   one job per line:  name|source-dir|command ...      (lines starting with # are ignored)
#   queue.env   KEY=VALUE settings for the jobs: LMP, NP, PREP, RESUMABLE, ... (see queue_job.sh)
# SLOTS (default 2) jobs run at once. With SYNC_DEST=<rsync destination> a background loop copies the
# queue directory there every SYNC_EVERY seconds (default 600; ckpt_sync.sh), so that a replacement
# machine can continue from the copy.
set -u
qd="$(cd "${1:?usage: queue_runner.sh <queue-dir>}" && pwd)"
HERE="$(cd "$(dirname "$0")" && pwd)"
SLOTS="${SLOTS:-2}"
exec 9> "$qd/.runner.lock"
flock -n 9 || { echo "queue_runner.sh: another runner holds $qd/.runner.lock" >&2; exit 0; }
[ -f "$qd/queue.env" ] && set -a && . "$qd/queue.env" && set +a
if [ -n "${SYNC_DEST:-}" ]; then
  "$HERE/ckpt_sync.sh" "$qd" "$SYNC_DEST" "${SYNC_EVERY:-600}" &
  sync_pid=$!
  trap 'kill "$sync_pid" 2> /dev/null' EXIT
fi
echo "queue_runner.sh: started $(date -u +%FT%TZ), $SLOTS slot(s)" >> "$qd/jobs.log"
while true; do
  pending="$(grep -v '^[[:space:]]*#' "$qd/queue.txt" | grep -v '^[[:space:]]*$' | while IFS='|' read -r name src cmd; do
    [ -e "$qd/$name/.job.done" ] || [ -e "$qd/$name/.job.failed" ] || printf '%s|%s|%s\n' "$name" "$src" "$cmd"
  done)"
  [ -n "$pending" ] || break
  printf '%s\n' "$pending" | xargs -d '\n' -P "$SLOTS" -I{} bash -c '
    IFS="|" read -r name src cmd <<< "$1"
    exec "$2/queue_job.sh" "$3" "$name" "$src" bash -c "$cmd"
  ' _ {} "$HERE" "$qd"
  # one pass over the pending jobs ended: jobs that failed are marked, the others finished
done
echo "queue_runner.sh: queue finished $(date -u +%FT%TZ)" >> "$qd/jobs.log"
