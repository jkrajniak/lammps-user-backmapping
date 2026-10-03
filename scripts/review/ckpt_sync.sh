#!/usr/bin/env bash
# Copy a queue directory to a persistent destination so that a replacement machine can continue.
#   ckpt_sync.sh <queue-dir> <destination> [seconds between copies; 0 = once]
# <destination> is anything rsync takes (user@host:/path, or a mounted bucket directory). Restart files,
# logs, markers and result files go along; trajectories (*.dcd) and temporary restart probes do not
# (EXCLUDE adds more patterns, space separated). Nothing is ever deleted at the destination, and a source
# without queue.txt is not synced, so a replacement machine that starts with an empty queue directory
# cannot wipe the copy before it has been restored. To continue on a new machine:
#   rsync -a <destination>/ <queue-dir>/ && queue_runner.sh <queue-dir>
set -u
src="${1:?usage: ckpt_sync.sh <queue-dir> <destination> [seconds]}"; dst="${2:?destination}"; every="${3:-0}"
ex=(--exclude '*.dcd' --exclude '.probe.*' --exclude '.runner.lock')
for pat in ${EXCLUDE:-}; do ex+=(--exclude "$pat"); done
while true; do
  if [ -f "$src/queue.txt" ]; then
    rsync -a --partial "${ex[@]}" "$src"/ "$dst"/ 2> "$src/.sync.err" || echo "ckpt_sync.sh: rsync failed $(date -u)" >> "$src/jobs.log"
  else
    echo "ckpt_sync.sh: $src has no queue.txt, not syncing $(date -u)" >> "$src/jobs.log"
  fi
  [ "$every" = 0 ] && exit 0
  sleep "$every"
done
