#!/usr/bin/env bash
# One job of the queue (called by queue_runner.sh, or by hand):
#   queue_job.sh <queue-dir> <name> <source-dir> <command...>
# The job directory <queue-dir>/<name> is created once by copying <source-dir>, with the relative paths
# of the settings files made absolute (the copy no longer sits next to the shared force-field and table
# directories). The command runs in the job directory and must be safe to repeat; the job is finished
# when it exits 0 (marker .job.done) and given up when it fails (marker .job.failed, remove it to retry).
# Settings of the queue (LMP, NP, PREP, RESUMABLE, ...) come from <queue-dir>/queue.env.
set -u
qd="$(cd "$1" && pwd)"; name="$2"; src="$3"; shift 3
[ -f "$qd/queue.env" ] && set -a && . "$qd/queue.env" && set +a
d="$qd/$name"
log="$qd/jobs.log"
[ -e "$d/.job.done" ] && exit 0
[ -e "$d/.job.failed" ] && { echo "$name has failed before ($d/.job.failed), skipped $(date -u)" >> "$log"; exit 0; }
if [ ! -e "$d" ]; then
  cp -r "$src" "$d"
  python3 - "$src" "$d" <<'PY'
import re, sys
from pathlib import Path
src, dst = Path(sys.argv[1]).resolve(), Path(sys.argv[2])
for st in dst.glob("settings*.yaml"):
    t = st.read_text()
    t = re.sub(r"^(\s*(?:data_dir|tables_dir|forcefield_dir|bakery_xml):\s*)([^\s#]+)",
               lambda m: m.group(1) + (m.group(2) if m.group(2).startswith("/") else str((src / m.group(2)).resolve())),
               t, flags=re.MULTILINE)
    st.write_text(t)
PY
fi
echo "$name start $(date -u)" >> "$log"
( cd "$d" && "$@" >> run.out 2>&1 )
rc=$?
if [ "$rc" = 0 ]; then touch "$d/.job.done"; else echo "rc=$rc $(date -u)" > "$d/.job.failed"; fi
echo "$name rc=$rc $(date -u)" >> "$log"
exit 0
