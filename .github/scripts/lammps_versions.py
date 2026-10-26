"""Matrix of LAMMPS versions for the integration workflow, as JSON for GITHUB_OUTPUT.

The oldest supported stable release is pinned; the latest stable release is resolved from the
tags of the LAMMPS repository at run time, so the matrix cannot go stale when LAMMPS releases.
Scheduled and manual runs add the develop branch as a non-blocking canary for API changes.

    python lammps_versions.py <github event name>
"""

from __future__ import annotations

import json
import re
import subprocess
import sys
from datetime import datetime

REPO = "https://github.com/lammps/lammps.git"
OLDEST = "stable_22Jul2025_update2"
TAG = re.compile(r"^stable_(\d{1,2}[A-Z][a-z]{2}\d{4})(?:_update(\d+))?$")


def latest_stable() -> str:
    out = subprocess.run(
        ["git", "ls-remote", "--tags", "--refs", REPO, "stable_*"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout
    best: tuple[datetime, int, str] | None = None
    for line in out.splitlines():
        name = line.split("refs/tags/")[-1]
        m = TAG.match(name)
        if not m:
            continue
        key = (datetime.strptime(m.group(1), "%d%b%Y"), int(m.group(2) or 0), name)
        if best is None or key > best:
            best = key
    if best is None:
        sys.exit("no stable_* tag found")
    return best[2]


def main() -> None:
    event = sys.argv[1] if len(sys.argv) > 1 else ""
    refs = {OLDEST: True, latest_stable(): True}
    if event in ("schedule", "workflow_dispatch"):
        refs["develop"] = False
    include = [{"ref": ref, "required": required} for ref, required in refs.items()]
    print(f"matrix={json.dumps({'include': include})}")


if __name__ == "__main__":
    main()
