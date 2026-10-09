"""Merge the per-segment output files of a resumed LAMMPS run.

A run that was interrupted and continued (resumable_lmp.sh) leaves one file per segment, named
``<name>.s<N>[.<ext>]``. This merges every such family into ``<name>[.<ext>]``:

- LAMMPS custom dump text files (``ITEM: TIMESTEP`` blocks): frames sorted by timestep, a timestep
  that occurs in several segments is taken from the latest segment.
- ``fix ave/time ... mode vector`` and ``fix ave/histo ... mode vector`` files (``timestep nrows``
  header lines followed by nrows data lines): blocks sorted and de-duplicated the same way.
- other formats (e.g. DCD) are left as segments, with a note.

Segments are never deleted. Run in the job directory (or pass it), safe to repeat.

    python merge_segments.py [directory]
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from pathlib import Path

SEGMENT = re.compile(r"^(?P<stem>.+?)\.s(?P<seg>\d+)(?P<ext>\.[A-Za-z0-9]+)?$")
SKIP_FAMILIES = ("log", "rst")


def families(directory: Path) -> dict[str, dict[int, Path]]:
    found: dict[str, dict[int, Path]] = defaultdict(dict)
    for path in sorted(directory.iterdir()):
        match = SEGMENT.match(path.name)
        if not match or not path.is_file() or path.name.startswith("."):
            continue
        stem, seg, ext = match.group("stem"), int(match.group("seg")), match.group("ext") or ""
        if stem.split(".")[0] in SKIP_FAMILIES or ext in (".a", ".b"):
            continue
        found[stem + ext][seg] = path
    return found


def sniff(path: Path) -> str:
    with path.open("rb") as fh:
        head = fh.read(64)
    if head.startswith(b"ITEM: TIMESTEP"):
        return "dump"
    if head.startswith(b"#"):
        return "avetime"
    return "other"


def dump_frames(path: Path) -> dict[int, list[str]]:
    frames: dict[int, list[str]] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].startswith("ITEM: TIMESTEP"):
            step = int(lines[i + 1])
            n_atoms = int(lines[i + 3])
            end = i + 9 + n_atoms
            if end > len(lines):
                break  # frame cut by the interruption
            frames[step] = lines[i:end]
            i = end
        else:
            i += 1
    return frames


def avetime_blocks(path: Path) -> tuple[list[str], dict[int, list[str]]]:
    """Comment lines before the first block, and the blocks keyed by timestep.

    A block is a header line whose first two fields are the timestep and the number of data rows
    (ave/time vector: ``step nrows``; ave/histo: ``step nbins total missing min max``), followed by
    that many rows. The line after a block is the next header.
    """
    header: list[str] = []
    blocks: dict[int, list[str]] = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        s = lines[i].strip()
        if not s or s.startswith("#"):
            if not blocks:
                header.append(lines[i])
            i += 1
            continue
        parts = s.split()
        step, nrows = int(parts[0]), int(parts[1])
        if i + 1 + nrows > len(lines):
            break  # block cut by the interruption
        blocks[step] = lines[i : i + 1 + nrows]
        i += 1 + nrows
    return header, blocks


def merge_family(name: str, segs: dict[int, Path], directory: Path) -> str:
    kind = sniff(segs[min(segs)])
    ordered = [segs[k] for k in sorted(segs)]
    if kind == "dump":
        merged: dict[int, list[str]] = {}
        for path in ordered:
            merged.update(dump_frames(path))
        out = [line for step in sorted(merged) for line in merged[step]]
        (directory / name).write_text("\n".join(out) + "\n")
        return f"{name}: {len(merged)} frames from {len(ordered)} segment(s)"
    if kind == "avetime":
        header: list[str] = []
        merged_blocks: dict[int, list[str]] = {}
        for path in ordered:
            head, blocks = avetime_blocks(path)
            header = header or head
            merged_blocks.update(blocks)
        out = header + [line for step in sorted(merged_blocks) for line in merged_blocks[step]]
        (directory / name).write_text("\n".join(out) + "\n")
        return f"{name}: {len(merged_blocks)} blocks from {len(ordered)} segment(s)"
    return f"{name}: format not mergeable, {len(ordered)} segment file(s) left as they are"


def main() -> int:
    directory = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    for name, segs in families(directory).items():
        print(merge_family(name, segs, directory))
    return 0


if __name__ == "__main__":
    sys.exit(main())
