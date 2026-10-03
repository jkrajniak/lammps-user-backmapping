"""Tests for scripts/review/merge_segments.py (merging the files of a resumed LAMMPS run)."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "review" / "merge_segments.py"


@pytest.fixture(scope="module")
def merge():
    spec = importlib.util.spec_from_file_location("merge_segments", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def frame(step: int, x: float) -> str:
    return (
        f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n2\nITEM: BOX BOUNDS pp pp pp\n"
        "0 10\n0 10\n0 10\nITEM: ATOMS id mol type x y z\n"
        f"1 1 1 {x} 0 0\n2 1 1 {x + 1} 0 0\n"
    )


def block(step: int, value: float, ncols: int = 2) -> str:
    return (
        f"{step} 2\n" + "".join(f"{i} {value} {i * value}\n" for i in range(1, 3))
        if ncols == 2
        else (f"{step} 2 100 0 0.0 180.0\n1 1.0 {value} 0.1\n2 3.0 {value} 0.2\n")
    )


def test_dump_frames_sorted_deduplicated_latest_segment_wins(
    tmp_path: Path, merge, monkeypatch: pytest.MonkeyPatch
) -> None:
    (tmp_path / "dump.at_prod.s1").write_text(frame(100, 1.0) + frame(200, 2.0) + frame(300, 3.0))
    (tmp_path / "dump.at_prod.s2").write_text(frame(300, 33.0) + frame(400, 4.0))
    monkeypatch.setattr(sys, "argv", ["merge_segments.py", str(tmp_path)])
    assert merge.main() == 0
    frames = merge.dump_frames(tmp_path / "dump.at_prod")
    assert sorted(frames) == [100, 200, 300, 400]
    assert "33.0" in "\n".join(frames[300])  # the later segment wins


def test_frame_cut_by_the_interruption_is_dropped(tmp_path: Path, merge) -> None:
    text = frame(100, 1.0) + frame(200, 2.0)
    (tmp_path / "dump.at_prod.s1").write_text(text[: len(text) - 25])
    frames = merge.dump_frames(tmp_path / "dump.at_prod.s1")
    assert sorted(frames) == [100]


@pytest.mark.parametrize("ncols", [2, 6])
def test_avetime_and_histogram_blocks_merge(
    tmp_path: Path, merge, ncols: int, monkeypatch: pytest.MonkeyPatch
) -> None:
    head = "# Time-averaged data for fix rdf_out\n# TimeStep Number-of-rows\n"
    (tmp_path / "rdf_backmap.s1.dat").write_text(
        head + block(20, 1.0, ncols) + block(40, 2.0, ncols)
    )
    (tmp_path / "rdf_backmap.s2.dat").write_text(
        head + block(40, 22.0, ncols) + block(60, 3.0, ncols)
    )
    monkeypatch.setattr(sys, "argv", ["merge_segments.py", str(tmp_path)])
    assert merge.main() == 0
    header, blocks = merge.avetime_blocks(tmp_path / "rdf_backmap.dat")
    assert sorted(blocks) == [20, 40, 60]
    assert header[0].startswith("# Time-averaged")
    assert "22.0" in "\n".join(blocks[40])


def test_logs_and_restart_files_are_not_merged(tmp_path: Path, merge) -> None:
    for name in ("log.at.s1.lammps", "log.at.s2.lammps", "rst.s1.a", "rst.s1.b"):
        (tmp_path / name).write_text("x")
    assert merge.families(tmp_path) == {}
