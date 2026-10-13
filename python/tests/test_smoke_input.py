"""``scripts/ci/make_smoke_input.py``: shortened input for the release smoke test."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "ci" / "make_smoke_input.py"

spec = importlib.util.spec_from_file_location("make_smoke_input", SCRIPT)
assert spec is not None
assert spec.loader is not None
smoke = importlib.util.module_from_spec(spec)
spec.loader.exec_module(smoke)

INPUT = """include sys.backmap.lmp
thermo 1000
minimize 1.0e-4 1.0e-6 1000 10000
run 10000
run 20000
run 10000
run 5000
run 5000
write_data sys_hybrid.data
run 10000
write_data sys_final.data
"""


def test_every_run_is_shortened_and_markers_added() -> None:
    text = smoke.shorten(INPUT, "sys.backmap.lmp", "sys.backmap.smoke.lmp", 10)
    runs = [ln for ln in text.splitlines() if ln.startswith("run ")]
    assert runs == ["run 1000", "run 2000", "run 1000", "run 500", "run 500", "run 1000"]
    lines = text.splitlines()
    assert lines[lines.index("run 1000") + 1] == 'print "EQ done"'
    assert lines[lines.index("run 2000") + 1] == 'print "RAMP done"'
    assert (
        lines[lines.index("run 1000", lines.index("run 500") + 1) + 1] == 'print "PRODUCTION done"'
    )
    assert "include sys.backmap.smoke.lmp" in lines
    assert "thermo 100" in lines


def test_short_runs_keep_a_minimum() -> None:
    text = smoke.shorten(INPUT, "sys.backmap.lmp", "x.lmp", 1000)
    assert min(int(ln.split()[1]) for ln in text.splitlines() if ln.startswith("run ")) == 50


def test_alpha_is_scaled_so_the_ramp_reaches_one() -> None:
    out = smoke.scale_alpha("fix bm all backmap cg_type 1 alpha 0.0001 lambda0 0.0\n", 10)
    assert "alpha 0.001 " in out


def test_inputs_without_the_expected_shape_are_refused() -> None:
    with pytest.raises(ValueError, match="three run"):
        smoke.shorten("run 10\n", "a", "b", 10)
    with pytest.raises(ValueError, match="alpha"):
        smoke.scale_alpha("fix bm all backmap\n", 10)


def test_the_example_input_has_the_stages_the_script_expects() -> None:
    example = Path(__file__).resolve().parents[2] / "examples" / "pet" / "large" / "in.pet"
    if not example.exists():
        pytest.skip("PET example input not found")
    text = smoke.shorten(example.read_text(), "pet.backmap.lmp", "pet.backmap.smoke.lmp", 10)
    assert text.count('print "') == 3
