"""Tests for scripts/review/network_tierc.py and the RDF core of rdf_vs_reference.py."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts" / "review"


@pytest.fixture(scope="module")
def nt():
    pytest.importorskip("scipy")  # review scripts only; not a dependency of the package
    sys.path.insert(0, str(SCRIPTS))
    return importlib.import_module("network_tierc")


def frame(rv, xyz_nm: np.ndarray, box_nm: float, step: int = 0):
    n = len(xyz_nm)
    return rv.Frame(
        step,
        np.full(3, box_nm * 10.0),
        np.arange(1, n + 1),
        np.ones(n, dtype=int),
        xyz_nm * 10.0,
    )


def test_ideal_gas_rdf_is_one(nt) -> None:
    rv = nt.rv
    rng = np.random.default_rng(1)
    frames = [frame(rv, rng.uniform(0, 5.0, size=(3000, 3)), 5.0, s) for s in range(4)]
    r, g = rv.rdf(frames, {1: "C"}, "C", "C", 1.0, 0.02, None)
    mask = (r > 0.3) & (r < 1.0)
    assert abs(g[mask].mean() - 1.0) < 0.01
    assert g[mask].std() < 0.05


def test_pair_counts_respect_periodic_boundary(nt) -> None:
    rv = nt.rv
    # two atoms 0.3 nm apart across the boundary of a 5 nm box
    fr = frame(rv, np.array([[0.1, 2.5, 2.5], [4.8, 2.5, 2.5]]), 5.0)
    group = rv.element_group({1: "C"}, "C")
    centers, counts, norm, _ = rv.rdf_per_frame([fr], group, group, 1.0, 0.002, None)
    assert counts.sum() == 2  # both ordered pairs
    assert centers[int(np.argmax(counts[0]))] == pytest.approx(0.3, abs=0.002)
    assert norm[0] == pytest.approx(2 * 2 / 125.0)


def test_exclusions_are_nrexcl_bonds_on_a_chain(nt) -> None:
    rv = nt.rv
    bonds = [(i, i + 1) for i in range(1, 6)]  # 6 atoms in a chain
    assert len(rv.excluded_pairs(bonds, 3)) == 5 + 4 + 3
    rng = np.random.default_rng(2)
    fr = frame(rv, rng.uniform(0, 2.0, size=(6, 3)), 5.0)
    group = rv.element_group({1: "C"}, "C")
    _, with_excl, _, _ = rv.rdf_per_frame(
        [fr], group, group, 6.0, 0.05, rv.excluded_pairs(bonds, 3)
    )
    _, without, _, _ = rv.rdf_per_frame([fr], group, group, 6.0, 0.05, None)
    assert without.sum() - with_excl.sum() == 2 * 12  # excluded pairs, counted in both directions


def test_residue_com_unwraps_across_the_box_and_weights_by_mass(nt, tmp_path: Path) -> None:
    rv = nt.rv
    top = tmp_path / "t.top"
    top.write_text(
        "[ atoms ]\n"
        "1 a 1 RES C1 1 0.0 12.0\n2 a 1 RES C2 2 0.0 12.0\n3 a 1 RES O1 3 0.0 16.0\n"
        "4 a 2 RES C1 4 0.0 12.0\n5 a 2 RES C2 5 0.0 12.0\n6 a 2 RES O1 6 0.0 16.0\n"
    )
    topology = nt.read_topology(top)
    # residue 1 straddles the periodic boundary of a 5 nm box, residue 2 does not
    xyz = np.array([[4.9, 1, 1], [0.1, 1, 1], [0.0, 1, 1], [2.0, 3, 3], [2.2, 3, 3], [2.4, 3, 3]])
    fr = frame(rv, xyz, 5.0)
    ids, com = nt.com_group(topology, "RES", "C?")(fr)
    assert list(ids) == [1, 2]
    assert com[0] == pytest.approx([0.0, 1.0, 1.0], abs=1e-9)  # (4.9 + 5.1) / 2 = 5.0 -> 0.0
    assert com[1] == pytest.approx([2.1, 3.0, 3.0], abs=1e-9)  # the oxygen is not in the selection
    with pytest.raises(ValueError, match="no atoms match"):
        nt.com_group(topology, "XXX", "C?")


def test_topology_must_match_the_data_file(nt, tmp_path: Path) -> None:
    rv = nt.rv
    top = tmp_path / "t.top"
    top.write_text("[ atoms ]\n1 a 1 R C1 1 0.0 12.011\n2 a 1 R O1 2 0.0 15.999\n")
    topology = nt.read_topology(top)
    fr = frame(rv, np.array([[1.0, 1, 1], [2.0, 1, 1]]), 5.0)
    fr.types = np.array([1, 2])
    nt.check_topology(topology, {1: "C", 2: "O"}, fr)
    with pytest.raises(ValueError, match="do not match"):
        nt.check_topology(topology, {1: "C", 2: "N"}, fr)


def test_first_peak_and_rms_ignore_the_bonded_range(nt) -> None:
    r = np.arange(0.0, 1.0, 0.01)
    g = np.where(r < 0.15, 50.0, 1.0)
    g[40] = 3.0
    assert nt.first_peak(r, g, 0.2) == (pytest.approx(0.4), 3.0)
    assert nt.rms_difference(r, g, g.copy(), 0.2) == 0.0
    assert nt.rms_difference(r, g + 0.5, g, 0.2) == pytest.approx(0.5)
