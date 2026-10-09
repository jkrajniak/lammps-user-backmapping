"""Fragment orientation on the CG beads (``prep.fragment_orientation``)."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import numpy as np
import pytest

from backmap_prep.cli import main
from backmap_prep.network.bakery.structures import rotation_aligning
from backmap_prep.schema import PrepConfig

EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples"


def _random_rotation(rng: np.random.Generator) -> np.ndarray:
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


def test_default_orientation_is_translate() -> None:
    assert PrepConfig().fragment_orientation == "translate"
    with pytest.raises(ValueError, match="fragment_orientation"):
        PrepConfig(fragment_orientation="rotate")  # type: ignore[arg-type]


def test_no_neighbors_gives_identity() -> None:
    assert np.allclose(rotation_aligning([], []), np.eye(3))


@pytest.mark.parametrize("n_vectors", [2, 3, 5])
def test_recovers_a_known_rotation(n_vectors: int) -> None:
    rng = np.random.default_rng(7)
    rot = _random_rotation(rng)
    template = rng.normal(size=(n_vectors, 3))
    target = template @ rot.T
    assert np.allclose(rotation_aligning(list(template), list(target)), rot, atol=1e-10)


def test_rotation_is_proper_for_a_mirrored_target() -> None:
    rng = np.random.default_rng(3)
    template = rng.normal(size=(4, 3))
    mirrored = template * np.array([-1.0, 1.0, 1.0])
    rot = rotation_aligning(list(template), list(mirrored))
    assert np.linalg.det(rot) == pytest.approx(1.0)
    assert np.allclose(rot @ rot.T, np.eye(3), atol=1e-12)


@pytest.mark.parametrize(
    ("template", "target"),
    [
        ([1.0, 0.0, 0.0], [0.0, 2.0, 0.0]),
        ([0.0, 0.0, 3.0], [0.0, 0.0, -1.0]),  # antiparallel
        ([1.0, 1.0, 0.0], [2.0, 2.0, 0.0]),  # parallel
    ],
)
def test_single_vector_alignment(template: list[float], target: list[float]) -> None:
    rot = rotation_aligning([np.array(template)], [np.array(target)])
    assert np.linalg.det(rot) == pytest.approx(1.0)
    got = rot @ (np.array(template) / np.linalg.norm(template))
    assert np.allclose(got, np.array(target) / np.linalg.norm(target), atol=1e-12)


def test_collinear_vectors_fall_back_to_single_vector_alignment() -> None:
    rot = rotation_aligning(
        [np.array([1.0, 0.0, 0.0]), np.array([-2.0, 0.0, 0.0])],
        [np.array([0.0, 1.0, 0.0]), np.array([0.0, -1.0, 0.0])],
    )
    assert np.allclose(rot @ np.array([1.0, 0.0, 0.0]), [0.0, 1.0, 0.0], atol=1e-12)


def _build_dodecane(tmp_path: Path, orientation: str) -> Path:
    src = EXAMPLES_DIR / "dodecane"
    if not src.exists():
        pytest.skip("dodecane example not found")
    run = tmp_path / orientation
    shutil.copytree(src, run)
    shutil.copytree(EXAMPLES_DIR / "forcefield", tmp_path / "forcefield", dirs_exist_ok=True)
    settings = run / "settings.yaml"
    text = settings.read_text()
    assert "\nprep:" not in text or "fragment_orientation" not in text
    text += f"\nprep:\n  fragment_orientation: {orientation}\n" if "\nprep:" not in text else ""
    settings.write_text(text)
    old = os.getcwd()
    try:
        os.chdir(run)
        main([str(settings)])
    finally:
        os.chdir(old)
    return run


def _gro_rows(path: Path) -> tuple[list[str], np.ndarray]:
    lines = path.read_text().splitlines()
    n = int(lines[1])
    names = [ln[10:15].strip() for ln in lines[2 : 2 + n]]
    xyz = np.array(
        [[float(ln[20:28]), float(ln[28:36]), float(ln[36:44])] for ln in lines[2 : 2 + n]]
    )
    return names, xyz


def test_align_rotates_fragments_rigidly(tmp_path: Path) -> None:
    """CG beads stay put, fragments keep their size and are oriented differently."""
    names_t, translate = _gro_rows(_build_dodecane(tmp_path, "translate") / "hyb_conf.gro")
    names_a, align = _gro_rows(_build_dodecane(tmp_path, "align") / "hyb_conf.gro")
    assert names_t == names_a
    is_bead = np.array([not n.startswith("C") for n in names_t])
    assert is_bead.any()
    assert np.allclose(translate[is_bead], align[is_bead])  # beads do not move
    assert not np.allclose(translate[~is_bead], align[~is_bead])  # atoms are rotated
    bead_of = np.cumsum(is_bead) - 1
    for atom in np.flatnonzero(~is_bead):
        bead = np.flatnonzero(is_bead)[bead_of[atom]]
        d_translate = np.linalg.norm(translate[atom] - translate[bead])
        d_align = np.linalg.norm(align[atom] - align[bead])
        assert d_align == pytest.approx(d_translate, abs=2e-3)  # .gro precision is 1e-3 nm
