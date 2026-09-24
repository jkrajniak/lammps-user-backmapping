"""Integration tests — verify backmap-prep generates correct output for each example."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import ClassVar

import pytest

from backmap_prep.cli import main
from backmap_prep.parsers import parse_top
from backmap_prep.schema import load_settings

EXAMPLES_DIR = Path(__file__).resolve().parents[2] / "examples"

EXAMPLES = [
    "dodecane",
    "dodecane-lammps-cg",
    "pe",
    "pe-lammps",
    "pe4",
    "pe_10",
    "pe_aa",
    "melamine",
]


@pytest.fixture(params=EXAMPLES)
def example_workdir(request, tmp_path: Path) -> Path:
    name = request.param
    src = EXAMPLES_DIR / name
    if not src.exists():
        pytest.skip(f"Example {name} not found at {src}")
    dst = tmp_path / name
    shutil.copytree(src, dst)
    ff_src = EXAMPLES_DIR / "forcefield"
    if ff_src.is_dir():
        shutil.copytree(ff_src, tmp_path / "forcefield", dirs_exist_ok=True)
    return dst


class TestExampleGeneration:
    def test_backmap_prep_succeeds(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        assert settings.exists(), f"Missing settings.yaml in {example_workdir}"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            result = main([str(settings)])
        finally:
            os.chdir(old_cwd)
        assert result == 0, f"backmap-prep failed for {example_workdir.name}"

    def test_data_file_generated(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            main([str(settings)])
        finally:
            os.chdir(old_cwd)
        data_files = list(example_workdir.glob("*.data"))
        assert len(data_files) >= 1, "Expected at least 1 .data file"

    def test_input_script_generated(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            main([str(settings)])
        finally:
            os.chdir(old_cwd)
        in_files = list(example_workdir.glob("in.*"))
        in_files = [f for f in in_files if not f.name.endswith(".yaml")]
        assert len(in_files) >= 1, f"Expected at least 1 input script, found {len(in_files)}"

    def test_data_file_has_atoms_and_bonds(self, example_workdir: Path) -> None:
        settings = example_workdir / "settings.yaml"
        old_cwd = os.getcwd()
        try:
            os.chdir(example_workdir)
            main([str(settings)])
        finally:
            os.chdir(old_cwd)
        data_files = list(example_workdir.glob("*.data"))
        content = data_files[0].read_text()
        assert "atoms" in content
        assert "bonds" in content
        assert "Masses" in content
        assert "Atoms" in content


class TestDeterministicOutput:
    @pytest.mark.parametrize("example_name", EXAMPLES)
    def test_identical_runs(self, example_name: str, tmp_path: Path) -> None:
        src = EXAMPLES_DIR / example_name
        if not src.exists():
            pytest.skip(f"Example {example_name} not found")

        ff_src = EXAMPLES_DIR / "forcefield"
        if ff_src.is_dir():
            shutil.copytree(ff_src, tmp_path / "forcefield", dirs_exist_ok=True)
        run1 = tmp_path / "run1"
        run2 = tmp_path / "run2"
        shutil.copytree(src, run1)
        shutil.copytree(src, run2)

        for workdir in (run1, run2):
            old_cwd = os.getcwd()
            try:
                os.chdir(workdir)
                main([str(workdir / "settings.yaml")])
            finally:
                os.chdir(old_cwd)

        for f1 in run1.glob("*.data"):
            f2 = run2 / f1.name
            assert f2.exists()
            assert f1.read_text() == f2.read_text(), (
                f"Non-deterministic output for {example_name}/{f1.name}"
            )

        for f1 in run1.glob("in.*"):
            if f1.name.endswith(".yaml"):
                continue
            f2 = run2 / f1.name
            assert f2.exists()
            assert f1.read_text() == f2.read_text(), (
                f"Non-deterministic output for {example_name}/{f1.name}"
            )


class TestLammpsNativeCgParity:
    """dodecane-lammps-cg must build to the same hybrid system as dodecane.

    dodecane-lammps-cg/dodecane_cg.data was derived (see its README) from a
    hybrid build of dodecane/'s current cg_conf.gro/topol_cg.top, so a fresh
    build of both examples should differ only in the cosmetic type-name
    convention (A/B vs the LAMMPS numeric type IDs "1"/"2") documented there.
    """

    _NAME_SUBS: ClassVar[dict[str, str]] = {
        "# A (CG)": "# 1 (CG)",
        "# B (CG)": "# 2 (CG)",
        "table_A_A": "table_1_1",
        "table_A_B": "table_1_2",
        "table_B_B": "table_2_2",
    }

    def _normalize(self, text: str) -> str:
        for old, new in self._NAME_SUBS.items():
            text = text.replace(old, new)
        return text

    def _build(self, name: str, tmp_path: Path) -> Path:
        src = EXAMPLES_DIR / name
        dst = tmp_path / name
        shutil.copytree(src, dst)
        old_cwd = os.getcwd()
        try:
            os.chdir(dst)
            main([str(dst / "settings.yaml")])
        finally:
            os.chdir(old_cwd)
        return dst

    def test_hybrid_data_matches_after_type_name_normalization(self, tmp_path: Path) -> None:
        gromacs_dir = self._build("dodecane", tmp_path)
        lammps_dir = self._build("dodecane-lammps-cg", tmp_path)

        gromacs_data = self._normalize((gromacs_dir / "dodecane.data").read_text())
        lammps_data = (lammps_dir / "dodecane.data").read_text()
        assert gromacs_data == lammps_data

    def test_input_script_matches_after_type_name_normalization(self, tmp_path: Path) -> None:
        gromacs_dir = self._build("dodecane", tmp_path)
        lammps_dir = self._build("dodecane-lammps-cg", tmp_path)

        gromacs_in = self._normalize((gromacs_dir / "in.dodecane").read_text())
        lammps_in = (lammps_dir / "in.dodecane").read_text()
        assert gromacs_in == lammps_in


def _split_data_sections(text: str) -> dict[str, list[str]]:
    """Split a LAMMPS data file into {section_name: [content lines]}, plus "header"."""
    sections: dict[str, list[str]] = {"header": []}
    current = "header"
    for line in text.splitlines():
        stripped = line.strip()
        if (
            stripped
            and not stripped[0].isdigit()
            and not any(kw in stripped for kw in ("xlo", "ylo", "zlo"))
        ):
            current = stripped.split("#", 1)[0].strip()
            sections.setdefault(current, [])
            continue
        if stripped:
            sections[current].append(stripped)
    return sections


def _atom_lines_without_image_flags(lines: list[str]) -> list[str]:
    """Drop trailing ix/iy/iz tokens (fields 8-10) from each Atoms line."""
    return [" ".join(line.split()[:7]) for line in lines]


class TestLammpsNativeAtFragmentParity:
    """pe-lammps must build to the same hybrid system as pe.

    pe-lammps/pe_at.data + in.pe_at were converted directly from pe/'s
    pe_single.gro + topol_aa.top (not extracted from a built hybrid — see
    openspec/changes/lammps-native-at-fragment-input/design.md), and
    pe-lammps/pe_cg.data from a fresh pe/ hybrid build via prepare_cg.py.
    A fresh build of both examples should differ only in the cosmetic
    type-name convention and per-atom PBC image-flag bookkeeping (positions
    are identical either way — see the design doc's Risks section for why
    the image-flag difference is a benign canonicalization artifact, not a
    correctness issue).
    """

    _NAME_SUBS: ClassVar[dict[str, str]] = {
        "# A (CG)": "# 1 (CG)",
        "# B (CG)": "# 2 (CG)",
        "# CH3": "# AT1",
        "# CH2": "# AT2",
        "table_A_A": "table_1_1",
        "table_A_B": "table_1_2",
        "table_B_B": "table_2_2",
    }

    def _normalize(self, text: str) -> str:
        for old, new in self._NAME_SUBS.items():
            text = text.replace(old, new)
        return text

    def _build(self, name: str, tmp_path: Path) -> Path:
        src = EXAMPLES_DIR / name
        dst = tmp_path / name
        shutil.copytree(src, dst)
        old_cwd = os.getcwd()
        try:
            os.chdir(dst)
            main([str(dst / "settings.yaml")])
        finally:
            os.chdir(old_cwd)
        return dst

    def test_hybrid_data_matches_ignoring_image_flags(self, tmp_path: Path) -> None:
        gromacs_dir = self._build("pe", tmp_path)
        lammps_dir = self._build("pe-lammps", tmp_path)

        gromacs_sections = _split_data_sections(
            self._normalize((gromacs_dir / "pe.data").read_text())
        )
        lammps_sections = _split_data_sections((lammps_dir / "pe.data").read_text())

        assert set(gromacs_sections) == set(lammps_sections)
        for name in gromacs_sections:
            if name == "Atoms":
                assert _atom_lines_without_image_flags(
                    gromacs_sections[name]
                ) == _atom_lines_without_image_flags(lammps_sections[name])
            else:
                assert gromacs_sections[name] == lammps_sections[name], f"section {name!r} differs"

    def test_input_script_matches_after_type_name_normalization(self, tmp_path: Path) -> None:
        gromacs_dir = self._build("pe", tmp_path)
        lammps_dir = self._build("pe-lammps", tmp_path)

        gromacs_in = self._normalize((gromacs_dir / "in.pe").read_text())
        lammps_in = (lammps_dir / "in.pe").read_text()
        assert gromacs_in == lammps_in


def _build_example(workdir: Path) -> Path:
    """Run backmap-prep in ``workdir`` and return the generated data file."""
    settings_path = workdir / "settings.yaml"
    old_cwd = os.getcwd()
    try:
        os.chdir(workdir)
        assert main([str(settings_path)]) == 0
    finally:
        os.chdir(old_cwd)
    prefix = load_settings(settings_path).output.prefix
    return workdir / f"{prefix}.data"


def _read_hybrid(
    path: Path,
) -> tuple[list[float], dict[int, dict], dict[int, float], dict[str, list[list[int]]]]:
    box: list[float] = []
    masses: dict[int, float] = {}
    cg_types: set[int] = set()
    atoms: dict[int, dict] = {}
    terms: dict[str, list[list[int]]] = {"Bonds": [], "Angles": [], "Dihedrals": []}
    section = None
    for line in path.read_text().splitlines():
        body = line.split("#")[0].strip()
        if not body:
            continue
        if line.strip().split("#")[0].strip() in (
            "Masses",
            "Atoms",
            "Bonds",
            "Angles",
            "Dihedrals",
            "Impropers",
            "Velocities",
        ):
            section = body
            continue
        parts = body.split()
        if len(parts) >= 4 and parts[2] in ("xlo", "ylo", "zlo"):
            box.append(float(parts[1]) - float(parts[0]))
        elif section == "Masses":
            masses[int(parts[0])] = float(parts[1])
            if "(CG)" in line:
                cg_types.add(int(parts[0]))
        elif section == "Atoms":
            atoms[int(parts[0])] = {
                "mol": int(parts[1]),
                "type": int(parts[2]),
                "x": [float(v) for v in parts[4:7]],
                "cg": int(parts[2]) in cg_types,
            }
        elif section in terms:
            terms[section].append([int(v) for v in parts[2:]])
    return box, atoms, masses, terms


class TestHybridInvariants:
    """Properties every generated hybrid system must have (OpenSpec unify-hybrid-engine)."""

    def test_every_bead_sits_on_its_fragment_com(self, example_workdir: Path) -> None:
        box, atoms, masses, _ = _read_hybrid(_build_example(example_workdir))
        by_mol: dict[int, list[dict]] = {}
        for atom in atoms.values():
            by_mol.setdefault(atom["mol"], []).append(atom)
        worst = 0.0
        for mol_atoms in by_mol.values():
            beads = [a for a in mol_atoms if a["cg"]]
            assert len(beads) == 1, "each molecule ID must hold one bead and its fragment"
            (bead,) = beads
            fragment = [a for a in mol_atoms if not a["cg"]]
            total = sum(masses[a["type"]] for a in fragment)
            shift = [0.0, 0.0, 0.0]
            for a in fragment:
                for k in range(3):
                    d = a["x"][k] - bead["x"][k]
                    d -= box[k] * round(d / box[k])
                    shift[k] += masses[a["type"]] * d / total
            worst = max(worst, sum(c * c for c in shift) ** 0.5)
        assert worst < 0.01, f"a bead is {worst:.3f} A from its fragment COM"

    def test_cg_bonded_terms_are_complete(self, example_workdir: Path) -> None:
        settings = load_settings(example_workdir / "settings.yaml")
        _, atoms, _, terms = _read_hybrid(_build_example(example_workdir))
        found = {
            name: sum(all(atoms[i]["cg"] for i in ids) for ids in entries)
            for name, entries in terms.items()
        }
        n_cg = sum(a["cg"] for a in atoms.values())
        n_beads = sum(len(mol.beads) for mol in settings.molecules)
        n_mol = n_cg // n_beads
        cg = settings.cg_system
        assert cg is not None
        if cg.format == "gromacs":
            assert cg.topology is not None
            top = parse_top(example_workdir / cg.topology, include_dirs=[example_workdir])
            per_mol = {
                "Bonds": sum(len(top.molecule_types[m].bonds) * n for m, n in top.molecules),
                "Angles": sum(len(top.molecule_types[m].angles) * n for m, n in top.molecules),
                "Dihedrals": sum(
                    len(top.molecule_types[m].dihedrals) * n for m, n in top.molecules
                ),
            }
            expected = per_mol
        else:
            ci = settings.cross_interactions
            expected = {
                "Bonds": n_mol * sum(len(e.pairs) for e in ci.bonds if e.cg_bonded),
                "Angles": n_mol * sum(len(e.triples) for e in ci.angles if e.cg_bonded),
                "Dihedrals": n_mol * sum(len(e.quadruples) for e in ci.dihedrals if e.cg_bonded),
            }
        assert found == expected
