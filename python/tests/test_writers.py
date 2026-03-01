"""Tests for backmap_prep.writers — LAMMPS data and input file generation."""

from __future__ import annotations

from typing import TYPE_CHECKING

from backmap_prep.builder import (
    AngleTypeInfo,
    AtomTypeInfo,
    BondTypeInfo,
    LammpsAngle,
    LammpsAtom,
    LammpsBond,
    PairTypeInfo,
    System,
)
from backmap_prep.schema import Settings
from backmap_prep.writers import write_lammps_data, write_lammps_input

if TYPE_CHECKING:
    from pathlib import Path


def _make_system() -> System:
    """Build a small test system."""
    return System(
        atoms=[
            LammpsAtom(1, 1, 1, 0.0, 1.0, 2.0, 3.0, "CG1", True),
            LammpsAtom(2, 1, 2, 0.0, 4.0, 5.0, 6.0, "CH2", False),
            LammpsAtom(3, 1, 2, 0.1, 7.0, 8.0, 9.0, "CH2", False),
        ],
        bonds=[
            LammpsBond(1, 1, 2, 3),
        ],
        angles=[
            LammpsAngle(1, 1, 1, 2, 3),
        ],
        atom_types=[
            AtomTypeInfo(1, "CG1", 72.0, True),
            AtomTypeInfo(2, "CH2", 14.0, False, sigma=3.95, epsilon=0.382),
        ],
        bond_types=[
            BondTypeInfo(1, "harmonic", "", [100.0, 1.54]),
        ],
        angle_types=[
            AngleTypeInfo(1, "harmonic", "", [50.0, 111.0]),
        ],
        pair_types=[
            PairTypeInfo(1, 1, "cg"),
            PairTypeInfo(1, 2, "none"),
            PairTypeInfo(2, 2, "atomistic", sigma=3.95, epsilon=0.382),
        ],
        box=(50.0, 50.0, 50.0),
        cg_type_id=1,
    )


def _make_settings() -> Settings:
    return Settings(
        molecules=[
            {
                "name": "M",
                "source": {"coordinates": "a.gro", "topology": "a.top"},
                "beads": [{"name": "B", "type": "C", "atoms": ["X"]}],
            }
        ],
        cg_system={"coordinates": "c.gro", "topology": "c.top"},
    )


class TestWriteLammpsData:
    def test_creates_file(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        assert p.exists()

    def test_header_counts(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "3 atoms" in content
        assert "1 bonds" in content
        assert "1 angles" in content
        assert "2 atom types" in content
        assert "1 bond types" in content
        assert "1 angle types" in content

    def test_box_dimensions(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "0.0 50.000000 xlo xhi" in content
        assert "0.0 50.000000 ylo yhi" in content
        assert "0.0 50.000000 zlo zhi" in content

    def test_masses_section(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "Masses" in content
        assert "72.000000" in content
        assert "14.000000" in content
        assert "(CG)" in content

    def test_atoms_section(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "Atoms # full" in content
        in_atoms = False
        atom_lines = []
        for line in content.split("\n"):
            if line.startswith("Atoms"):
                in_atoms = True
                continue
            if in_atoms and line.strip() == "":
                if atom_lines:
                    break
                continue
            if in_atoms:
                atom_lines.append(line)
        assert len(atom_lines) == 3

    def test_bonds_section(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "Bonds" in content
        assert "1 1 2 3" in content

    def test_angles_section(self, tmp_path: Path) -> None:
        system = _make_system()
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "Angles" in content
        assert "1 1 1 2 3" in content

    def test_no_bonds_no_section(self, tmp_path: Path) -> None:
        system = _make_system()
        system.bonds = []
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        assert "0 bonds" in content
        assert "\nBonds\n" not in content


class TestWriteLammpsInput:
    def test_creates_file(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        assert p.exists()

    def test_header(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "units real" in content
        assert "atom_style full" in content
        assert "boundary p p p" in content

    def test_read_data(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "mydata.data")
        content = p.read_text()
        assert "read_data mydata.data" in content

    def test_pair_style(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "pair_style backmap" in content
        assert "pair_coeff 1 1 cg" in content
        assert "pair_coeff 1 2 none" in content
        assert "pair_coeff 2 2 atomistic" in content

    def test_fix_backmap(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "fix bm all backmap" in content
        assert "cg_type 1" in content

    def test_three_phases(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "Phase 1" in content
        assert "Phase 2" in content
        assert "Phase 3" in content
        assert "fix_modify bm active no" in content
        assert "fix_modify bm active yes" in content

    def test_special_bonds_nrexcl3(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0" in content

    def test_hybrid_bond_style(self, tmp_path: Path) -> None:
        system = _make_system()
        system.bond_types.append(BondTypeInfo(2, "backmap/harmonic", "at", [50.0, 1.5]))
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "bond_style hybrid" in content
