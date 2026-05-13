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

    def test_coordinates_wrapped_into_box(self, tmp_path: Path) -> None:
        system = _make_system()
        system.atoms = [
            LammpsAtom(1, 1, 1, 0.0, 1.0, 2.0, 3.0, "CG1", True),
            LammpsAtom(2, 1, 2, 0.0, -1.5, 55.0, 60.0, "CH2", False),
            LammpsAtom(3, 1, 2, 0.0, 7.0, -3.0, 9.0, "CH2", False),
        ]
        system.box = (50.0, 50.0, 50.0)
        system.bonds = []
        system.angles = []
        p = tmp_path / "test.data"
        write_lammps_data(system, p)
        content = p.read_text()
        lines = [
            ln
            for ln in content.split("\n")
            if ln.strip() and ln[0].isdigit() and len(ln.split()) == 7
        ]
        # atom 2: x=-1.5 wraps to 48.5, z=60 wraps to 10
        assert any("48.500000" in ln and "10.000000" in ln for ln in lines)
        # atom 3: y=-3 wraps to 47
        assert any("47.000000" in ln for ln in lines)

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

    def test_integration_before_fix_backmap(self, tmp_path: Path) -> None:
        """Integration fix must appear before fix backmap for correct ordering."""
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        integrate_pos = content.find("fix integrate")
        backmap_pos = content.find("fix bm all backmap")
        assert integrate_pos < backmap_pos, "fix integrate must come before fix backmap"

    def test_nose_hoover_nvt(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        settings.simulation.thermostat = "nose_hoover"
        settings.simulation.thermostat_tdamp = 0.1
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "fix integrate at_atoms nvt temp" in content
        assert "fix thermo at_atoms langevin" not in content

    def test_nose_hoover_npt(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        settings.simulation.thermostat = "nose_hoover"
        settings.simulation.ensemble = "npt"
        settings.simulation.pressure = 1.0
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "fix integrate at_atoms npt temp" in content
        assert "iso" in content

    def test_langevin_backward_compat(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        settings.simulation.thermostat = "langevin"
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "fix integrate at_atoms nve" in content
        assert "fix thermo at_atoms langevin" in content

    def test_fix_backmap_multi_cg_types(self, tmp_path: Path) -> None:
        system = _make_system()
        system.atom_types = [
            AtomTypeInfo(1, "A", 29.0, True),
            AtomTypeInfo(2, "B", 28.0, True),
            AtomTypeInfo(3, "CH3", 15.0, False),
            AtomTypeInfo(4, "CH2", 14.0, False),
        ]
        system.pair_types = [
            PairTypeInfo(1, 1, "cg"),
            PairTypeInfo(1, 2, "cg"),
            PairTypeInfo(2, 2, "cg"),
            PairTypeInfo(3, 3, "atomistic", sigma=3.75, epsilon=0.21),
            PairTypeInfo(3, 4, "none"),
            PairTypeInfo(4, 4, "atomistic", sigma=3.91, epsilon=0.12),
        ]
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "cg_type 1 2 " in content

    def test_default_backmap_without_cg_equil(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.test"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "# Backmapping: λ 0 → 1" in content
        assert "write_data test_hybrid.data" in content
        assert "# AT production" not in content
        assert "fix_modify bm active no" not in content
        assert "fix_modify bm active yes" in content

    def test_three_phases_when_equilibration_enabled(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        settings.simulation.equilibration_steps = 10000
        settings.simulation.production_steps = 10000
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


class TestRestartGeneration:
    def _make_restart_settings(self) -> Settings:
        s = _make_settings()
        s.simulation.restart_interval = 5000
        s.simulation.production_steps = 10000
        return s

    def test_no_restart_by_default(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "write_restart" not in content
        assert "restart " not in content
        assert "write_data test_hybrid.data" in content
        assert not (tmp_path / "in.backmap.setup").exists()

    def test_restart_commands_present(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "restart 5000 restart.backmap restart.backmap2" in content
        assert "write_restart restart.backmap" in content

    def test_sentinel_commands(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        content = p.read_text()
        assert "phase_1.done" in content
        assert "phase_2.done" in content
        assert "phase_3.done" not in content
        assert "write_data test_hybrid.data" in content

    def test_per_phase_scripts_generated(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        assert (tmp_path / "in.backmap.setup").exists()
        assert (tmp_path / "in.backmap.phase1").exists()
        assert (tmp_path / "in.backmap.phase2").exists()
        assert not (tmp_path / "in.backmap.phase3").exists()

    def test_restart_backmap_only_no_at_production_phase_file(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        settings.simulation.production_steps = 0
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        assert (tmp_path / "in.backmap.phase1").exists()
        assert not (tmp_path / "in.backmap.phase2").exists()
        p1 = (tmp_path / "in.backmap.phase1").read_text()
        assert "write_data test_hybrid.data" in p1

    def test_phase2_uses_read_restart_for_at_production(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        content = (tmp_path / "in.backmap.phase2").read_text()
        assert "read_restart restart.backmap" in content
        assert "include in.backmap.setup" in content
        assert "AT production" in content

    def test_phase1_backmapping_has_active_yes(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        content = (tmp_path / "in.backmap.phase1").read_text()
        assert "fix_modify bm active yes" in content

    def test_three_phase_restart_scripts_when_equilibration_enabled(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = self._make_restart_settings()
        settings.simulation.equilibration_steps = 500
        settings.simulation.production_steps = 10000
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        assert (tmp_path / "in.backmap.phase3").exists()
        p3 = (tmp_path / "in.backmap.phase3").read_text()
        assert "read_restart restart.backmap" in p3
        assert "Phase 3" in p3

    def test_no_per_phase_without_restart(self, tmp_path: Path) -> None:
        system = _make_system()
        settings = _make_settings()
        p = tmp_path / "in.backmap"
        write_lammps_input(system, settings, p, "test.data")
        assert not (tmp_path / "in.backmap.phase1").exists()
        assert not (tmp_path / "in.backmap.phase2").exists()
        assert not (tmp_path / "in.backmap.phase3").exists()
