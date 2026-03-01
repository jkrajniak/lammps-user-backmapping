"""Tests for backmap_prep.schema — Pydantic settings models and validation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
import yaml

if TYPE_CHECKING:
    from pathlib import Path

from backmap_prep.schema import (
    BeadDef,
    CGSystem,
    CrossInteractions,
    MoleculeDef,
    OutputConfig,
    Settings,
    SimulationParams,
    SourceFiles,
    load_settings,
)


class TestSourceFiles:
    def test_string_values_accepted(self) -> None:
        sf = SourceFiles(coordinates="at.gro", topology="at.top")
        assert sf.coordinates == "at.gro"
        assert sf.topology == "at.top"

    def test_list_coordinates_rejected(self) -> None:
        with pytest.raises(ValueError, match="not yet implemented"):
            SourceFiles(coordinates=["a.gro", "b.gro"], topology="at.top")

    def test_list_topology_rejected(self) -> None:
        with pytest.raises(ValueError, match="not yet implemented"):
            SourceFiles(coordinates="at.gro", topology=["a.top", "b.top"])


class TestSimulationParams:
    def test_defaults(self) -> None:
        sp = SimulationParams()
        assert sp.alpha == 0.001
        assert sp.temperature == 300.0
        assert sp.thermostat == "langevin"

    def test_alpha_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="alpha must be positive"):
            SimulationParams(alpha=0.0)

    def test_negative_alpha_rejected(self) -> None:
        with pytest.raises(ValueError, match="alpha must be positive"):
            SimulationParams(alpha=-0.5)

    def test_temperature_must_be_positive(self) -> None:
        with pytest.raises(ValueError, match="temperature must be positive"):
            SimulationParams(temperature=0.0)

    def test_custom_values(self) -> None:
        sp = SimulationParams(
            alpha=0.01,
            temperature=450.0,
            thermostat="velocity_rescaling",
        )
        assert sp.alpha == 0.01
        assert sp.temperature == 450.0


class TestOutputConfig:
    def test_defaults(self) -> None:
        oc = OutputConfig()
        assert oc.prefix == "system"
        assert oc.format == "lammps"
        assert oc.units == "real"


class TestBeadDef:
    def test_basic_bead(self) -> None:
        bd = BeadDef(name="B1", type="CG1", atoms=["C1", "C2"])
        assert bd.name == "B1"
        assert bd.atoms_by_degree is None
        assert bd.remove is None


class TestCGSystem:
    def test_basic(self) -> None:
        cg = CGSystem(coordinates="cg.gro", topology="cg.top")
        assert cg.format == "gromacs"
        assert cg.predefined_active_sites is None


class TestSettings:
    def test_minimal_settings(self, minimal_settings_dict: dict) -> None:
        s = Settings(**minimal_settings_dict)
        assert len(s.molecules) == 1
        assert s.molecules[0].name == "TestMol"
        assert isinstance(s.cross_interactions, CrossInteractions)

    def test_deferred_atoms_by_degree_rejected(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["molecules"][0]["beads"][0]["atoms_by_degree"] = ["C1"]
        with pytest.raises(ValueError, match=r"atoms_by_degree.*not yet implemented"):
            Settings(**minimal_settings_dict)

    def test_deferred_remove_rejected(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["molecules"][0]["beads"][0]["remove"] = ["H1"]
        with pytest.raises(ValueError, match=r"remove.*not yet implemented"):
            Settings(**minimal_settings_dict)

    def test_deferred_charge_management_rejected(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["molecules"][0]["charge_management"] = {"mode": "split"}
        with pytest.raises(ValueError, match=r"charge_management.*not yet implemented"):
            Settings(**minimal_settings_dict)

    def test_deferred_two_phase_rejected(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["simulation"] = {"two_phase": True}
        with pytest.raises(ValueError, match=r"two_phase.*not yet implemented"):
            Settings(**minimal_settings_dict)


class TestLoadSettings:
    def test_load_from_yaml(self, tmp_path: Path, minimal_settings_dict: dict) -> None:
        p = tmp_path / "settings.yaml"
        p.write_text(yaml.dump(minimal_settings_dict))
        s = load_settings(p)
        assert s.molecules[0].name == "TestMol"
        assert s.simulation.alpha == 0.001

    def test_load_with_custom_simulation(self, tmp_path: Path, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["simulation"] = {
            "alpha": 0.01,
            "temperature": 400.0,
        }
        p = tmp_path / "settings.yaml"
        p.write_text(yaml.dump(minimal_settings_dict))
        s = load_settings(p)
        assert s.simulation.alpha == 0.01
        assert s.simulation.temperature == 400.0

    def test_missing_file_raises(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            load_settings(tmp_path / "nonexistent.yaml")


class TestMoleculeDef:
    def test_ident_defaults_to_none(self) -> None:
        md = MoleculeDef(
            name="Mol",
            source=SourceFiles(coordinates="a.gro", topology="a.top"),
            beads=[BeadDef(name="B1", type="CG1", atoms=["C1"])],
        )
        assert md.ident is None

    def test_with_ident(self) -> None:
        md = MoleculeDef(
            name="Mol",
            ident="DOD",
            source=SourceFiles(coordinates="a.gro", topology="a.top"),
            beads=[BeadDef(name="B1", type="CG1", atoms=["C1"])],
        )
        assert md.ident == "DOD"
