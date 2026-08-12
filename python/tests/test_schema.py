"""Tests for backmap_prep.schema — Pydantic settings models and validation."""

from __future__ import annotations

from typing import TYPE_CHECKING

import pytest
import yaml
from pydantic import ValidationError

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

    def test_degree_dependent_sources_accepted(self) -> None:
        sf = SourceFiles(
            coordinates=[
                {"file": "a.gro", "molecule_degree": 0},
                {"file": "b.gro", "molecule_degree": 1, "when": "A1"},
            ],
            topology=[
                {"file": "a.itp", "molecule_degree": 0},
                {"file": "b.itp", "molecule_degree": 1, "when": "A1"},
            ],
        )
        assert isinstance(sf.coordinates, list)
        assert len(sf.coordinates) == 2
        assert sf.coordinates[1].when == "A1"

    def test_lammps_format(self) -> None:
        sf = SourceFiles(format="lammps", data="at.data", input_script="in.at")
        assert sf.format == "lammps"
        assert sf.data == "at.data"
        assert sf.input_script == "in.at"

    def test_gromacs_missing_coordinates_raises(self) -> None:
        with pytest.raises(ValidationError, match="coordinates"):
            SourceFiles(topology="at.top")

    def test_gromacs_missing_topology_raises(self) -> None:
        with pytest.raises(ValidationError, match="topology"):
            SourceFiles(coordinates="at.gro")

    def test_lammps_missing_data_raises(self) -> None:
        with pytest.raises(ValidationError, match="'data'"):
            SourceFiles(format="lammps", input_script="in.at")

    def test_lammps_missing_input_script_raises(self) -> None:
        with pytest.raises(ValidationError, match="'input_script'"):
            SourceFiles(format="lammps", data="at.data")

    def test_lammps_rejects_degree_dependent_coordinates(self) -> None:
        with pytest.raises(ValidationError, match="degree-dependent"):
            SourceFiles(
                format="lammps",
                data="at.data",
                input_script="in.at",
                coordinates=[{"file": "a.gro", "molecule_degree": 0}],
            )


class TestPrepConfig:
    def test_network_engine_with_bakery_xml(self) -> None:
        s = Settings(
            prep={"engine": "network", "bakery_xml": "settings.xml"},
        )
        assert s.prep.engine == "network"
        assert s.prep.bakery_xml == "settings.xml"

    def test_linear_requires_molecules(self) -> None:
        with pytest.raises(ValueError, match="linear prep requires"):
            Settings(prep={"engine": "linear"})


class TestSimulationParams:
    def test_defaults(self) -> None:
        sp = SimulationParams()
        assert sp.alpha == 0.001
        assert sp.temperature == 300.0
        assert sp.thermostat == "langevin"
        assert sp.ensemble == "nvt"
        assert sp.equilibration_steps == 0
        assert sp.production_steps == 0

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

    def test_nose_hoover_thermostat(self) -> None:
        sp = SimulationParams(thermostat="nose_hoover", thermostat_tdamp=0.2)
        assert sp.thermostat == "nose_hoover"
        assert sp.thermostat_tdamp == 0.2

    def test_npt_ensemble(self) -> None:
        sp = SimulationParams(ensemble="npt", pressure=2.0, barostat_pdamp=0.5)
        assert sp.ensemble == "npt"
        assert sp.pressure == 2.0
        assert sp.barostat_pdamp == 0.5


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

    def test_lammps_format(self) -> None:
        cg = CGSystem(format="lammps", data="cg_system.data")
        assert cg.format == "lammps"
        assert cg.data == "cg_system.data"

    def test_gromacs_missing_coordinates_raises(self) -> None:
        with pytest.raises(ValidationError, match="coordinates"):
            CGSystem(topology="cg.top")

    def test_gromacs_missing_topology_raises(self) -> None:
        with pytest.raises(ValidationError, match="topology"):
            CGSystem(coordinates="cg.gro")

    def test_lammps_missing_data_raises(self) -> None:
        with pytest.raises(ValidationError, match="'data'"):
            CGSystem(format="lammps")


class TestSettings:
    def test_minimal_settings(self, minimal_settings_dict: dict) -> None:
        s = Settings(**minimal_settings_dict)
        assert len(s.molecules) == 1
        assert s.molecules[0].name == "TestMol"
        assert isinstance(s.cross_interactions, CrossInteractions)

    def test_network_engine_rejects_lammps_at_fragment(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["prep"] = {"engine": "network", "bakery_xml": "settings.xml"}
        minimal_settings_dict["molecules"][0]["source"] = {
            "format": "lammps",
            "data": "at.data",
            "input_script": "in.at",
        }
        with pytest.raises(ValidationError, match="network"):
            Settings(**minimal_settings_dict)

    def test_atoms_by_degree_accepted(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["molecules"][0]["beads"][0]["atoms_by_degree"] = [
            {"degree": 1, "molecule_degree": "0", "atoms": ["1:MOL:C1"]},
        ]
        s = Settings(**minimal_settings_dict)
        bead = s.molecules[0].beads[0]
        assert bead.atoms_by_degree is not None
        assert bead.atoms_by_degree[0].degree == 1

    def test_remove_list_accepted(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["molecules"][0]["beads"][0]["remove"] = [
            {"active_site": "MOL:ATOM", "atoms": ["1:MOL:H1"]},
        ]
        s = Settings(**minimal_settings_dict)
        assert s.molecules[0].beads[0].remove is not None

    def test_charge_management_accepted(self, minimal_settings_dict: dict) -> None:
        minimal_settings_dict["molecules"][0]["charge_management"] = {
            "equilibrate": True,
            "transfers": [
                {
                    "when": "IPD:N1:2",
                    "from_atom": "IPD:H8",
                    "to_atoms": "EPO:C1#H25",
                },
            ],
        }
        s = Settings(**minimal_settings_dict)
        assert s.molecules[0].charge_management is not None
        assert s.molecules[0].charge_management.equilibrate is True

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


class TestRestartInterval:
    def test_default_is_none(self) -> None:
        sp = SimulationParams()
        assert sp.restart_interval is None

    def test_positive_value_accepted(self) -> None:
        sp = SimulationParams(restart_interval=5000)
        assert sp.restart_interval == 5000

    def test_zero_rejected(self) -> None:
        with pytest.raises(ValueError, match="restart_interval must be a positive integer"):
            SimulationParams(restart_interval=0)

    def test_negative_rejected(self) -> None:
        with pytest.raises(ValueError, match="restart_interval must be a positive integer"):
            SimulationParams(restart_interval=-100)


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
