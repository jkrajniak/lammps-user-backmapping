"""Pydantic models for the YAML settings file."""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from pathlib import Path

import yaml
from pydantic import BaseModel, Field, field_validator, model_validator


class BeadDef(BaseModel):
    """A single CG bead within a molecule."""

    name: str
    type: str
    atoms: list[str]

    # Deferred Phase 3 fields
    atoms_by_degree: list[str] | None = None
    active_site: str | None = None
    remove: list[str] | None = None


class SourceFiles(BaseModel):
    """Source AT files for a molecule type."""

    coordinates: str | list[str]
    topology: str | list[str]

    @field_validator("coordinates", "topology")
    @classmethod
    def check_deferred_list(cls, v: str | list[str], info: object) -> str | list[str]:
        if isinstance(v, list):
            field_name = getattr(info, "field_name", "unknown")
            raise ValueError(
                f"Feature 'degree-dependent source files' in "
                f"'{field_name}' is not yet implemented (planned for Phase 3)"
            )
        return v


class MoleculeDef(BaseModel):
    """CG molecule definition with AT mapping."""

    name: str
    ident: str | None = None
    source: SourceFiles
    beads: list[BeadDef]

    # Deferred Phase 3 fields
    charge_management: dict[str, str] | None = None
    charge_map: dict[str, str] | None = None
    type_map: dict[str, str] | None = None


class CGSystem(BaseModel):
    """CG configuration files."""

    coordinates: str
    topology: str
    format: Literal["gromacs"] = "gromacs"
    predefined_active_sites: str | None = None


class CrossBond(BaseModel):
    """A cross-CG bond interaction."""

    params: str
    pairs: list[list[str]]
    table: str | None = None
    cg_bonded: bool = False


class CrossAngle(BaseModel):
    """A cross-CG angle interaction."""

    params: str
    triples: list[list[str]]
    table: str | None = None
    cg_bonded: bool = False


class CrossDihedral(BaseModel):
    """A cross-CG dihedral interaction."""

    params: str
    quadruples: list[list[str]]
    table: str | None = None
    cg_bonded: bool = False


class CrossInteractions(BaseModel):
    """Cross-CG bonded interactions (lambda-weighted)."""

    bonds: list[CrossBond] = Field(default_factory=list)
    angles: list[CrossAngle] = Field(default_factory=list)
    dihedrals: list[CrossDihedral] = Field(default_factory=list)
    pairs: list[list[str]] = Field(default_factory=list)


class SimulationParams(BaseModel):
    """Backmapping simulation parameters (in GROMACS units)."""

    alpha: float = 0.001
    initial_resolution: float = 0.0
    nonuniform_lambda: bool = False

    timestep: float = 0.001
    timestep_backmapping: float = 0.001

    equilibration_steps: int = 10000
    production_steps: int = 10000

    temperature: float = 300.0
    thermostat: Literal["langevin", "velocity_rescaling"] = "langevin"
    thermostat_gamma: float = 0.5
    thermostat_target: Literal["atomistic", "all", "cg_only"] = "atomistic"

    lj_cutoff: float = 1.2
    cg_cutoff: float = 1.4
    coulomb_cutoff: float = 0.9

    table_groups: list[str] = Field(default_factory=list)

    exclusion_nrexcl: int = 3

    energy_interval: int = 1000
    trajectory_interval: int = 1000

    rng_seed: int = -1

    # Deferred Phase 2 fields
    alpha2: float | None = None
    two_phase: bool = False
    second_phase_em: bool = False
    cap_force: float | None = None
    cap_force_ramp: float | None = None
    em_steps: int = 0
    em_gamma: float = 0.0001
    em_ftol: float = 10.0
    disable_angles: bool = False
    disable_dihedrals: bool = False
    coulomb_epsilon1: float = 1.0
    coulomb_epsilon2: float = 78.0

    @field_validator("alpha")
    @classmethod
    def alpha_positive(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("alpha must be positive")
        return v

    @field_validator("temperature")
    @classmethod
    def temperature_positive(cls, v: float) -> float:
        if v <= 0:
            raise ValueError("temperature must be positive")
        return v


class OutputConfig(BaseModel):
    """Output file configuration."""

    prefix: str = "system"
    format: Literal["lammps"] = "lammps"
    units: Literal["real"] = "real"


class Settings(BaseModel):
    """Top-level settings model."""

    molecules: list[MoleculeDef]
    cg_system: CGSystem
    cross_interactions: CrossInteractions = Field(default_factory=CrossInteractions)
    simulation: SimulationParams = Field(default_factory=SimulationParams)
    output: OutputConfig = Field(default_factory=OutputConfig)

    class Config:
        extra = "allow"

    @model_validator(mode="after")
    def check_deferred_features(self) -> Settings:
        for mol in self.molecules:
            for bead in mol.beads:
                if bead.atoms_by_degree is not None:
                    raise ValueError(
                        "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"
                    )
                if bead.remove is not None:
                    raise ValueError(
                        "Feature 'remove' (atom removal) is not yet implemented "
                        "(planned for Phase 3)"
                    )
            if mol.charge_management is not None:
                raise ValueError(
                    "Feature 'charge_management' is not yet implemented (planned for Phase 3)"
                )
        if self.simulation.two_phase:
            raise ValueError(
                "Feature 'two_phase' backmapping is not yet implemented (planned for Phase 2)"
            )
        return self


def load_settings(path: Path) -> Settings:
    """Load and validate a YAML settings file."""
    with open(path) as f:
        raw = yaml.safe_load(f)
    return Settings(**raw)
