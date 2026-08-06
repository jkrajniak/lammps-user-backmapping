"""Pydantic models for the YAML settings file."""

from __future__ import annotations

from pathlib import Path
from typing import Literal

import yaml
from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator


class DegreeSourceFile(BaseModel):
    """One degree-dependent AT source file (coordinates or topology)."""

    file: str
    molecule_degree: int | str = 0
    when: str | None = None


class AtomsByDegree(BaseModel):
    """AT atom list for a CG bead at a given bonding degree."""

    degree: int
    atoms: list[str]
    molecule_degree: str | None = None
    active_site: str | None = None
    charge_map: list[str | float] | None = None
    type_map: list[str] | None = None


class ChargeTransfer(BaseModel):
    """Charge redistribution rule when a bond forms at an active site."""

    when: str
    from_atom: str
    to_atoms: str


class ChargeManagement(BaseModel):
    """Charge equilibration settings for reactive networks."""

    equilibrate: bool = False
    transfers: list[ChargeTransfer] = Field(default_factory=list)


class ActiveSiteSpec(BaseModel):
    """Reactive site on a bead (v2 mapping)."""

    atom: str
    max_degree: int


class BeadMappingEntry(BaseModel):
    """Degree-dependent bead atom list with optional delta resolution (v2)."""

    degree: int
    molecule_degree: str | None = None
    atoms: list[str] | None = None
    base: int | None = None
    add: list[str] = Field(default_factory=list)
    remove: list[str] = Field(default_factory=list)
    active_sites: list[ActiveSiteSpec] = Field(default_factory=list)
    charge_map: list[str | float] | None = None
    type_map: list[str] | None = None

    @model_validator(mode="after")
    def check_atoms_or_delta(self) -> BeadMappingEntry:
        if self.atoms is None and self.base is None:
            raise ValueError(f"mapping degree {self.degree} must define 'atoms' or 'base'")
        return self


class BeadDef(BaseModel):
    """A single CG bead within a molecule."""

    name: str
    type: str
    atoms: list[str] = Field(default_factory=list)
    atoms_by_degree: list[AtomsByDegree] | None = None
    mapping: list[BeadMappingEntry] | None = None
    active_site: str | None = None
    remove: list[dict[str, str | list[str]]] | None = None

    @model_validator(mode="after")
    def check_atoms_or_degree_mapping(self) -> BeadDef:
        if not self.atoms and not self.atoms_by_degree and not self.mapping:
            raise ValueError(
                f"bead '{self.name}' must define 'atoms', 'atoms_by_degree', or 'mapping'"
            )
        return self


def _normalize_source_files(
    v: str | list[str] | list[DegreeSourceFile] | list[dict[str, object]],
) -> str | list[DegreeSourceFile]:
    if isinstance(v, str):
        return v
    if not v:
        raise ValueError("source file list must not be empty")
    first = v[0]
    if isinstance(first, DegreeSourceFile):
        return v  # type: ignore[return-value]
    if isinstance(first, dict):
        return [DegreeSourceFile(**entry) for entry in v]  # type: ignore[arg-type]
    if isinstance(first, str):
        raise ValueError(
            "use DegreeSourceFile entries (file + molecule_degree) for multi-file sources"
        )
    raise ValueError(f"unsupported source file entry type: {type(first)!r}")


class UnifiedSourceRow(BaseModel):
    """v2: coordinates and topology for one molecule degree."""

    coordinates: str
    topology: str
    molecule_degree: int | str = 0
    when: str | None = None


class SourceFiles(BaseModel):
    """Source AT files for a molecule type."""

    coordinates: str | list[DegreeSourceFile]
    topology: str | list[DegreeSourceFile]

    @field_validator("coordinates", "topology", mode="before")
    @classmethod
    def normalize_source(cls, v: object) -> str | list[DegreeSourceFile]:
        if isinstance(v, (str, list)):
            return _normalize_source_files(v)
        raise ValueError("coordinates and topology must be a path string or a list of file entries")

    @classmethod
    def from_unified_rows(cls, rows: list[UnifiedSourceRow]) -> SourceFiles:
        coordinates = [
            DegreeSourceFile(
                file=row.coordinates,
                molecule_degree=row.molecule_degree,
                when=row.when,
            )
            for row in rows
        ]
        topology = [
            DegreeSourceFile(
                file=row.topology,
                molecule_degree=row.molecule_degree,
                when=row.when,
            )
            for row in rows
        ]
        return cls(coordinates=coordinates, topology=topology)


class MoleculeDef(BaseModel):
    """CG molecule definition with AT mapping."""

    name: str
    ident: str | None = None
    source: SourceFiles
    beads: list[BeadDef]

    charge_management: ChargeManagement | None = None
    charge_map: dict[str, list[str | float]] | None = None
    type_map: dict[str, list[str]] | None = None

    @model_validator(mode="before")
    @classmethod
    def normalize_v2_source(cls, data: object) -> object:
        if not isinstance(data, dict):
            return data
        source = data.get("source")
        if isinstance(source, list):
            rows = [UnifiedSourceRow(**row) if isinstance(row, dict) else row for row in source]
            data = {**data, "source": SourceFiles.from_unified_rows(rows)}
        return data


class CGSystem(BaseModel):
    """CG configuration files."""

    coordinates: str
    topology: str
    format: Literal["gromacs"] = "gromacs"
    predefined_active_sites: str | None = None


class CrossBond(BaseModel):
    """A cross-CG bond interaction."""

    params: str | None = None
    comment: str | None = None
    pairs: list[list[str]]
    table: str | None = None
    cg_bonded: bool = False

    @model_validator(mode="after")
    def default_params(self) -> CrossBond:
        if self.params is None:
            self.params = f"; {self.comment}" if self.comment else ";"
        return self


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
    angles_file: str | None = Field(
        default=None,
        description="External YAML file with CrossAngle entries (v2)",
    )
    dihedrals_file: str | None = Field(
        default=None,
        description="External YAML file with CrossDihedral entries (v2)",
    )


class HybridOutputConfig(BaseModel):
    """GROMACS hybrid output paths (network prep)."""

    coordinates: str = "hyb_conf.gro"
    topology: str = "hyb_topol.top"
    topology_includes: list[str] = Field(default_factory=list)
    includes: list[str] = Field(default_factory=list)
    exclusion: int = 3
    molecule_type_name: str = "HYB"
    system_name: str = "HYB"

    @model_validator(mode="before")
    @classmethod
    def normalize_hybrid_aliases(cls, data: object) -> object:
        if not isinstance(data, dict):
            return data
        if "includes" in data and "topology_includes" not in data:
            data = {**data, "topology_includes": data["includes"]}
        mol_type = data.get("molecule_type")
        if isinstance(mol_type, dict):
            data = {
                **data,
                "molecule_type_name": mol_type.get("name", "HYB"),
                "exclusion": mol_type.get("exclusion", data.get("exclusion", 3)),
            }
        return data

    @model_validator(mode="after")
    def merge_includes(self) -> HybridOutputConfig:
        if self.includes and not self.topology_includes:
            self.topology_includes = list(self.includes)
        elif self.includes:
            merged = list(self.topology_includes)
            for entry in self.includes:
                if entry not in merged:
                    merged.append(entry)
            self.topology_includes = merged
        return self


class PrepConfig(BaseModel):
    """Which preparation engine to use."""

    engine: Literal["linear", "network"] = "linear"
    bakery_xml: str | None = Field(
        default=None,
        description="Bakery XML settings for network prep (Phase 3 bridge)",
    )
    data_dir: str | None = Field(
        default=None,
        description="Directory for relative asset paths (resolved from YAML location)",
    )
    tables_dir: str | None = Field(
        default=None,
        description="Directory for CG tabulated potentials (.xvg), resolved from YAML location",
    )
    forcefield_dir: str | None = Field(
        default=None,
        description="Directory for OPLS/forcefield .itp includes (or set GMXDATA)",
    )
    allow_no_bonds: bool = False
    chain_rng_seed: int | None = Field(
        default=None,
        description="RNG seed for random chain selection in multi-chain AT sources",
    )


class SimulationParams(BaseModel):
    """Backmapping simulation parameters (in GROMACS units)."""

    alpha: float = 0.001
    initial_resolution: float = 0.0
    nonuniform_lambda: bool = False

    timestep: float = 0.001
    timestep_backmapping: float = 0.001

    equilibration_steps: int = Field(
        default=0,
        ge=0,
        description=(
            "In-hybrid equilibration steps with λ frozen (fix_modify active no; "
            "CG-AT coupling stays on). "
            "Use 0 when the CG melt was equilibrated before building the hybrid (recommended)."
        ),
    )
    omp_threads: int = Field(
        default=0,
        ge=0,
        description=(
            "OpenMP thread count for generated CG equilibration scripts "
            "(writes package omp + suffix omp; 0 disables)."
        ),
    )
    production_steps: int = Field(
        default=0,
        ge=0,
        description=(
            "Optional extra MD steps after λ reaches 1 in the same input file. "
            "Use 0 for backmapping-only output (write_data …_hybrid.data); run pure "
            "AT production (e.g. RDF) in a separate input using the extracted "
            "AT system."
        ),
    )

    temperature: float = 300.0
    ensemble: Literal["nvt", "npt"] = "nvt"
    thermostat: Literal["langevin", "nose_hoover", "velocity_rescaling"] = "langevin"
    thermostat_gamma: float = 0.5
    thermostat_tdamp: float = 0.1  # ps — Nosé-Hoover damping time
    thermostat_target: Literal["atomistic", "all", "cg_only"] = "atomistic"

    pressure: float = 1.0  # bar (used when ensemble = npt)
    barostat_pdamp: float = 1.0  # ps — barostat damping time

    lj_cutoff: float = 1.2
    cg_cutoff: float = 1.4
    coulomb_cutoff: float = 0.9

    table_groups: list[str] = Field(default_factory=list)

    exclusion_nrexcl: int = 3

    energy_interval: int = 1000
    trajectory_interval: int = 1000

    rng_seed: int = -1

    restart_interval: int | None = None

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

    @field_validator("restart_interval")
    @classmethod
    def restart_interval_positive(cls, v: int | None) -> int | None:
        if v is not None and v <= 0:
            raise ValueError("restart_interval must be a positive integer")
        return v


class OutputConfig(BaseModel):
    """Output file configuration."""

    prefix: str = "system"
    format: Literal["lammps"] = "lammps"
    units: Literal["real"] = "real"


class Settings(BaseModel):
    """Top-level settings model."""

    version: int | None = None
    prep: PrepConfig = Field(default_factory=PrepConfig)
    molecules: list[MoleculeDef] = Field(default_factory=list)
    cg_system: CGSystem | None = None
    hybrid: HybridOutputConfig | None = None
    cross_interactions: CrossInteractions = Field(default_factory=CrossInteractions)
    simulation: SimulationParams = Field(default_factory=SimulationParams)
    output: OutputConfig = Field(default_factory=OutputConfig)

    model_config = ConfigDict(extra="allow")

    @model_validator(mode="after")
    def check_engine_requirements(self) -> Settings:
        if self.prep.engine == "linear":
            if not self.molecules:
                raise ValueError("linear prep requires at least one molecule definition")
            if self.cg_system is None:
                raise ValueError("linear prep requires cg_system")
        elif self.prep.engine == "network":
            if self.prep.bakery_xml is None:
                if not self.molecules:
                    raise ValueError(
                        "network prep requires prep.bakery_xml or full molecule definitions"
                    )
                if self.cg_system is None:
                    raise ValueError("network v2 prep requires cg_system")
                if self.hybrid is None:
                    raise ValueError("network v2 prep requires hybrid output config")
        if self.simulation.two_phase:
            raise ValueError(
                "Feature 'two_phase' backmapping is not yet implemented (planned for Phase 2)"
            )
        return self


def _merge_angles_file(settings: Settings, yaml_dir: Path) -> Settings:
    """Load external cross-angle definitions referenced by angles_file."""
    angles_file = settings.cross_interactions.angles_file
    if not angles_file:
        return settings
    path = Path(angles_file)
    if not path.is_absolute():
        path = (yaml_dir / path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"cross_interactions.angles_file not found: {path}")
    with open(path) as f:
        raw_angles = yaml.safe_load(f)
    if not isinstance(raw_angles, list):
        raise ValueError(f"angles_file must contain a YAML list: {path}")
    extra = [CrossAngle(**entry) if isinstance(entry, dict) else entry for entry in raw_angles]
    merged = list(settings.cross_interactions.angles) + extra
    cross = settings.cross_interactions.model_copy(update={"angles": merged, "angles_file": None})
    return settings.model_copy(update={"cross_interactions": cross})


def _merge_dihedrals_file(settings: Settings, yaml_dir: Path) -> Settings:
    """Load external cross-dihedral definitions referenced by dihedrals_file."""
    dihedrals_file = settings.cross_interactions.dihedrals_file
    if not dihedrals_file:
        return settings
    path = Path(dihedrals_file)
    if not path.is_absolute():
        path = (yaml_dir / path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"cross_interactions.dihedrals_file not found: {path}")
    with open(path) as f:
        raw_dihedrals = yaml.safe_load(f)
    if not isinstance(raw_dihedrals, list):
        raise ValueError(f"dihedrals_file must contain a YAML list: {path}")
    extra = [
        CrossDihedral(**entry) if isinstance(entry, dict) else entry for entry in raw_dihedrals
    ]
    merged = list(settings.cross_interactions.dihedrals) + extra
    cross = settings.cross_interactions.model_copy(
        update={"dihedrals": merged, "dihedrals_file": None}
    )
    return settings.model_copy(update={"cross_interactions": cross})


def resolve_data_dir(settings_path: Path, settings: Settings) -> Path:
    """Return the working directory for relative network asset paths."""
    if settings.prep.data_dir:
        return (settings_path.parent / settings.prep.data_dir).resolve()
    return settings_path.parent.resolve()


def resolve_tables_dir(settings_path: Path, settings: Settings) -> Path | None:
    """Return directory containing CG .xvg tables, if configured."""
    if not settings.prep.tables_dir:
        return None
    return (settings_path.parent / settings.prep.tables_dir).resolve()


def resolve_forcefield_dir(settings_path: Path, settings: Settings) -> Path | None:
    """Return bundled OPLS/GROMACS forcefield directory, if configured."""
    if settings.prep.forcefield_dir:
        return (settings_path.parent / settings.prep.forcefield_dir).resolve()
    for candidate in (
        settings_path.parent / "forcefield",
        settings_path.parent.parent / "forcefield",
    ):
        if (candidate / "oplsaa.ff" / "forcefield.itp").is_file():
            return candidate.resolve()
    return None


def resolve_bakery_xml(settings_path: Path, settings: Settings) -> Path:
    """Resolve a bakery ``settings.xml`` path from ``prep.bakery_xml`` or a
    sibling ``<settings_stem>.xml`` file next to the YAML settings."""
    if settings.prep.bakery_xml:
        xml_path = (settings_path.parent / settings.prep.bakery_xml).resolve()
        if not xml_path.is_file():
            raise FileNotFoundError(f"bakery XML not found: {xml_path}")
        return xml_path
    sibling = settings_path.with_suffix(".xml")
    if sibling.is_file():
        return sibling.resolve()
    raise ValueError(
        f"network prep requires prep.bakery_xml in settings.yaml or a sibling file {sibling.name}"
    )


def load_settings(path: Path) -> Settings:
    """Load and validate a YAML settings file."""
    with open(path) as f:
        raw = yaml.safe_load(f)
    settings = Settings(**raw)
    settings = _merge_angles_file(settings, path.parent)
    return _merge_dihedrals_file(settings, path.parent)
