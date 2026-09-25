"""Data classes for a LAMMPS hybrid system.

The system itself is built by the single hybrid engine
(``network.api.build_network_lammps``) and written by ``writers``.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class LammpsAtom:
    atom_id: int
    mol_id: int
    type_id: int
    charge: float
    x: float  # Angstrom
    y: float  # Angstrom
    z: float  # Angstrom
    type_name: str = ""
    is_cg: bool = False
    ix: int = 0
    iy: int = 0
    iz: int = 0


@dataclass
class LammpsBond:
    bond_id: int
    type_id: int
    i: int
    j: int


@dataclass
class LammpsAngle:
    angle_id: int
    type_id: int
    i: int
    j: int
    k: int


@dataclass
class LammpsDihedral:
    dihedral_id: int
    type_id: int
    i: int
    j: int
    k: int
    l: int


@dataclass
class LammpsImproper:
    improper_id: int
    type_id: int
    i: int
    j: int
    k: int
    l: int


@dataclass
class LammpsCrossPair:
    i: int
    j: int
    sigma: float
    epsilon: float
    keyword: str = "at"
    qq_scale: float = 0.0  # 1-4 Coulomb scale (GROMACS fudgeQQ); 0 = LJ only


@dataclass
class BondTypeInfo:
    type_id: int
    style: str  # "harmonic", "backmap/harmonic", "backmap/table"
    keyword: str  # "at" or "cg" for backmap/* styles, "" for static
    params: list[float]
    table_file: str | None = None
    table_keyword: str | None = None


@dataclass
class AngleTypeInfo:
    type_id: int
    style: str
    keyword: str
    params: list[float]
    table_file: str | None = None
    table_keyword: str | None = None


@dataclass
class DihedralTypeInfo:
    type_id: int
    style: str
    keyword: str
    params: list[float]
    table_file: str | None = None
    table_keyword: str | None = None


@dataclass
class ImproperTypeInfo:
    type_id: int
    style: str  # "backmap/harmonic"
    keyword: str
    params: list[float]  # K (energy/rad^2, E = K dchi^2), chi0 (deg)


@dataclass
class AtomTypeInfo:
    type_id: int
    name: str
    mass: float
    is_cg: bool
    sigma: float = 0.0  # Angstrom
    epsilon: float = 0.0  # kcal/mol


@dataclass
class PairTypeInfo:
    """Which pair kind for each type pair: 'atomistic', 'cg', 'none'."""

    itype: int
    jtype: int
    kind: str
    sigma: float = 0.0
    epsilon: float = 0.0
    table_file: str | None = None
    table_keyword: str | None = None


@dataclass
class System:
    """Complete LAMMPS system representation."""

    atoms: list[LammpsAtom] = field(default_factory=list)
    bonds: list[LammpsBond] = field(default_factory=list)
    angles: list[LammpsAngle] = field(default_factory=list)
    dihedrals: list[LammpsDihedral] = field(default_factory=list)
    impropers: list[LammpsImproper] = field(default_factory=list)
    atom_types: list[AtomTypeInfo] = field(default_factory=list)
    bond_types: list[BondTypeInfo] = field(default_factory=list)
    angle_types: list[AngleTypeInfo] = field(default_factory=list)
    dihedral_types: list[DihedralTypeInfo] = field(default_factory=list)
    improper_types: list[ImproperTypeInfo] = field(default_factory=list)
    pair_types: list[PairTypeInfo] = field(default_factory=list)
    box: tuple[float, float, float] = (0.0, 0.0, 0.0)  # Angstrom
    cg_type_id: int = 0
    has_cross_bonds: bool = False
    has_cross_angles: bool = False
    has_cross_dihedrals: bool = False
    has_cross_pairs: bool = False
    cross_pairs: list[LammpsCrossPair] = field(default_factory=list)
    cross_pairs_file: str = "pairs.dat"
    fudge_lj: float = 1.0  # 1-4 LJ scale of the AT force field (GROMACS fudgeLJ)
    fudge_qq: float = 1.0  # 1-4 Coulomb scale (GROMACS fudgeQQ)
    write_image_flags: bool = False  # data file carries image flags

    # Atoms-per-bead by CG type ID (fix backmap apb); empty when every bead
    # and its fragment form one molecule, as the hybrid engine writes them.
    apb_by_cg_type: dict[int, int] = field(default_factory=dict)

    # Table files to convert: (src, dst) pairs
    table_files: list[tuple[str, str]] = field(default_factory=list)  # bond tables
    angle_table_files: list[tuple[str, str]] = field(default_factory=list)  # angle tables
    dihedral_table_files: list[tuple[str, str]] = field(default_factory=list)
    pair_table_files: list[tuple[str, str]] = field(default_factory=list)  # pair tables
