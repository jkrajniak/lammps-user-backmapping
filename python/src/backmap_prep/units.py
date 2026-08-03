"""Unit conversion: GROMACS → LAMMPS real units."""

from __future__ import annotations

# GROMACS → LAMMPS real conversion factors
NM_TO_ANGSTROM = 10.0
KJ_TO_KCAL = 0.239006

# Derived factors
DISTANCE = NM_TO_ANGSTROM  # nm → Å
ENERGY = KJ_TO_KCAL  # kJ/mol → kcal/mol
FORCE = KJ_TO_KCAL / NM_TO_ANGSTROM  # kJ/(mol·nm) → kcal/(mol·Å)
TIME = 1000.0  # ps → fs
CHARGE = 1.0  # e → e
MASS = 1.0  # g/mol → g/mol
# GROMACS harmonic bond/angle potentials use E = (k/2)x²; LAMMPS's `harmonic`
# and `backmap/harmonic` bond/angle styles use E = Kx² (force factor -2K), so
# k must be halved on top of the unit conversion. See
# research/decisions/2026-07-21-at-harmonic-k-half-for-epplus-parity.md.
SPRING_BOND = KJ_TO_KCAL / (NM_TO_ANGSTROM**2) / 2.0  # kJ/(mol·nm²) → kcal/(mol·Å²), k→K
SPRING_ANGLE = KJ_TO_KCAL / 2.0  # kJ/(mol·rad²) → kcal/(mol·rad²), k→K
PRESSURE = 0.986923  # bar → atm


def distance(val: float) -> float:
    """nm → Angstrom"""
    return val * DISTANCE


def energy(val: float) -> float:
    """kJ/mol → kcal/mol"""
    return val * ENERGY


def force(val: float) -> float:
    """kJ/(mol·nm) → kcal/(mol·Å)"""
    return val * FORCE


def time(val: float) -> float:
    """ps → fs"""
    return val * TIME


def spring_bond(val: float) -> float:
    """GROMACS bond k (kJ/(mol·nm²), E=(k/2)x²) → LAMMPS K (kcal/(mol·Å²), E=Kx²)"""
    return val * SPRING_BOND


def spring_angle(val: float) -> float:
    """GROMACS angle k (kJ/(mol·rad²), E=(k/2)x²) → LAMMPS K (kcal/(mol·rad²), E=Kx²)"""
    return val * SPRING_ANGLE


def lj_pair_params(
    eps_i: float,
    eps_j: float,
    sig_i: float,
    sig_j: float,
    *,
    combination_rule: int = 2,
    fudge_lj: float = 1.0,
) -> tuple[float, float]:
    """GROMACS func-1 pair lookup → LAMMPS real sigma (Å), epsilon (kcal/mol)."""
    epsilon = fudge_lj * (eps_i * eps_j) ** 0.5 if eps_i > 0 and eps_j > 0 else 0.0
    if combination_rule == 2:
        sigma = 0.5 * (sig_i + sig_j)
    else:
        sigma = (sig_i * sig_j) ** 0.5 if sig_i > 0 and sig_j > 0 else 0.0
    return distance(sigma), energy(epsilon)


def gromacs_rb_to_lammps(coeffs: list[float]) -> list[float]:
    """Convert GROMACS func-3 RB coefficients to LAMMPS ryckaert C0..C5."""
    return [((-1.0) ** index) * energy(value) for index, value in enumerate(coeffs[:6])]


def gromacs_charmm_dihedral(k_kj: float, multiplicity: int, phi0_deg: float) -> list[float]:
    """Convert GROMACS func-1 dihedral to LAMMPS charmm K, n, delta (degrees)."""
    return [energy(k_kj), float(multiplicity), phi0_deg]


def gromacs_harmonic_dihedral(k_kj: float, multiplicity: int, phi0_deg: float) -> list[float]:
    """Convert GROMACS func-1 improper to LAMMPS harmonic dihedral K, sign, n."""
    sign = -1.0 if abs(abs(phi0_deg) - 180.0) < 1.0 else 1.0
    return [energy(k_kj), sign, float(multiplicity)]


def sigma(val: float) -> float:
    """LJ sigma: nm → Angstrom"""
    return val * DISTANCE


def epsilon(val: float) -> float:
    """LJ epsilon: kJ/mol → kcal/mol"""
    return val * ENERGY


def pressure(val: float) -> float:
    """bar → atm"""
    return val * PRESSURE
