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
SPRING_BOND = KJ_TO_KCAL / (NM_TO_ANGSTROM**2)  # kJ/(mol·nm²) → kcal/(mol·Å²)
SPRING_ANGLE = KJ_TO_KCAL  # kJ/(mol·rad²) → kcal/(mol·rad²)


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
    """kJ/(mol·nm²) → kcal/(mol·Å²)"""
    return val * SPRING_BOND


def spring_angle(val: float) -> float:
    """kJ/(mol·rad²) → kcal/(mol·rad²)"""
    return val * SPRING_ANGLE


def sigma(val: float) -> float:
    """LJ sigma: nm → Angstrom"""
    return val * DISTANCE


def epsilon(val: float) -> float:
    """LJ epsilon: kJ/mol → kcal/mol"""
    return val * ENERGY
