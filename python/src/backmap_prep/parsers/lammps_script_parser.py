"""Bounded LAMMPS input-script parser for AT-fragment force-field coefficients.

Covers exactly the four coefficient families `builder.py` consumes for AT
fragments (``molecules[].source.format: lammps``): ``bond_coeff`` (only
``bond_style harmonic``), ``angle_coeff`` (only ``angle_style harmonic``),
``dihedral_coeff`` (only ``dihedral_style ryckaert`` — the package's own
native style, matching the form ``units.gromacs_rb_to_lammps`` already
converts GROMACS RB coefficients into, so a LAMMPS-native fragment supplies
these directly with no conversion), and ``pair_coeff i i`` diagonal (self)
entries (cross-type LJ parameters are always computed in Python via mixing,
never read from a source file). Nothing else this codebase would consume is
in scope: AT fragments are always analytic here, never tabulated, and
cross-type ``pair_coeff`` lines are tolerated but ignored rather than
rejected, since a real production input script (e.g. ``in.dodecane_at``) may
already contain them and should be reusable as-is.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pathlib import Path

_SUPPORTED_BOND_STYLES = frozenset({"harmonic"})
_SUPPORTED_ANGLE_STYLES = frozenset({"harmonic"})
_SUPPORTED_DIHEDRAL_STYLES = frozenset({"ryckaert"})
# Force-field style commands this reader understands and validates. Any other
# command ending in "_style" (e.g. atom_style, improper_style) is unrelated to
# the four coefficient families above -- ignored, not rejected, except
# improper_style, which is a force-field style this reader has no consumer for.
_FORCEFIELD_STYLE_COMMANDS = frozenset(
    {"bond_style", "angle_style", "dihedral_style", "pair_style", "improper_style"}
)


@dataclass
class AtFragmentCoefficients:
    """Coefficients extracted from a bounded AT-fragment input script.

    ``bond``/``angle`` map type ID to ``(k, r0)``/``(k, theta0)`` in LAMMPS's
    own coefficient order (matching ``bond_coeff``/``angle_coeff`` argument
    order). ``dihedral`` maps type ID to the six raw ``dihedral_coeff``
    values, unconverted. ``pair`` maps atom type ID to ``(epsilon, sigma)``
    from that type's diagonal ``pair_coeff i i`` entry.
    """

    bond: dict[int, tuple[float, float]] = field(default_factory=dict)
    angle: dict[int, tuple[float, float]] = field(default_factory=dict)
    dihedral: dict[int, list[float]] = field(default_factory=dict)
    pair: dict[int, tuple[float, float]] = field(default_factory=dict)


def parse_at_fragment_script(path: Path) -> AtFragmentCoefficients:
    """Parse a bounded LAMMPS input-script fragment for AT-fragment coefficients.

    Requires a literal ``units real`` declaration (no conversion is applied
    anywhere on this path). Any ``*_style`` directive outside the supported
    bond/angle/dihedral families aborts with a named error. All other
    commands (``run``, ``thermo``, ``fix``, ``pair_style``, ...) are ignored.
    """
    text = path.read_text()

    coeffs = AtFragmentCoefficients()
    units_seen: str | None = None

    for raw_line in text.splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line:
            continue
        tokens = line.split()
        cmd = tokens[0]

        if cmd == "units":
            units_seen = tokens[1] if len(tokens) > 1 else None
        elif cmd == "bond_style":
            style = tokens[1] if len(tokens) > 1 else None
            if style not in _SUPPORTED_BOND_STYLES:
                raise ValueError(
                    f"unsupported bond_style {style!r} in {path} "
                    f"(supported: {sorted(_SUPPORTED_BOND_STYLES)})"
                )
        elif cmd == "angle_style":
            style = tokens[1] if len(tokens) > 1 else None
            if style not in _SUPPORTED_ANGLE_STYLES:
                raise ValueError(
                    f"unsupported angle_style {style!r} in {path} "
                    f"(supported: {sorted(_SUPPORTED_ANGLE_STYLES)})"
                )
        elif cmd == "dihedral_style":
            style = tokens[1] if len(tokens) > 1 else None
            if style not in _SUPPORTED_DIHEDRAL_STYLES:
                raise ValueError(
                    f"unsupported dihedral_style {style!r} in {path} "
                    f"(supported: {sorted(_SUPPORTED_DIHEDRAL_STYLES)})"
                )
        elif cmd == "pair_style":
            pass  # any pair_style is accepted; only pair_coeff values matter
        elif cmd == "bond_coeff":
            if len(tokens) < 4:
                raise ValueError(f"malformed bond_coeff line in {path}: {raw_line!r}")
            coeffs.bond[int(tokens[1])] = (float(tokens[2]), float(tokens[3]))
        elif cmd == "angle_coeff":
            if len(tokens) < 4:
                raise ValueError(f"malformed angle_coeff line in {path}: {raw_line!r}")
            coeffs.angle[int(tokens[1])] = (float(tokens[2]), float(tokens[3]))
        elif cmd == "dihedral_coeff":
            if len(tokens) < 8:
                raise ValueError(f"malformed dihedral_coeff line in {path}: {raw_line!r}")
            coeffs.dihedral[int(tokens[1])] = [float(t) for t in tokens[2:8]]
        elif cmd == "pair_coeff":
            if len(tokens) < 5:
                raise ValueError(f"malformed pair_coeff line in {path}: {raw_line!r}")
            i, j = int(tokens[1]), int(tokens[2])
            if i != j:
                continue  # cross-type LJ is always Python-computed via mixing
            coeffs.pair[i] = (float(tokens[3]), float(tokens[4]))
        elif cmd in _FORCEFIELD_STYLE_COMMANDS:
            raise ValueError(f"unsupported style directive {line!r} in {path}")
        # everything else (run, thermo, fix, neighbor, special_bonds, ...) is ignored

    if units_seen != "real":
        raise ValueError(
            f"LAMMPS input script {path} must declare 'units real' (found: {units_seen!r})"
        )

    return coeffs
