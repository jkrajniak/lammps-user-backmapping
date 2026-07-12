# Dihedral styles (USER-BACKMAP)

## `ryckaert` (static intra-bead)

For dihedrals whose four atoms belong to the same CG bead:

```
dihedral_style ryckaert
dihedral_coeff N C0 C1 C2 C3 C4 C5
```

Coefficients are in kcal/mol with LAMMPS polymer φ convention (trans = 180°).
GROMACS func-3 RB entries convert via `C_lammps[n] = (-1)^n × energy(C_gromacs[n])`.

## `backmap/ryckaert` (cross-bead AT)

Weight uses the single global lambda value and whether all four dihedral
atoms (**i, j, k, l**) map to the same CG bead (see
[theory: Force Weighting](../theory.md#force-weighting)): full strength once
&lambda;<sub>global</sub> > 0 if all four are in the same bead, otherwise
&lambda;<sub>global</sub> (linear fade-in).

```
dihedral_style backmap/ryckaert
dihedral_coeff N at C0 C1 C2 C3 C4 C5
```

Use `cg` instead of `at` for CG-only cross dihedrals.

## `backmap/table` (tabulated CG)

```
dihedral_style backmap/table linear 1000
dihedral_coeff M cg table_d1.table ENTRY
```

Tables use φ in degrees (-180..180), energy in kcal/mol, force in kcal/(mol·deg).
