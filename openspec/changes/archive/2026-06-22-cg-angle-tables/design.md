# Design — CG angle tables

## Context

Bond tabulated CG cross-interactions already work via `bond_style backmap/table`.
Angle func-8 entries in hybrid TOP use the same GROMACS convention: tablenr in
params[0], file `table_a{tablenr}.xvg` from `prep.tables_dir`.

Espresso++ reads angles with `TabulatedAngular` (theta in radians internally after
`convertTable` from degrees in xvg).

## Goals

- All 200 rim135 func-8 cross-angles use real `backmap/table cg` coefficients.
- At λ=0 (pure CG), CG table angles exert full tabulated torque (`w = 1 − 0 = 1`).
- Step-0 LAMMPS parse succeeds with rebuilt USER-BACKMAP package.

## Decisions

### 1. Separate C++ style (not bond reuse)

`angle_style backmap/table` mirrors `bond_backmap_table` for I/O and interpolation
but applies forces using the i–j–k geometry from `angle_backmap_harmonic`.

Lambda weight uses endpoints i and k (not j).

### 2. Tablenr resolution

`int(angle.params[0])` → `table_a{index}.xvg` only (never `table_b*`).

### 3. Angle XVG unit conversion

GROMACS `table_a*.xvg` column 0 is **degrees** (0–180). Current `_convert_xvg`
incorrectly applies `units.distance` (nm→Å).

New `_convert_angle_xvg`:
- col0: degrees (pass-through to LAMMPS table, which uses degrees for angles)
- col1: kJ/mol → kcal/mol
- col2: `-dV/dθ` in kJ/(mol·rad) → kcal/(mol·deg) via `force * (π/180)`

This matches bakery `convertTable` angle branch (`fd * 180/π` when converting
derivative w.r.t. radians to degrees).

### 4. Hybrid writer

When both harmonic and table backmap angles exist:

```
angle_style hybrid backmap/harmonic backmap/table linear 1000
angle_coeff N backmap/table cg table_a1.table ENTRY
```

### 5. Validation gate

Before Tier B re-run:
- No `angle_coeff … cg 0.000000 0.0000` in generated `in.rim135`
- `table_a1.table` and `table_a2.table` present
- Unit tests green

## Follow-up (out of scope)

Intra-bead AT angles exported as `backmap/harmonic` instead of static `harmonic`
may affect weighting; track separately from tabulated CG angles.
