# Bond Styles

The backmapping package provides two lambda-weighted bond styles for cross-CG
bonded interactions.

## backmap/harmonic

### Syntax

```
bond_style backmap/harmonic
bond_coeff N at/cg K r0
```

- **N** -- bond type number
- **at** or **cg** -- weighting mode
    - `at`: weight = &lambda;<sub>i</sub> &times; &lambda;<sub>j</sub>
      (fades in during backmapping)
    - `cg`: weight = 1 - &lambda;<sub>i</sub> &times; &lambda;<sub>j</sub>
      (fades out during backmapping)
- **K** -- force constant (energy/distance&sup2; units)
- **r0** -- equilibrium bond length (distance units)

### Description

A harmonic bond potential scaled by the lambda weight:

\[
E = w \times \frac{1}{2} K (r - r_0)^2
\]

\[
F = -w \times K (r - r_0)
\]

where \( w \) is computed from the lambda values of the two bonded atoms.

This style is used for AT cross bonds (harmonic bonds between atoms in
different CG beads) and, less commonly, for CG cross bonds when a harmonic
approximation is acceptable.

### Example

```
bond_style hybrid harmonic backmap/harmonic

# Intra-bead AT bond (standard, no lambda weighting)
bond_coeff 1 harmonic 800.0 1.53

# AT cross bond (fades in with lambda)
bond_coeff 2 backmap/harmonic at 800.0 1.53
```

---

## backmap/table {: #backmap-table }

### Syntax

```
bond_style backmap/table style N
bond_coeff M at/cg filename keyword
```

**Style command:**

- **style** -- interpolation style: `linear` or `spline`
- **N** -- number of table entries for interpolation

**Coeff command:**

- **M** -- bond type number
- **at** or **cg** -- weighting mode (same as backmap/harmonic)
- **filename** -- path to the table file
- **keyword** -- section name within the table file

### Description

A tabulated bond potential scaled by the lambda weight:

\[
E = w \times E_\text{table}(r)
\]

\[
F = w \times F_\text{table}(r)
\]

This style is primarily used for CG cross bonds where the interaction
potential was derived from coarse-graining methods (e.g., iterative
Boltzmann inversion) and is provided as a numerical table.

### Table File Format

The table file contains sections identified by a keyword. Each section has
a header line with the keyword and number of entries, followed by data lines:

```
ENTRY 500
1  0.200  100.5  -502.5
2  0.204  95.3   -476.5
3  0.208  90.2   -451.0
...
```

Columns: index, distance, energy, force.

### Example

```
bond_style backmap/table linear 1000

# CG cross bond from tabulated potential (fades out with lambda)
bond_coeff 1 cg table_b1.table ENTRY
```

### Interpolation

- **linear**: piecewise linear interpolation between table entries
- **spline**: natural cubic spline interpolation (smoother but slower)

The table is re-sampled to `N` evenly-spaced points during initialization.

!!! warning
    Bond distances outside the table range \[r_min, r_max\] will cause a
    fatal error. Make sure your table covers the full range of expected
    distances.

---

## Using with hybrid bond_style

In practice, backmapping simulations use a hybrid bond style combining
standard and lambda-weighted bonds:

```
bond_style hybrid harmonic backmap/harmonic backmap/table linear 1000

# Intra-bead AT bonds (always active)
bond_coeff 1 harmonic 800.0 1.53

# CG cross bonds (tabulated, fade out)
bond_coeff 2 backmap/table cg table_b1.table ENTRY

# AT cross bonds (harmonic, fade in)
bond_coeff 3 backmap/harmonic at 800.0 1.53
```

## Requirements

Both bond styles require [`fix backmap`](fix-backmap.md) to be defined.
They read per-atom lambda values from the fix at each timestep.

## Restart

Bond coefficients (K, r0, is_cg) are written to and read from restart files.

## Related

- [fix backmap](fix-backmap.md) -- provides per-atom lambda values
- [angle_style backmap/harmonic](angle-styles.md) -- lambda-weighted angles
