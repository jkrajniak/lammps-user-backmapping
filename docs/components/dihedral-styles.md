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
[theory: Force Weighting](../theory.md#force-weighting)): always full
strength if all four are in the same bead, otherwise
&lambda;<sub>global</sub> (linear fade-in).

```
dihedral_style backmap/ryckaert
dihedral_coeff N at C0 C1 C2 C3 C4 C5
```

Use `cg` instead of `at` for CG-only cross dihedrals.

## `backmap/fourier` (multi-term periodic, arbitrary phase)

```
dihedral_style backmap/fourier
dihedral_coeff N at|cg m K1 n1 d1 ... Km nm dm
```

\[
E = w \sum_{i=1}^{m} K_i \left[ 1 + \cos(n_i \phi - d_i) \right]
\]

Same form and kernel as LAMMPS `dihedral_style fourier`, with the phase
\( d_i \) in degrees (any value, e.g. 60°) and multiplicities from 0 up. The
weight is the same as for `backmap/ryckaert`. GROMACS dihedral functions 1, 4
and 9 (\( \phi_s, k, n \), one line per term) map term by term to
\( d = \phi_s \), \( K = k \), \( n \). Use it for CHARMM/AMBER-type
dihedrals that do not reduce to Ryckaert-Bellemans (multiplicity 6, phases
other than 0°/180°).

## `improper_style backmap/harmonic`

```
improper_style backmap/harmonic
improper_coeff N at|cg K chi0
```

\[
E = w\, K (\chi - \chi_0)^2
\]

\( \chi \) is the angle between the planes (I,J,K) and (J,K,L), as in LAMMPS
`improper_style harmonic`; \( \chi_0 \) in degrees. GROMACS improper function 2
(\( \xi_0, k_\xi \), \( \tfrac12 k \) convention) maps to \( K = k_\xi/2 \),
\( \chi_0 = \xi_0 \).

## `backmap/table` (tabulated CG)

```
dihedral_style backmap/table linear 1000
dihedral_coeff M cg table_d1.table ENTRY
```

Tables use φ in degrees (-180..180), energy in kcal/mol, force in kcal/(mol·deg).

## Data files

`write_data` writes the coefficients of `ryckaert`, `backmap/ryckaert` and
`backmap/harmonic` in `dihedral_coeff` argument order, so a written data file
can be read back directly. `backmap/table` coefficients are not written.
