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

## `backmap/opls` (OPLS form)

```
dihedral_style backmap/opls
dihedral_coeff N at|cg K1 K2 K3 K4
```

\[
E = w \left[ \tfrac{K_1}{2}(1 + \cos\phi) + \tfrac{K_2}{2}(1 - \cos 2\phi)
+ \tfrac{K_3}{2}(1 + \cos 3\phi) + \tfrac{K_4}{2}(1 - \cos 4\phi) \right]
\]

Same form, coefficients and \( \phi \) convention (trans = 180°) as LAMMPS
`dihedral_style opls`, so OPLS coefficients from a LAMMPS force-field file can
be used unchanged. The series is exact in \( \cos\phi \); the style converts
K1..K4 once to the `backmap/ryckaert` coefficients

\[
C_0 = \tfrac{K_1}{2} + K_2 + \tfrac{K_3}{2},\;
C_1 = \tfrac{K_1}{2} - \tfrac{3K_3}{2},\;
C_2 = -K_2 + 4K_4,\;
C_3 = 2K_3,\;
C_4 = -4K_4,\;
C_5 = 0
\]

and uses the Ryckaert-Bellemans kernel. The weight is the same as for
`backmap/ryckaert`. `write_data` and restart files keep K1..K4.

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
