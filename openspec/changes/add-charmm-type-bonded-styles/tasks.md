# Tasks

## 1. C++ styles
- [x] 1.1 `angle_style backmap/charmm` (compute, coeff, single, restart, write_data)
- [x] 1.2 `dihedral_style backmap/fourier` (compute, coeff, restart, write_data)
- [x] 1.3 `improper_style backmap/harmonic` (compute, coeff, restart, write_data)
- [x] 1.4 Regression tests against the stock styles: at lambda = 1 (at) and
      lambda = 0 (cg) the energy and forces equal the stock style's; at lambda = 0.5
      they are half; intra-bead terms are unweighted.
- [x] 1.5 Docs: `docs/components/` pages, `src/README`, CHANGELOG.

## 2. backmap-prep (follow-up PR, stacked on the builder chain)
- [ ] 2.1 Angle func 5 -> backmap/charmm; dihedral func 1/4/9 (multi-line
      [ dihedraltypes ], X wildcards) -> backmap/fourier; improper func 2 ->
      backmap/harmonic; AT-only file uses stock charmm / fourier / harmonic.
- [ ] 2.2 GROMACS energy check (`tests/reference/gmx_energy_check.py`) on a
      Slipids POPC system: every term agrees.
