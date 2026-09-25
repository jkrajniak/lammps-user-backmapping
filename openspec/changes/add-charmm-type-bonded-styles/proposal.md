# Add CHARMM-type bonded backmap styles (Urey-Bradley, Fourier, harmonic improper)

## Why

Two needs of the COMPHY-D-26-00942 revision meet here:

- **R1-7** asks for impropers and OPLS-type dihedrals.
- **MARTINI 3 POPC example** (research decision `2026-09-25-martini-popc-example.md`):
  the atomistic Slipids POPC (and equally CHARMM36) uses Urey-Bradley angles
  (GROMACS angle func 5; 223 of 256 POPC angles have a nonzero UB term),
  multi-term periodic dihedrals with arbitrary phase (func 9; multiplicities up
  to 6, phases 0/60/180, so not exactly Ryckaert-Bellemans) and harmonic
  impropers (func 2).

The package has only harmonic angles, Ryckaert-Bellemans / single-term harmonic
dihedrals (`backmap/harmonic` allows d = +-1 only) and no impropers, so these
molecules cannot be backmapped.

## What Changes

- New `angle_style backmap/charmm`: `angle_coeff N at|cg K theta0 K_ub r_ub`,
  E = w [K (theta - theta0)^2 + K_ub (r13 - r_ub)^2] (LAMMPS `angle charmm` form).
- New `dihedral_style backmap/fourier`: `dihedral_coeff N at|cg m K1 n1 d1 ...`,
  E = w sum_i K_i [1 + cos(n_i phi - d_i)] (LAMMPS `dihedral fourier` form,
  arbitrary phase d_i in degrees).
- New `improper_style backmap/harmonic`: `improper_coeff N at|cg K chi`,
  E = w K (chi - chi0)^2 (LAMMPS `improper harmonic` form).
- Weight w as for every `backmap/*` style (`compute_weight3`): 1 intra-bead,
  lambda for `at` terms across beads, 1 - lambda for `cg` terms.
- Each style supports `write_data`, restart, and forces are applied with the
  standard `newton_bond || i < nlocal` guard.
- A follow-up change in `backmap-prep` maps GROMACS angle func 5, dihedral
  func 1/4/9 and improper func 2 onto these styles (and onto the stock
  `charmm` / `fourier` / `harmonic` styles for the AT-only force field).

## Impact

- Specs: `bonded-backmap-styles` (ADDED requirements below).
- Code: `src/angle_backmap_charmm.*`, `src/dihedral_backmap_fourier.*`,
  `src/improper_backmap_harmonic.*`; tests in `python/tests/test_lammps_styles_charmm.py`.
- Kernels are the stock LAMMPS ones (`angle_charmm.cpp`, `dihedral_fourier.cpp`,
  `improper_harmonic.cpp`) with the weight applied to energy, force and virial;
  no new physics.
