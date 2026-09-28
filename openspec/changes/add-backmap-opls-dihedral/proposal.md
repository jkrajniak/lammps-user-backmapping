# Add the OPLS dihedral form (`dihedral_style backmap/opls`, native `opls` fragments)

## Why

Referee 1, comment 7 (COMPHY-D-26-00942) asks for the OPLS dihedral form
(`dihedral_style opls` in LAMMPS). The response plan decided (D2) to add a
concrete `backmap/opls` style. The OPLS form is an exact four-term cosine
series, so it can already be written as `backmap/ryckaert` or
`backmap/fourier`, but OPLS force fields for LAMMPS are commonly distributed
with `dihedral_style opls` coefficients (K1..K4), and asking users to convert
them by hand invites sign and factor-of-two errors.

## What Changes

- New `dihedral_style backmap/opls`: `dihedral_coeff N at|cg K1 K2 K3 K4`,
  E = w [K1/2 (1 + cos phi) + K2/2 (1 - cos 2phi) + K3/2 (1 + cos 3phi)
  + K4/2 (1 - cos 4phi)] (stock `dihedral_style opls` form, same phi
  convention, trans = 180 deg). Implemented as a subclass of
  `backmap/ryckaert` that converts the coefficients once in `coeff()`:
  C0 = K1/2 + K2 + K3/2, C1 = K1/2 - 3K3/2, C2 = -K2 + 4K4, C3 = 2K3,
  C4 = -4K4, C5 = 0. The compute kernel is the validated RB kernel, unchanged.
  `write_data` and restarts keep K1..K4.
- backmap-prep native AT fragments (`source.format: lammps`): the bounded
  script parser accepts `dihedral_style opls` and converts each
  `dihedral_coeff` to the six Ryckaert coefficients with the same formulas,
  so the rest of the pipeline is unchanged.

## Impact

- Specs: `bonded-backmap-styles` (ADDED), `backmap-input-generator` (MODIFIED:
  native fragment dihedral styles).
- Code: `src/dihedral_backmap_opls.{h,cpp}`,
  `python/src/backmap_prep/parsers/lammps_script_parser.py`.
- Tests: stock-equivalence test (`test_lammps_bonded_charmm.py` pattern),
  finite-difference entry in `test_lammps_force_fd.py`, parser unit tests.
- Additive: no existing style or generated input changes, so results produced
  with v1.3.0-rc5 are unaffected.
