# pair_style backmap: special-bond factors and CG-specific exclusions

## Why

`pair_style backmap` called its sub-styles with factor 1 for every pair in
the neighbor list and ignored the special bits. That is correct only while
every `special_bonds` factor is 0 (pairs dropped from the list), which is what
all generated inputs use. The MARTINI 3 POPC example (research decision
`2026-09-25-martini-popc-example.md`) needs the CG model to exclude only bonded
neighbours (MARTINI nrexcl = 1) while the Slipids AT force field excludes up to
1-4 (nrexcl = 3). bakery supported this (`exclusion_at` / `exclusion_cg`); the
LAMMPS port had one global `special_bonds`.

## What Changes

- Each pair is evaluated with `force->special_lj/coul[sbmask(j)]`; factors
  below 1e-30 count as excluded.
- New trailing keyword `cg_special w12 w13 w14`: weights for CG-CG pairs that
  replace the `special_bonds` factors. Combined with tiny nonzero 1-3/1-4
  `special_bonds` factors this gives nrexcl 1 for CG and 3 for AT.
- Restart stores the new settings.

## Impact

- Code: `src/pair_backmap.{h,cpp}`; test `python/tests/test_lammps_pair_special.py`.
- No change for `special_bonds lj 0 0 0 coul 0 0 0` (every existing example).
