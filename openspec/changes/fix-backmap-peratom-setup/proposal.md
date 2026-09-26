# fix backmap: per-atom COM/force output and CG force distribution at setup

## Why

Review comment R2-3 asks for a static decomposition-parity test: identical
coordinates, `run 0`, compare bead COMs, redistributed CG forces and total
forces across rank counts and processor grids. `fix backmap` exposed only a
per-atom lambda vector (response plan item A1).

While adding it: `FixBackmap::setup()` never called `post_force()`. LAMMPS
evaluates forces in setup at the start of every run and uses them for the
first velocity half-kick; the CG forces then stayed on the beads and the AT
atoms got none.

## What Changes

- New keyword `peratom lambda|full`; `full` gives a 7-column per-atom array
  (lambda, bead position, CG force share / bead CG force), filled in
  `post_force()`. Default (`lambda`) unchanged.
- `setup()` calls `post_force()` after building the bead map.

## Impact

- Code: `src/fix_backmap.{h,cpp}`; test `python/tests/test_lammps_backmap_setup.py`.
- Results change only in the first step of each run with active CG forces
  (start of the lambda ramp; frozen-CG and lambda = 1 stages unaffected).
