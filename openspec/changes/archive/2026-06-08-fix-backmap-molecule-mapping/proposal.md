## Why

The `fix backmap` CG-to-AT molecule mapping is fundamentally broken for molecules with multiple CG beads. It assumes one CG atom per molecule (keyed by molecule ID), but real systems like dodecane have 6 CG beads per molecule (A-B-B-B-B-A). This causes incorrect force distribution, wrong COM tracking, mass mismatch warnings, and simulation crashes ("Bond atoms missing") within the first few thousand steps of the backmapping phase.

## What Changes

- **Redesign the CG→AT mapping** in `fix_backmap.cpp` to support multiple CG beads per molecule, where each CG bead maps to its own subset of AT atoms (not one CG per molecule).
- **Support multiple CG atom types** — change `cg_type` from a single integer to a set of CG types (e.g., `cg_types 1 2`), so both A and B beads are recognized as CG. **BREAKING**: The `cg_type T` keyword syntax changes to `cg_types T1 T2 ...` (or a compatible alternative).
- **Fix force distribution** in `post_force()` to distribute each CG bead's forces only to its own mapped AT atoms.
- **Fix COM tracking** in `end_of_step()` to update each CG bead's position to the COM of only its own mapped AT atoms.
- **Fix ghost atom mass accounting** in `build_molecule_map()` to avoid double-counting masses when molecules span MPI domains.
- **Update the Python writer** (`writers.py`) to emit the corrected fix syntax with multiple CG types.
- **Fix the dodecane data file** to ensure AT atom coordinates are wrapped inside the simulation box, resolving "Inconsistent image flags" warnings.

## Capabilities

### New Capabilities

(none — this is a bug fix to existing capabilities)

### Modified Capabilities

- `fix-backmap-resolution`: The CG-AT molecule mapping requirement changes from "one CG atom per molecule identified by `cg_type`" to "multiple CG beads per molecule, each mapping to a subset of AT atoms, identified by multiple CG types via `cg_types`". The force distribution and COM tracking requirements change to operate per-CG-bead rather than per-molecule.
- `backmap-input-generator`: The Python writer must emit the updated `cg_types` keyword (plural) with all CG type IDs, replacing the single `cg_type` value.

## Impact

- **C++ source**: `fix_backmap.cpp`, `fix_backmap.h` — significant refactor of `build_molecule_map()`, `post_force()`, `end_of_step()`, constructor, and `MolMap` data structure.
- **Python source**: `writers.py` — update the fix command generation to use `cg_types` with multiple type IDs.
- **Data file**: `examples/dodecane/dodecane.data` — wrap out-of-box AT atom coordinates.
- **Input script**: `examples/dodecane/in.dodecane` — use updated fix syntax.
- **Documentation**: Fix syntax docs, settings reference, and dodecane tutorial must reflect the new `cg_types` keyword.
- **Tests**: Python writer tests must be updated for the new keyword.
- **Breaking change**: Existing input scripts using `cg_type T` will need to be updated to `cg_types T` (or `cg_types T1 T2 ...`).
