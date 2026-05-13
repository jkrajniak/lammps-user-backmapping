## 1. C++ fix_backmap — data structures and parsing

- [x] 1.1 Replace `int cg_type` with `std::set<int> cg_types` in `fix_backmap.h`; update constructor to parse multiple type IDs after `cg_type` keyword (read until next recognized keyword or non-integer)
- [x] 1.2 Replace `struct MolMap` / `std::map<tagint, MolMap> mol_map` with `struct BeadMap` / `std::vector<BeadMap> bead_map` in `fix_backmap.h`
- [x] 1.3 Update `initial_integrate()` to zero velocities/forces for all atoms matching any type in `cg_types` (not just one type)

## 2. C++ fix_backmap — build_molecule_map redesign

- [x] 2.1 Rewrite `build_molecule_map()`: collect CG and AT atoms per molecule, sort by tag, validate `num_AT % num_CG == 0`, assign AT atoms to CG beads sequentially
- [x] 2.2 Fix ghost atom handling: use `atom->mass[atom->type[i]]` (per-type mass) for all atoms including ghosts; do not double-count local+ghost for the same physical atom
- [x] 2.3 Rewrite `validate_masses()` to check per-bead mass consistency (CG bead type mass vs sum of its mapped AT atom type masses)

## 3. C++ fix_backmap — force distribution and COM tracking

- [x] 3.1 Update `post_force()` to iterate over `bead_map` and distribute each CG bead's forces only to its own mapped AT atoms
- [x] 3.2 Update `end_of_step()` to compute COM per CG bead using only its mapped AT atoms; update lambda increment to iterate over all local atoms

## 4. C++ fix_backmap — diagnostics and init

- [x] 4.1 Update `init()` log message to report number of bead mappings (not molecule mappings) and list all CG types
- [x] 4.2 Add error for molecules where `num_AT % num_CG != 0`

## 5. Python writer — fix command and data file

- [x] 5.1 Update `writers.py` to emit `cg_type T1 T2 ...` with all CG type IDs (sorted ascending) from `system.atom_types` where `is_cg=True`; remove `system.cg_type_id` single-value reference
- [x] 5.2 Add coordinate wrapping in `write_lammps_data()`: wrap atom positions into `[0, L)` before writing the Atoms section
- [x] 5.3 Update Python writer tests to verify multi-CG-type fix command and wrapped coordinates

## 6. Example and validation

- [x] 6.1 Regenerate `examples/dodecane/dodecane.data` with wrapped coordinates using the updated generator
- [x] 6.2 Update `examples/dodecane/in.dodecane` to use `cg_type 1 2` syntax
- [ ] 6.3 Run the dodecane example end-to-end (all 3 phases) and verify no crashes, no mass mismatch warnings, no bond-missing errors

## 7. Documentation and changelog

- [x] 7.1 Update fix backmap documentation (`docs/`) with new `cg_type` multi-type syntax and atom ordering convention
- [x] 7.2 Update CHANGELOG.md with the fix under `[Unreleased]` → `Fixed`
- [x] 7.3 Update README.md if it references the fix syntax
