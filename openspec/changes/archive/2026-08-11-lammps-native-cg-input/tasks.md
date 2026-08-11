## 1. Schema

- [x] 1.1 Extend `CGSystem` in `schema.py`: `format: Literal["gromacs", "lammps"]`, make `coordinates`/`topology` optional, add an optional `data` field. Validate `gromacs` requires `coordinates`+`topology` and `lammps` requires `data` (matching the existing `check_engine_requirements` validation style).
- [x] 1.2 Add a schema-level error message for an unrecognized `cg_system.format` value matching the spec's "Unsupported format" scenario wording.

## 2. LAMMPS data-file reader

- [x] 2.1 Create `parsers/lammps_data_parser.py`: parse box bounds (`xlo xhi` etc.), `Masses`, and `Atoms # full` sections. Tolerate and skip any other section (`Bonds`, `Angles`, `Velocities`, ...) without error — `builder.py` never reads CG bonded topology from `cg_top` (confirmed by exhaustive grep; see `design.md` correction note), so these sections are inert either way.
- [x] 2.2 Return objects that satisfy the exact shape `builder.py` consumes from `gro_parser.GroFile`/`top_parser.Topology`: reuse those dataclasses directly rather than inventing parallel ones — a `GroFile` (title, atoms, box) and a `Topology` with one `MoleculeType` whose `atoms` list is ordered exactly as the `data` file's `Atoms` section. Use each atom's numeric LAMMPS type ID, as a string, for `TopAtom.type`/`AtomType` keys (no atom-name field exists in a `data` file; `builder.py` never needs one for the CG side).
- [x] 2.3 Document the atom-ordering assumption (contiguous per-molecule blocks, matching `cg_start = mol_idx * cg_atom_count` in `builder.py:690`) directly in the parser's module docstring, since there's no ordering metadata in a `data` file to validate against automatically.
- [x] 2.4 Unit tests: minimal single-molecule and multi-molecule `data` file fixtures; a missing-`Masses`-section case and a missing-`Atoms`-section case, each producing a clear, named error.

## 3. Builder integration

- [x] 3.1 In `build_system()`, branch construction of `cg_gro`/`cg_top`/`cg_mol` on `settings.cg_system.format`: `parse_gro`/`parse_top` for `gromacs`, the new LAMMPS data parser for `lammps`.
- [x] 3.2 At the three CG-specific `units.distance()` call sites (`builder.py:207-209` box, `builder.py:700-702` and `builder.py:728-730` CG positions), use a local identity-vs-real closure keyed on `cg_system.format == "lammps"` instead of always calling `units.distance()`. Do not touch any AT-side conversion call site.
- [x] 3.3 Confirm zero changes needed elsewhere in `builder.py` — the CG positional/type-keyed atom-matching logic, `cross_interactions` processing, and `table_groups` processing are unaffected by `cg_system.format`. This is the acceptance bar for "scoped to CG side only."
- [x] 3.4 Confirm AT-fragment (`molecules[].source`) parsing is untouched and still GROMACS-only.

## 4. `.table` pair-lookup passthrough (format-independent fix)

- [x] 4.1 In the `table_groups` CG pair-table lookup (`builder.py:656-676`), check for `table_<a>_<b>.table` (and the reversed name) before falling back to the existing `table_<a>_<b>.xvg` construction. When a `.table` file is found, set `pt.table_file`/`pt.table_keyword` directly without adding it to `sys.pair_table_files` (no conversion needed — nothing to run through `table_converter.py`).
- [x] 4.2 Unit test: given both `table_A_A.xvg` and `table_A_A.table` present, the `.table` file is used and no `.xvg` conversion is triggered.

## 5. Verification example (`examples/dodecane-lammps-cg/`)

`examples/dodecane/prepare_cg.py` already extracts a standalone pure-CG `data` file (box, `Masses`, `Atoms # full`, plus an unused-by-us `Bonds` section) from the hybrid output — today only for a standalone post-hoc re-equilibration check (`in.dodecane_cg_equil`), one direction away from what this change needs. This task group closes the loop: feed that extracted CG system back into `backmap-prep` as input.

- [x] 5.1 Add `examples/dodecane-lammps-cg/`, seeded from `examples/dodecane/`: reuse `dodecane_single.gro`/`topol_aa.top` (AT fragment, stays GROMACS — out of scope for this change), the existing `cross_interactions`/`simulation` blocks from `settings.yaml` unchanged, and the existing `table_A_A.table`/`table_A_B.table`/`table_B_B.table`/`table_b1.table` files (exercises the task 4 `.table`-first lookup).
- [x] 5.2 Derive `dodecane_cg.data` the same way `prepare_cg.py` does today, from `examples/dodecane/dodecane.data`. Confirm its `Atoms` section ordering (contiguous per-molecule blocks) matches what `builder.py:683-729` expects (task 2.3).
- [x] 5.3 Write `settings.yaml` for this example identical to `examples/dodecane/settings.yaml` except `cg_system: {format: lammps, data: dodecane_cg.data}` and `beads[].type` values changed from `A`/`B` to their numeric LAMMPS type-ID strings (`"1"`/`"2"`, from `dodecane_cg.data`'s `Masses` section) — and `simulation.table_groups` updated to match (`["1", "2"]`).
- [x] 5.4 Run `backmap-prep` on both `examples/dodecane/settings.yaml` (GROMACS CG path) and `examples/dodecane-lammps-cg/settings.yaml` (LAMMPS-native CG path); diff the resulting hybrid `.data`/`.in` output. Positions/box/coefficients must match to numerical tolerance.
- [x] 5.5 Add a short `README.md` to `examples/dodecane-lammps-cg/` explaining it's a parity example, pointing back at `examples/dodecane/README.md` for the physical system description, and noting the `type: "1"`/`"2"` convention.
- [x] 5.6 ~~Record the round-trip result as a research experiment note~~ — decided during implementation not to: `research/experiments/` is for paper-facing physics/simulation runs (SHA-256-pinned evidence), not code-level parity checks. The acceptance evidence lives instead as permanent, automated regression tests (`tests/test_examples.py::TestLammpsNativeCgParity`), which is stronger than a one-off note since it re-verifies on every test run.

## 6. Docs

- [x] 6.1 Add a "LAMMPS-native CG system" section to `docs/cli/backmap-prep.md`: the `data`-file requirements (`Masses` + `Atoms # full`, `units real` assumed), the numeric-type-ID bead convention, and a minimal worked example pointing at `examples/dodecane-lammps-cg/`.
- [x] 6.2 Update the settings-format reference to show the `cg_system.format: lammps` variant alongside the existing `gromacs` example.

## 7. Spec/paper follow-up

- [x] 7.1 Confirm `openspec archive` picks up the "CG system section", "Source file parsing", "Unit conversion", and "Feature phasing" requirement changes cleanly against `openspec/specs/backmap-input-generator/spec.md`.
- [x] 7.2 If this lands before the next paper revision, consider whether the manuscript's input-format description (`paper/main.tex:105-113`) should mention the LAMMPS-native CG option; otherwise note it as a documented future-work item consistent with the current text.
