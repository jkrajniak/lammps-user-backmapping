## 1. Schema

- [x] 1.1 Extend `SourceFiles` in `schema.py`: `format: Literal["gromacs", "lammps"]`, `data: str | None`, `input_script: str | None`. Validate `gromacs` requires `coordinates`+`topology`, `lammps` requires `data`+`input_script`.
- [x] 1.2 Reject `format: lammps` combined with a degree-dependent (`list[DegreeSourceFile]`) `coordinates`/`topology` value, and reject it when `prep.engine == "network"` — same validation layer `check_engine_requirements` already uses for `cg_system.format: lammps` + network engine.
- [x] 1.3 Error message wording matching the spec's "Unsupported format" and "LAMMPS-native source rejected for degree-dependent mapping" scenarios.

## 2. LAMMPS data-file connectivity reader (AT fragment)

- [x] 2.1 Extend `parsers/lammps_data_parser.py` with an AT-fragment variant: parse `Masses`, `Atoms # full` (as for CG), plus `Bonds`, `Angles`, `Dihedrals` sections this time (unlike the CG reader, these are load-bearing here — see `design.md`). Produce `TopBond`/`TopAngle`/`TopDihedral` objects referencing type IDs.
- [x] 2.2 Use each atom's numeric LAMMPS atom ID, as a string, for `TopAtom.name` (mirrors `parse_cg_system`'s use of the numeric type ID for `TopAtom.type`) — this is what `beads[].atoms` entries will reference.
- [x] 2.3 Unit tests: minimal AT-fragment `data` file fixture with bonds, angles, and dihedrals; missing-section error cases.

## 3. LAMMPS input-script coefficient reader (AT fragment)

- [x] 3.1 Create `parsers/lammps_script_parser.py` (the module the CG-side change scoped out) supporting exactly: `units` (must be `real`), `bond_style`/`bond_coeff` (`harmonic` only), `angle_style`/`angle_coeff` (`harmonic` only), `dihedral_style`/`dihedral_coeff` (`ryckaert` only), `pair_style`/`pair_coeff` (diagonal `i i` entries only; `i != j` lines ignored, not rejected).
- [x] 3.2 Ignore recognized-but-irrelevant commands (`pair_modify`, `special_bonds`, `neighbor`, `fix`, `run`, `thermo`, `compute`, output commands) rather than erroring on them — `in.dodecane_at` is a real production script with several of these; reusing it as-is (task 5) is the acceptance bar.
- [x] 3.3 Abort with a named "unsupported style" error for any `*_style` directive outside the supported set (e.g. `bond_style morse`).
- [x] 3.4 Abort if `units real` is not declared.
- [x] 3.5 Cross-check every type ID referenced by the `data` file's `Bonds`/`Angles`/`Dihedrals` sections has a corresponding coefficient entry; abort naming the missing type ID otherwise.
- [x] 3.6 Decide the `include`-support open question from `design.md` consistently with however the CG-side change resolved it.
- [x] 3.7 Unit tests: bond/angle/dihedral/pair coefficient extraction, missing-units error, unsupported-style error, missing-coefficient error, ignored cross-pair-coeff line.

## 4. Builder integration

- [x] 4.1 Branch AT-fragment loading (`at_gro`/`at_top`) in `build_system()` on `mol_def.source.format`, dispatching to the GROMACS or LAMMPS-native reader pair.
- [x] 4.2 Confirm zero changes needed to the bead-mapping/placement logic itself, and to `resolve_dihedral_params`/`gromacs_rb_to_lammps` call sites (the LAMMPS-native path supplies already-resolved RB coefficients directly, bypassing `resolve_dihedral_params`'s wildcard lookup entirely — verify this bypass is clean, not just untested).
- [x] 4.3 Confirm CG-system loading (either format) is untouched.
- [x] 4.4 No unit conversion applied to AT-fragment values read via the LAMMPS-native path (positions, bond/angle/dihedral coefficients, pair sigma/epsilon) — same principle as the CG-side change, extended to every AT-fragment call site that currently calls `units.distance`/`units.spring_bond`/`units.spring_angle`/`units.sigma`/`units.epsilon`/`units.gromacs_rb_to_lammps` on AT-sourced values. Enumerate these call sites explicitly (they are currently intermixed with CG-side ones in several of the same code blocks, unlike the CG change where CG and AT conversions were already cleanly separated) before writing the identity-vs-real closures.

## 5. Verification example (`examples/pe-lammps/`)

**Revised during implementation** — see `design.md`'s "IMPORTANT CORRECTION": `extract_at_frame.py`-based extraction (the original plan below) produces a *wrong* AT-fragment template, because AT placement uses a rigid single-molecule offset-from-atom-1 template, and a hybrid instance's actual positions reflect that template applied against a non-ideal CG backbone, not the template itself. The corrected approach converts `pe_single.gro`/`topol_aa.top` directly.

- [x] 5.1 ~~Add an AT-fragment LAMMPS-native variant of the dodecane example... Derive `dodecane_at.data` via `extract_at_frame.py`...~~ Built `examples/pe-lammps/` instead (fully LAMMPS-native: both `cg_system` and `molecules[].source`), chosen over dodecane both because the user asked for a PE example specifically and because PE's richer topology (angles, 2 AT types) exercises more of the reader. `pe_at.data`/`in.pe_at` generated by converting `pe_single.gro`/`topol_aa.top` directly via `backmap_prep.parsers.parse_gro`/`parse_top` + `backmap_prep.units`, not via extraction. `pe_cg.data` derived via `prepare_cg.py` (legitimate for CG, per the CG-side change) from a freshly-built hybrid.
- [x] 5.2 `beads[].atoms` and `cross_interactions` atom references (not bead references) rewritten from `PE:Cn` to `PE:n` via a one-off transform script (regex on the parsed YAML, not hand-editing the ~600-line settings file).
- [x] 5.3 Ran `backmap-prep` on both `examples/pe/settings.yaml` and `examples/pe-lammps/settings.yaml`; diffed the resulting hybrid output. Result: 0/1500 atom position/type/charge mismatches, `in.pe` coefficients match exactly modulo the expected `.table` filename convention difference. 79/1500 atoms have a different (but self-consistent, physically-equivalent) PBC image-flag assignment — see `design.md` Risks for the root cause and why it's not a defect.
- [x] 5.4 **Dihedral coverage gap — confirmed structural, not fixable by choosing a different existing example** (see `design.md` correction note): PE's beads are 2 atoms each, same as dodecane, so no dihedral is ever fully intra-bead for either. No linear-engine example in this repo has a 3+-atom bead. Coverage for `dihedral_coeff`/`ryckaert` parsing stays at the unit/integration level (task 2.3/3.7's synthetic fixtures) and is disclosed as a known limitation in `examples/pe-lammps/README.md`, not silently left uncovered.
- [x] 5.5 Added `tests/test_examples.py::TestLammpsNativeAtFragmentParity`, normalizing the type-name/table-filename difference and comparing the `Atoms` section with trailing image-flag tokens stripped (documented as benign, not swept under the rug); added `pe-lammps` to the shared `EXAMPLES` list.

## 6. Docs

- [x] 6.1 Update `docs/settings-reference.md`: `molecules[].source.format: lammps`, the numeric-atom-ID convention for `beads[].atoms`, the bounded `input_script` command subset.
- [x] 6.2 Add a note to `docs/tutorial-new-system.md` Step 1/2, parallel to the CG-side note already added in Step 3.
- [x] 6.3 Update `docs/cli/backmap-prep.md`, `README.md`, `docs/index.md`, `CHANGELOG.md` — same set of files the CG-side change had to catch up on a turn after the initial implementation; do them together this time instead of in a follow-up pass.

## 7. Spec/paper follow-up

- [x] 7.1 Merged the "Molecules section — CG-AT mapping", "Source file parsing", and "Feature phasing" requirement changes into `openspec/specs/backmap-input-generator/spec.md`. `openspec archive` hit the same pre-existing header-matching issue as the CG-side change (confirmed same root cause: base spec file missing the `## Purpose`/`## Requirements` wrapper); merged manually as before.
- [x] 7.2 Considered — deferred, consistent with the CG-side decision: this is an implementation detail the paper's current text doesn't need to enumerate for this revision; worth a one-clause future-work mention if that section is revisited for other reasons.
