## Context

`backmap-prep` builds a hybrid CG/AT LAMMPS system from four logically independent inputs: CG coordinates, CG force field, AT fragment coordinates, and AT fragment force field. Today all four are read from GROMACS `.gro`/`.top` via `parsers/gro_parser.py` and `parsers/top_parser.py`. The output side (`writers.py`) already writes LAMMPS-native `data` files and input scripts using a small, fixed set of styles (`bond_style harmonic|table`, `angle_style harmonic|table`, `dihedral_style ryckaert|harmonic|charmm|table`, `pair_style backmap ... lj/cut/coul/cut`, plus CG tables converted from `.xvg` to `.table`).

The CG side of `builder.py` matches atoms by position and type, not by name: `cg_start = mol_idx * cg_atom_count`, then walks `cg_gro.atoms[gro_idx]` in molecule order, looking up type/mass/charge in `cg_top.atom_types` (`builder.py:683-729`). This is the key fact that makes CG-side LAMMPS-native input tractable without new naming infrastructure — a LAMMPS `data` file's `Atoms` section is exactly this shape already (ordered, per-atom type ID, charge, position).

The AT-fragment side is matched by free-text atom name (`atom.name`, used at `builder.py:322,343,457,564,745,813` etc.) sourced from `.gro`'s name column and `.top`'s atom name field — a LAMMPS `data` file has no equivalent field. Supporting LAMMPS-native AT fragments would require a new identity-mapping sidecar and is explicitly out of scope for this change (tracked as a future follow-up).

The existing `backmap-input-generator` spec already reserves `cg_system.format: lammps` as an aspirational "Phase 4" scenario, but it was never scoped precisely enough to implement (no decision on mixing rules, units enforcement, or how force-field coefficients — which live in a LAMMPS *input script*, not the `data` file — get supplied at all).

> **IMPORTANT CORRECTION (found during implementation):** the paragraph above, and the original Goals/Decisions below, assumed CG-CG bonded/nonbonded coefficients flow through `cg_system.topology` the same way AT-fragment coefficients flow through `molecules[].source.topology` (`at_mol.bonds`/`.angles`/`.dihedrals` are actively read at `builder.py:340-620,776-940`). Tracing the actual CG-side code path shows this is false: `cg_mol.bonds`/`cg_mol.angles`/`cg_mol.dihedrals` (from `cg_top`) are **never read anywhere in `builder.py`** (confirmed by exhaustive grep for `cg_mol.`/`cg_top.` usage — only `cg_top.molecules`, `cg_mol.atoms`, and `cg_top.atom_types` for a mass fallback are used). CG-CG *bonded* terms come entirely from `cross_interactions` in the YAML settings (already format-independent — e.g. dodecane's `settings.yaml` hand-specifies `cg_bonded: true` pairs with GROMACS-func-style inline params, sourced from neither `.top` nor `.data`). CG-CG *nonbonded* terms come from the `table_groups` naming convention (`table_{type_a}_{type_b}.xvg`, `builder.py:656-676`), keyed on `AtomTypeInfo.name` and files present in `base_dir` — also independent of `cg_system.topology`'s own (parsed-but-unused) `[atomtypes]` sigma/epsilon fields.
>
> Consequence: a LAMMPS-native CG `data` file alone — box, `Masses`, `Atoms # full` — supplies everything `builder.py` actually consumes from the CG side. No separate input-script coefficient reader is needed, and the "two input files" decision, the bounded input-script grammar, the `units real`-in-script requirement, and the mixing-rule decision below are **not implemented as originally written**. They're kept below, marked, rather than deleted, because they document a real alternative that was seriously considered and could become relevant if `builder.py` is later extended to actually consume CG-CG bonded coefficients from `cg_system` directly (see revised Open Questions). The *positions/box* half of the "no unit conversion" decision still stands and is implemented — see "Decisions (as implemented)" below.

## Goals / Non-Goals

**As implemented** (see correction above) — narrower than originally stated:
- Let a user supply CG *coordinates, box, per-atom type, mass, and charge* as a LAMMPS `data` file, instead of GROMACS `.gro`/`.top`.
- No unit conversion on these values — a `data` file produced by `write_data` under `units real` is already in the target units.
- CG-CG bonded terms (via `cross_interactions`) and nonbonded terms (via `table_groups`) are configured exactly as they are today, unchanged, regardless of `cg_system.format`.

**Original goals (not implemented — see correction):**
- Let a user supply the CG system as a LAMMPS `data` file (coordinates, box, per-atom type/charge/mass, bond/angle connectivity) plus a LAMMPS input-script fragment (styles + coefficients), instead of GROMACS `.top`/`.gro`.
- Support both analytic (`harmonic`, `lj/cut`) and tabulated (`table`) CG potentials, since tabulated IBI-derived potentials are the common case in practice.
- Pass `.table` files through unchanged — no `.xvg`-equivalent conversion needed for a format that's already native.
- Keep the existing GROMACS CG path, `builder.py`'s positional/type-keyed matching logic, and AT-fragment handling completely unchanged.

**Non-Goals:**
- LAMMPS-native AT fragment input (needs the atom-identity naming sidecar — separate change).
- Support for arbitrary LAMMPS pair/bond/angle/dihedral styles beyond the bounded set the writer already emits.
- Support for LAMMPS unit styles other than `real`.
- PDB/XYZ CG input (also nominally in the old Phase 4 bucket, not addressed here).

## Decisions

### As implemented

**Decision: one input file (`data`), not two.** Since `builder.py` never reads `cg_top`'s bonded/pair sections (see correction above), a `data` file's `Masses` + `Atoms # full` sections are sufficient. `cg_system` for `format: lammps` takes a single `data` field, mirroring `format: gromacs`'s `coordinates`+`topology` pair only in spirit, not in shape.

**Decision: `data`-file values are read as-is, no unit conversion.** `builder.py` applies `units.distance()` at exactly three CG-specific call sites: the box (`builder.py:207-209`) and the two CG-position blocks (`builder.py:700-702,728-730`). For `format: lammps` these three sites use an identity pass-through instead, since a `data` file written under `units real` (the convention this project's own `write_data` output already uses) needs no conversion. Mass and charge already have an identity conversion factor in `units.py` (`MASS = CHARGE = 1.0`) so they need no special-casing.

**Decision: extend the `table_groups` CG pair-table lookup to try `.table` before `.xvg`.** This is unrelated to `cg_system.format` — it's a small, independently-useful fix to `builder.py:656-676`, which today only ever constructs an `.xvg` filename even when a `.table` file of the same name already exists (as it does in `examples/dodecane/`). Format-agnostic: benefits GROMACS-sourced and LAMMPS-native CG systems equally.

**Decision: `cg_atom.type` for a LAMMPS-native CG system is the LAMMPS numeric type ID, as a string (e.g. `"1"`, `"2"`).** `builder.py` only ever uses `atom.type` as an internal, distinct dict key (`type_map`) and a human-readable label (`AtomTypeInfo.name`) — it is never matched against an external naming convention for CG atoms (unlike the AT-fragment side). So the numeric type ID from the `data` file's `Masses`/`Atoms` sections can be used directly, with no naming sidecar. The one consequence: `molecules[].beads[].type` in `settings.yaml` must then reference that same numeric ID as a string (e.g. `type: "1"`) instead of a symbolic name — document this in `docs/cli/backmap-prep.md`.

### Original decisions (not implemented — kept for record, see correction above)

**Decision: two input files, not one.** A LAMMPS `data` file alone cannot supply force-field coefficients (those are `pair_coeff`/`bond_coeff`/`angle_coeff` commands, conventionally kept in the input script, not the data file, in most LAMMPS workflows this project's users will already have). Rather than requiring users to restructure their coefficients into the `data` file's optional `Pair Coeffs`/`Bond Coeffs` sections, `cg_system` for `format: lammps` takes two paths: a `data` file and an `input_script` fragment. This matches how LAMMPS users actually keep their files, and it mirrors the existing `coordinates`/`topology` split already used for `format: gromacs`.
*Alternative considered*: require coefficients embedded in the `data` file's optional sections only. Rejected — forces an unnatural file layout on users and doesn't support tabulated potentials (`pair_style table` coefficients are file paths, not embeddable in a `data` file section).

**Decision: bound the input-script grammar to exactly the styles the writer already emits.** The reader only needs to understand `units`, `mass`, `pair_style`/`pair_coeff` (`lj/cut`, `table`, or `hybrid` combining them), `bond_style`/`bond_coeff` (`harmonic`, `table`, `hybrid`), `angle_style`/`angle_coeff` (`harmonic`, `table`, `hybrid`), and `pair_modify mix`. Any other command is either ignored (e.g. `run`, `thermo`, output commands common in a full input script a user might reuse as-is) or rejected with a clear "unsupported style" error if it appears as a `*_style` directive outside this set.
*Alternative considered*: a general LAMMPS input-script interpreter (handling variables, loops, `include`). Rejected as unbounded scope for no payoff — nothing downstream can consume styles beyond the ones already listed, so parsing them would produce data nothing can use.

**Decision: require `units real`, reject everything else outright.** The existing GROMACS path always converts nm/kJ/ps into LAMMPS `real` units (`units.py`); the LAMMPS-native path skips that conversion entirely and trusts the input values as-is. If the input script declares any unit style other than `real` (or omits `units` — LAMMPS itself would default to `lj`), abort immediately rather than silently misinterpreting values.
*Alternative considered*: auto-convert from other unit styles. Rejected — adds a second full unit-conversion table for a case with no evidenced user need, and risks producing plausible-but-wrong numbers by silent conversion.

**Decision: cross-type LJ mixing is either explicit `pair_coeff` for every unlike pair, or a single `pair_modify mix geometric|arithmetic` command — no numeric "combination rule" concept.** This mirrors how LAMMPS itself resolves mixing, rather than reusing GROMACS's `combination_rule: 1|2|3` integer, which has no LAMMPS equivalent.

**Decision: `.table`-extension files referenced from `pair_coeff`/`bond_coeff`/`angle_coeff table` are passed through byte-for-byte.** `table_converter.py` already special-cases `.xvg → .table`; this adds a no-op passthrough branch keyed on the `.table` extension, consistent with the spec's existing (unimplemented) passthrough scenario.

**Decision: scope this to the CG system only, not AT fragments.** See Context — the CG side has no naming dependency, the AT-fragment side does. Bundling both into one change would gate CG-side value on solving a harder, separate problem.

## Verification example

`examples/dodecane/` already contains a working instance of a LAMMPS-native CG `data` file, produced by an existing but differently-directed workflow: `prepare_cg.py` extracts a standalone pure-CG `data` file (box, `Masses`, `Atoms # full`, plus a `Bonds` section that the CG reader will tolerate but not read) from the hybrid output. Today this only feeds a standalone post-backmap re-equilibration check (`in.dodecane_cg_equil`), one direction away from what this change needs. `in.dodecane_cg_equil` itself is *not* consumed by `backmap-prep` under the as-implemented design — it remains useful only as independent evidence that the extracted `data` file is physically valid (LAMMPS can `read_data` it and run), and as documentation for what a real CG-equilibration script for this system looks like.

The verification example (see `tasks.md` §6) closes the loop the *implemented* way: feed the extracted `data` file back into `backmap-prep` as `cg_system.data`, alongside the unchanged `cross_interactions`/`table_groups` configuration, and confirm the resulting hybrid build is numerically identical to the GROMACS-path build for the same system.

## Risks / Trade-offs

- **[Risk] A `data` file with a different atom/molecule ordering convention than the one `builder.py`'s positional matching assumes (contiguous per-molecule blocks) silently produces a scrambled hybrid system** → Mitigation: this is the same risk the GROMACS `.gro` path already has (no name-based cross-check exists there either); document the ordering assumption explicitly in `docs/cli/backmap-prep.md`; the round-trip verification example is the concrete check that catches a regression here.
- **[Risk] `units real` assumption is wrong for some future user (e.g. `units metal`)** → Mitigation: since there's no input script to read a `units` directive from, this is simply documented as a hard assumption of the `data` file's contents (matching what `write_data` under `units real` already produces); revisit only if a real user needs a different unit style.
- **[Trade-off] No AT-fragment LAMMPS support yet** → users mixing a LAMMPS-native CG system with LAMMPS-native AT fragments still can't go fully GROMACS-free in this change. Explicitly disclosed as a follow-up in the proposal rather than silently implied as "done."
- **[Trade-off] `molecules[].beads[].type` must reference the LAMMPS numeric type ID as a string when `cg_system.format: lammps`**, instead of a symbolic name as with GROMACS — a small but real difference in how the same YAML file section is written depending on CG format. Documented, not hidden.

## Open Questions

- Should `builder.py` eventually read CG-CG bonded coefficients directly from `cg_system` (GROMACS or LAMMPS) instead of requiring them to be duplicated in `cross_interactions`? Out of scope here, but this is where the original "input-script coefficient reader" idea (see "Original decisions" above) would become relevant again — if pursued, revisit that design rather than the `cross_interactions`-only description above, since the underlying question (bounded style grammar, `units real`, mixing rules) still applies.
- Exact wording for the "ordering assumption" documentation note, and whether to add a best-effort sanity check (e.g. atom count divisible by the expected per-molecule count) rather than only documenting the assumption.
