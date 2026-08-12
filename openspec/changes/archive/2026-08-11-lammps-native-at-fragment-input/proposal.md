## Why

`lammps-native-cg-input` (archived `2026-08-11-lammps-native-cg-input`) let the CG side of `backmap-prep` be supplied as a native LAMMPS `data` file instead of GROMACS `.gro`/`.top`. AT fragments (`molecules[].source`) were explicitly deferred at the time, blocked on an assumption that AT-fragment bead-mapping's name-keyed atom references (`atoms: ["1:DOD:C1", ...]`) would need a new identity-mapping sidecar, since a LAMMPS `data` file has no atom-name field.

Re-examining that assumption (see `design.md`): `builder.py` only ever uses the last colon-separated segment of an `atoms:` entry as an opaque lookup key — it never requires that key to be a symbolic name. The same numeric-ID convention already adopted for `cg_system.format: lammps`'s bead `type` values (`type: "1"` instead of `type: A`) applies here too: a LAMMPS-native AT fragment can reference atoms by their numeric LAMMPS atom ID (`atoms: ["1:DOD:1", "1:DOD:2"]`), with no new sidecar file needed.

Unlike the CG side, the AT fragment's bonded/pair force field is genuinely consumed from the topology, not just its coordinates — `at_mol.bonds`/`.angles`/`.dihedrals` are read to build the actual intra-bead bond/angle/dihedral coefficients (`builder.py:340-620,776-940`), and `at_top.atom_types[...].sigma/epsilon` drive AT-AT pair coefficients. So a LAMMPS-native AT-fragment reader needs a real (if narrowly bounded) force-field-coefficient reader, which the CG-side change correctly found unnecessary for `cg_system`.

## What Changes

- Add `molecules[].source.format: lammps` as an alternative to `format: gromacs` (default), alongside the existing `coordinates`/`topology` fields.
- Add a LAMMPS `data`-file reader for AT fragments: box (unused — AT fragment positions are relative to a single template, not a periodic system), per-atom type/charge/mass (`Masses`, `Atoms # full`), and `Bonds`/`Angles`/`Dihedrals` connectivity by type ID — unlike the CG reader, these ARE read here.
- Add a bounded LAMMPS input-script reader covering exactly the four coefficient families `builder.py` consumes for AT fragments: `bond_coeff` (`bond_style harmonic`), `angle_coeff` (`angle_style harmonic`), `dihedral_coeff` (`dihedral_style ryckaert`, 6 Ryckaert–Bellemans-form coefficients — the package's own native style, so no unit-conversion math is needed on read, unlike the GROMACS RB→LAMMPS conversion path), and `pair_coeff <i> <i>` self-terms (LJ epsilon/sigma per AT atom type; cross-terms remain computed in Python via the existing arithmetic/geometric mixing already used for the GROMACS path — no mixing-rule concept needed on the input side).
- Adopt the numeric-atom-ID convention for `beads[].atoms` entries when the owning molecule's `source.format: lammps` (mirrors the CG side's numeric type-ID convention for `beads[].type`).
- Assume `units real`, no conversion — same convention as the CG-side change.
- Explicitly out of scope: GROMACS `[virtual_sites3]`-equivalent support (no LAMMPS `data`-file construct maps onto it directly), degree-dependent multi-file AT sources (`atoms_by_degree`/Phase 3 reactive networks), and the network engine generally — mirrors the CG-side change's linear-engine-only scope.

## Capabilities

### New Capabilities
(none — this extends an existing capability)

### Modified Capabilities
- `backmap-input-generator`: `molecules[].source` gains a `format: lammps` option (`data` + `input_script` fields); "Molecules section — CG-AT mapping" requirement's atom-reference format gains the numeric-ID convention for LAMMPS-native fragments; "Source file parsing" requirement gains a scoped LAMMPS AT-fragment reader; "Feature phasing" requirement is updated to reflect AT-fragment LAMMPS input as implemented (linear engine, single-file sources only).

## Impact

- `python/src/backmap_prep/schema.py`: `SourceFiles` gains `format: Literal["gromacs", "lammps"]`, `data: str | None`, `input_script: str | None`; validation requiring `data`+`input_script` for `lammps`, `coordinates`+`topology` for `gromacs`; reject `format: lammps` combined with degree-dependent (`list[DegreeSourceFile]`) sources or `prep.engine: network`.
- `python/src/backmap_prep/parsers/lammps_data_parser.py`: extend with an AT-fragment variant reading `Bonds`/`Angles`/`Dihedrals` sections (which the CG reader intentionally skips), producing `TopBond`/`TopAngle`/`TopDihedral` objects.
- New `python/src/backmap_prep/parsers/lammps_script_parser.py` (the module the CG-side change scoped out and didn't need): a bounded reader for `bond_coeff`/`angle_coeff`/`dihedral_coeff`/`pair_coeff`.
- `python/src/backmap_prep/builder.py`: branch AT-fragment loading on `mol_def.source.format`; no change to the AT bond/angle/dihedral coefficient *consumption* logic itself, since the reader produces the same `TopBond`/`TopAngle`/`TopDihedral`/`AtomType` shapes.
- Docs: `docs/settings-reference.md`, `docs/tutorial-new-system.md`, `docs/cli/backmap-prep.md`.
- Tests: a new `examples/dodecane-lammps-cg-fragment/` (or similarly named) example — likely built the same way as `examples/dodecane-lammps-cg/`, by extracting a LAMMPS-native AT fragment `data`/input-script pair from an already-built hybrid/AT-reference system — plus round-trip parity tests against `examples/dodecane/`.
