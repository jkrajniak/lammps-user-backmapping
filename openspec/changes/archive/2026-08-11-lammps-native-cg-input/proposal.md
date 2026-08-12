## Why

`backmap-prep` currently requires the CG system (`cg_system` in `settings.yaml`) to be supplied as GROMACS `.top`/`.gro`. The paper's own motivation for porting to LAMMPS is that many target users already run production workloads in LAMMPS and don't want to install/maintain a GROMACS-adjacent toolchain — but those same users, if their CG model was built with LAMMPS-native tooling, currently have no way to feed it into `backmap-prep` without first producing a GROMACS-style topology they may never otherwise need. `backmap-input-generator`'s existing "Source file parsing" requirement already reserves a "Phase 4 — Format flexibility" slot for `cg_system.format: lammps`, including a passthrough scenario for native `.table` files, but it was never implemented (`schema.py` still pins `CGSystem.format` to `Literal["gromacs"]`) and was never scoped in enough detail to build.

This change scopes and implements that deferred slot for the CG side only. The CG side of the pipeline already matches atoms positionally and by type (`builder.py:683-729`), not by the free-text atom names GROMACS `.gro`/`.top` carry — unlike the AT-fragment side, which is driven by name-keyed bead mappings (`atoms: [C1, H1, ...]`) that a LAMMPS `data` file cannot supply without inventing a new naming sidecar. Scoping this change to the CG system avoids that unresolved design question entirely and delivers a self-contained, independently useful capability.

## What Changes

**Revised during implementation** (see `design.md` "Correction" note): tracing `builder.py` shows the CG topology's own bonded/pair sections (`cg_mol.bonds`/`.angles`) are never read — CG-CG bonded terms come entirely from `cross_interactions` in the YAML, format-independent already. So no separate LAMMPS input-script coefficient reader is needed; a `data`-file reader alone is sufficient and this change is smaller than originally scoped.

- Add `cg_system.format: lammps` as a supported alternative to `format: gromacs`.
- Add a LAMMPS `data`-file reader for the CG system: box, per-atom type/mass (from `Masses`)/charge/position (from `Atoms # full`), matched positionally per molecule exactly as the existing GROMACS `cg_gro`/`cg_top` path is today. `Bonds`/`Angles` sections, if present, are tolerated but not read (consistent with `cg_mol.bonds`/`.angles` already being unused for the GROMACS path).
- No unit conversion is applied to LAMMPS-native CG values (positions, box): the `data` file is assumed to already be in LAMMPS `real` units, matching how `write_data` produces it.
- Extend the existing CG pair-table lookup (`table_groups`) to also try a `.table`-extension file directly before falling back to `.xvg` — a small, format-agnostic improvement that also benefits the GROMACS path, since dodecane's example already ships both `.xvg` and pre-converted `.table` files.
- Update the "Feature phasing" requirement: this change delivers the CG-input half of the previously-unscoped Phase 4 bucket; AT-fragment LAMMPS-native input remains explicitly deferred pending a separate design for atom-identity naming.

## Capabilities

### New Capabilities
(none — this extends an existing capability)

### Modified Capabilities
- `backmap-input-generator`: `cg_system` section gains a `format: lammps` option; "Source file parsing" requirement gains a scoped (not aspirational) LAMMPS CG `data`-file reader; "Feature phasing" requirement is updated to reflect this slice of Phase 4 as implemented, with AT-fragment LAMMPS input still deferred.

## Impact

- `python/src/backmap_prep/schema.py`: `CGSystem.format` becomes `Literal["gromacs", "lammps"]`; `coordinates`/`topology` become optional; new optional `data` field for the LAMMPS-native path.
- `python/src/backmap_prep/parsers/`: new `lammps_data_parser.py`, producing `gro_parser.GroFile`/`top_parser.Topology`-compatible objects (reusing those dataclasses directly) so `builder.py` needs minimal changes.
- `python/src/backmap_prep/builder.py`: branch `cg_gro`/`cg_top` construction on `cg_system.format`; skip unit conversion at the (exactly three) CG-position/box call sites for the LAMMPS-native path; extend the `table_groups` pair-table lookup to try `.table` before `.xvg`.
- `python/src/backmap_prep/units.py`: not used on the LAMMPS-native CG path (values already in `real` units) — no changes to the module itself, just skipped at the relevant call sites.
- Docs: `docs/cli/backmap-prep.md` and the settings-format reference need a new worked example (a LAMMPS-native CG system paired with GROMACS AT fragments, since AT fragments stay GROMACS-only in this change).
- Tests: new `examples/dodecane-lammps-cg/` directory, built from the LAMMPS-native CG system that `examples/dodecane/prepare_cg.py` + `in.dodecane_cg_equil` already produce today (currently used only for a standalone post-hoc re-equilibration check, not as `backmap-prep` input). Round-trip validation runs `backmap-prep` on both `examples/dodecane/settings.yaml` (GROMACS CG path) and `examples/dodecane-lammps-cg/settings.yaml` (LAMMPS-native CG path) and diffs the resulting hybrid output.
