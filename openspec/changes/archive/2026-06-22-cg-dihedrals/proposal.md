# CG dihedrals (PR2b)

## Why

Rim135 hybrid topology contains **29,821** `[ dihedrals ]` and **3,600** `[ cross_dihedrals ]`
(func 3 Ryckaert–Bellemans, no inline coeffs). Export writes `0 dihedrals`. PE examples
define tabulated CG dihedrals (`table_d*.xvg`) but also export zero dihedrals.

Without dihedral torques, Tier B dynamics remain unphysical even after PR2a angle tables.

## What changes

- Add `dihedral_style ryckaert`, `backmap/ryckaert`, `backmap/table` C++ styles.
- Extend `backmap-prep` to resolve func-3 RB from OPLS dihedraltypes, export hybrid dihedrals.
- Add `table_d*.xvg` conversion for PE CG tabulated dihedrals.

## Non-goals

- PR2c cross 1–4 pairs (`[ cross_pairs ]`)

## Impact

- **C++:** `dihedral_ryckaert`, `dihedral_backmap_ryckaert`, `dihedral_backmap_table`
- **Python:** `top_parser`, `builder`, `network/lammps_builder`, `table_converter`, `writers`
- **Blocks:** phase3 task 3.4 Tier B until complete (with PR2a angles)
