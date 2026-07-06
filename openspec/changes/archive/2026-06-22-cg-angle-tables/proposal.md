# CG angle tables (PR2a)

## Why

Rim135 hybrid topology contains 200 `[ cross_angles ]` entries with GROMACS func 8
(`cg_bonded`, tabulated CG potentials). Espresso++ maps these to `table_a{tablenr}.xvg`
via `TabulatedAngleInteractionType`.

The current LAMMPS export assigns all 200 angles to a single type with
`angle_coeff … cg 0.0 0.0` — zero restoring torque during the backmapping ramp.
Tier B dynamics blowups are expected until real tabulated angle forces are wired.

## What changes

- Add `angle_style backmap/table` C++ style (lambda-weighted tabulated angles).
- Extend `backmap-prep` to resolve func-8 angles to `table_a*.xvg` and emit
  hybrid `angle_style` input.
- Add angle-specific XVG → LAMMPS table conversion (degrees, not bond distance).

## Capabilities

### New capabilities

- `angle-backmap-table` — C++ style + generator routing for CG tabulated angles.

### Modified capabilities

- `bonded-backmap-styles` — document `angle_style backmap/table`.
- `backmap-input-generator` — func-8 angle table routing and hybrid writer.
- `integration-testing` — rim135 assertions for angle tables.

## Non-goals

- Dihedrals (PR2b)
- 1–4 cross pairs (PR2c)
- Reclassifying intra-bead AT angles as static `harmonic` (follow-up)

## Impact

- **C++:** `src/angle_backmap_table.{h,cpp}` (new)
- **Python:** `builder.py`, `network/lammps_builder.py`, `table_converter.py`, `writers.py`
- **Tests:** `test_table_converter.py`, `test_network.py`, `test_writers.py`
- **Blocks:** phase3 task 3.4 (Tier B LAMMPS smoke) until complete
