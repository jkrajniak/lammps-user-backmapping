# Tasks — CG angle tables (PR2a)

## 1. OpenSpec

- [x] 1.1 Create `openspec/changes/cg-angle-tables/` with proposal, design, tasks, delta specs
- [x] 1.2 Link blocker from phase3 task 3.4

## 2. C++ angle style

- [x] 2.1 Add `src/angle_backmap_table.h` and `src/angle_backmap_table.cpp`
- [x] 2.2 Register `AngleStyle(backmap/table, AngleBackmapTable)`
- [x] 2.3 Implement restart read/write for is_cg flags

## 3. Python generator

- [x] 3.1 Extend `AngleTypeInfo` with `table_file` / `table_keyword`
- [x] 3.2 Wire func-8 angles in `network/lammps_builder.py`
- [x] 3.3 Add `_convert_angle_xvg` in `table_converter.py`
- [x] 3.4 Update `writers.py` for hybrid angle styles

## 4. Tests and validation

- [x] 4.1 Unit tests for angle XVG conversion and hybrid writer
- [x] 4.2 Rim135 integration: no placeholder coeffs, table files emitted
- [x] 4.3 Rebuild rim135 and verify generated input (LAMMPS step-0 parse with rebuilt binary)
