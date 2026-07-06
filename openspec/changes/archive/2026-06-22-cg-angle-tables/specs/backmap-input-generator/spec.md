## ADDED Requirements

### Requirement: Func-8 CG angle table routing

When hybrid TOP angles have GROMACS func 8 (tabulated CG), the generator SHALL
resolve `table_a{tablenr}.xvg` from `prep.tables_dir` and assign
`angle_style backmap/table cg` types.

#### Scenario: Rim135 func-8 cross-angles

- **WHEN** `backmap-prep build examples/epoxy/settings.v2.yaml` is run
- **THEN** all func-8 cross-angles SHALL reference `backmap/table cg` types
- **AND** the generated input SHALL NOT contain `angle_coeff … cg 0.0 0.0`

#### Scenario: Missing angle table file

- **WHEN** a func-8 angle references tablenr N but `table_aN.xvg` is not found
- **THEN** the generator SHALL fail with a clear error

### Requirement: Angle table XVG conversion

Angle table files (`table_a*.xvg`) SHALL be converted to LAMMPS `.table` format
with angle-specific units: degrees for the independent variable, kcal/mol for
energy, and kcal/(mol·deg) for force.

#### Scenario: table_a1.xvg conversion

- **WHEN** `table_a1.xvg` is converted for rim135
- **THEN** column 0 SHALL remain in degrees (not nm→Å)
- **AND** energy and force columns SHALL use kcal/mol and kcal/(mol·deg)
