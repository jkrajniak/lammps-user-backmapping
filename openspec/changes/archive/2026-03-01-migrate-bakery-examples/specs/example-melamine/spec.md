## ADDED Requirements

### Requirement: Melamine example with triangular CG topology

The `examples/melamine/` directory SHALL contain a melamine-formaldehyde backmapping example with a simple 3-bead triangular CG topology. The example SHALL include:

- `settings.yaml` defining the MF molecule with 3 A-type CG beads, each mapping 9 AT atoms
- AT source files: `single_mf.gro` (single-molecule reference), `topol_aa.top` (AT topology)
- CG source files: `cg_conf.gro` (CG coordinates for multiple molecules), `topol_cg.top` (CG topology)
- Tabulated CG pair potentials: `table_A_A.xvg`
- Tabulated CG bond potentials: `table_b1.xvg`
- Generated LAMMPS files: `melamine.data`, `in.melamine`
- `README.md` with system description, mapping table, and run instructions

The MF molecule has 3 CG beads arranged in a triangle (A1-A2-A3, with bonds A1-A2, A2-A3, A3-A1). Each bead maps 9 atomistic atoms (N, C, O, H atoms of the melamine-formaldehyde monomer unit).

#### Scenario: Melamine settings.yaml is valid
- **WHEN** `backmap-prep examples/melamine/settings.yaml` is run
- **THEN** it SHALL produce `melamine.data` and `in.melamine` without errors

#### Scenario: Melamine triangular topology
- **WHEN** the melamine settings define 3 beads (A1, A2, A3) with cross bonds [A1-A2], [A2-A3], [A3-A1]
- **THEN** the LAMMPS data file SHALL contain 3 CG-CG cross bonds per molecule forming a triangle

#### Scenario: Melamine molecule mapping
- **WHEN** the melamine settings define each bead with 9 atoms
- **THEN** each molecule in the data file SHALL contain 3 CG atoms + 27 AT atoms = 30 atoms total

### Requirement: Reduced melamine system for testing

The melamine example SHALL use a reduced system size suitable for quick testing. The CG coordinates file SHALL contain 10-50 molecules (vs. the original bakery 500). The `README.md` SHALL document the original system size and how to scale up.

#### Scenario: Small melamine system
- **WHEN** the melamine example CG coordinates contain 50 molecules
- **THEN** the total system SHALL have 150 CG beads + 1350 AT atoms = 1500 atoms

#### Scenario: README documents scaling
- **WHEN** a user reads `examples/melamine/README.md`
- **THEN** it SHALL state the original bakery system size (500 molecules) and explain that larger systems can be created by replacing `cg_conf.gro`

### Requirement: Melamine directory follows standard structure

The melamine example SHALL follow the same directory structure convention as the PE examples:

```
examples/melamine/
├── README.md
├── settings.yaml
├── single_mf.gro
├── topol_aa.top
├── cg_conf.gro
├── topol_cg.top
├── table_A_A.xvg
├── table_b1.xvg
├── melamine.data
└── in.melamine
```

#### Scenario: Melamine example is self-contained
- **WHEN** a user clones the repository and navigates to `examples/melamine/`
- **THEN** all files needed to run `backmap-prep settings.yaml` and then `lmp -in in.melamine` SHALL be present
