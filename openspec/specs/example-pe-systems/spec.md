## ADDED Requirements

### Requirement: PE example with 2:1 united-atom mapping

The `examples/pe/` directory SHALL contain a complete polyethylene backmapping example using OPLS united-atom force field with 2 UA atoms per CG bead. The example SHALL include:

- `settings.yaml` defining the PE molecule with 2 bead types (A for terminal, B for backbone), each mapping 2 UA carbon atoms
- AT source files: `pe_single.gro` (single-chain reference), `topol_aa.top` (AT topology with bond/angle params)
- CG source files: `cg_conf.gro` (CG coordinates), `topol_cg.top` (CG topology with tabulated bonds)
- Tabulated potential files: `table_b1.xvg` and `table_b2.xvg` (CG bonds), `table_a1.xvg` and `table_a2.xvg` (CG angles), `table_A_A.xvg`, `table_A_B.xvg`, `table_B_B.xvg` (CG pair tables)
- Generated LAMMPS files: `pe.data`, `in.pe`
- `README.md` with system description, mapping table, and run instructions

The CG model SHALL have 50 beads per chain (2 A + 48 B) with each bead mapping 2 UA carbon atoms. Cross interactions SHALL include AT harmonic bonds between adjacent beads and CG tabulated bonds.

#### Scenario: PE settings.yaml is valid
- **WHEN** `backmap-prep examples/pe/settings.yaml` is run
- **THEN** it SHALL produce `pe.data` and `in.pe` without errors

#### Scenario: PE molecule mapping
- **WHEN** the PE settings define bead A1 with atoms `[1:PE:C1, 1:PE:C2]` and bead B1 with atoms `[1:PE:C3, 1:PE:C4]`
- **THEN** each molecule in the data file SHALL contain 50 CG atoms + 100 AT atoms = 150 atoms total

#### Scenario: PE cross interactions include dihedrals
- **WHEN** the PE settings define cross bonds, cross angles, and cross dihedrals between adjacent beads
- **THEN** the generated LAMMPS input SHALL use `bond_style hybrid`, `angle_style hybrid`, and `dihedral_style hybrid` with `backmap/*` sub-styles

### Requirement: PE4 example with 4:1 united-atom mapping

The `examples/pe4/` directory SHALL contain a polyethylene backmapping example with 4 UA atoms per CG bead. The example SHALL follow the same structure as the PE example but with coarser mapping:

- 25 CG beads per chain (2 A + 23 B), each bead mapping 4 UA carbon atoms
- AT source files with 100-atom single chain reference
- CG source files with 25-bead chains
- Tabulated CG bond/angle/dihedral potentials
- `settings.yaml`, generated LAMMPS files, `README.md`

#### Scenario: PE4 settings.yaml is valid
- **WHEN** `backmap-prep examples/pe4/settings.yaml` is run
- **THEN** it SHALL produce `pe4.data` and `in.pe4` without errors

#### Scenario: PE4 molecule has 4 atoms per bead
- **WHEN** the PE4 settings define bead A1 with 4 atoms and bead B1 with 4 atoms
- **THEN** each molecule SHALL contain 25 CG atoms + 100 AT atoms = 125 atoms total

### Requirement: PE-10 example with all-atom 10:1 mapping

The `examples/pe_10/` directory SHALL contain a polyethylene backmapping example using OPLS all-atom force field with approximately 30-36 atoms per CG bead (10 monomers of CH2 per bead, with explicit hydrogens). The example SHALL include:

- `settings.yaml` defining the PE molecule with A terminal beads (~30 atoms each including H) and B backbone beads (~36 atoms each)
- AT source files with all-atom topology including C-H bonds and H-C-H angles
- CG source files with 10-bead chains
- Tabulated CG bond/angle potentials
- Generated LAMMPS files, `README.md`

This example validates `backmap-prep` with a high AT-to-CG ratio and OPLS/AA force field parameters.

#### Scenario: PE-10 settings.yaml is valid
- **WHEN** `backmap-prep examples/pe_10/settings.yaml` is run
- **THEN** it SHALL produce `pe_10.data` and `in.pe_10` without errors

#### Scenario: PE-10 molecule has many atoms per bead
- **WHEN** the PE-10 settings define bead B1 with ~36 atoms
- **THEN** the data file SHALL contain the correct number of intra-CG bonds and angles for all C-C, C-H, C-C-C, H-C-H, and C-C-H interactions within each bead

### Requirement: PE-AA example with explicit hydrogen 2:1 mapping

The `examples/pe_aa/` directory SHALL contain a polyethylene all-atom backmapping example using OPLS/AA force field where each bead maps 2 heavy carbons plus their attached hydrogens (6-7 atoms per bead). The example SHALL include:

- `settings.yaml` defining the PE molecule with A terminal beads (7 atoms: 2C + 5H) and B backbone beads (6 atoms: 2C + 4H)
- AT source files with full OPLS/AA topology (opls_136/opls_140 types)
- CG source files with 50-bead chains
- Tabulated CG potentials
- Generated LAMMPS files, `README.md`

#### Scenario: PE-AA settings.yaml is valid
- **WHEN** `backmap-prep examples/pe_aa/settings.yaml` is run
- **THEN** it SHALL produce `pe_aa.data` and `in.pe_aa` without errors

#### Scenario: PE-AA includes hydrogen atoms in mapping
- **WHEN** the PE-AA settings define bead A1 with atoms including H1, H2, H3 (methyl hydrogens) and C1, C2
- **THEN** each bead's atom list SHALL include both carbon and hydrogen atoms with correct OPLS/AA types and charges

### Requirement: Consistent example directory structure

All PE examples SHALL follow a consistent directory structure:

```
examples/<name>/
├── README.md              # System description, mapping table, how to run
├── settings.yaml          # backmap-prep configuration
├── <single>_single.gro    # Single-chain AT reference coordinates
├── topol_aa.top           # AT topology (bonds, angles, types)
├── cg_conf.gro            # CG coordinates (multi-chain)
├── topol_cg.top           # CG topology (tabulated bonds)
├── table_*.xvg            # CG tabulated potentials (GROMACS format)
├── <name>.data            # Generated LAMMPS data file
└── in.<name>              # Generated/reference LAMMPS input script
```

#### Scenario: Example directory is self-contained
- **WHEN** a user clones the repository and navigates to `examples/pe/`
- **THEN** all files needed to run `backmap-prep settings.yaml` and then `lmp -in in.pe` SHALL be present in that directory

#### Scenario: README documents the system
- **WHEN** a user reads `examples/pe/README.md`
- **THEN** it SHALL contain: system description (molecule, CG model, AT model), a mapping table (bead → atoms), quick-start commands, and a file listing with descriptions

### Requirement: Source files copied from bakery with consistent naming

Source GROMACS files SHALL be copied from the bakery `examples/` directories and renamed to follow consistent conventions:

| Bakery name | New name | Purpose |
|-------------|----------|---------|
| `pe_single.gro` or equiv. | `pe_single.gro` | Single-chain AT reference |
| `topol.top` or `topol_at.top` | `topol_aa.top` | AT topology |
| `cg_conf.gro` or `conf_cg.gro` | `cg_conf.gro` | CG coordinates |
| `cg_topol.top` | `topol_cg.top` | CG topology |
| `table_*.xvg` | `table_*.xvg` | CG potentials (keep original names) |

#### Scenario: File renaming for consistency
- **WHEN** bakery has `topol.top` (PE example) and `topol_at.top` (PE-10 example)
- **THEN** both SHALL be renamed to `topol_aa.top` in their respective example directories
