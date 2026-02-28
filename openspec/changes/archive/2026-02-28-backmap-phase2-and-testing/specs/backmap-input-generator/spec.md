## MODIFIED Requirements

### Requirement: Cross-interactions section

The `cross_interactions` section SHALL define cross-CG bonded interactions — bonds, angles, dihedrals, and 1-4 pairs that connect atoms in different CG beads. These are the interactions that receive lambda weighting.

```yaml
cross_interactions:
  bonds:
    - params: "1 0.153 7150000.0"
      pairs:
        - [DOD:C2, DOD:C3]
    - params: "8 1 1.0"
      table: table_b1.xvg
      pairs:
        - [DOD:A1, DOD:B1]
      cg_bonded: true
  angles:
    - params: "1 111.56 585.76"
      triples:
        - [DOD:C4, DOD:C5, DOD:C6]
  dihedrals:
    - params: "3 9.2048 12.552 -13.120 -33.472 6.2726 33.464"
      quadruples:
        - [DOD:C2, DOD:C3, DOD:C4, DOD:C5]
        - [DOD:C3, DOD:C4, DOD:C5, DOD:C6]
    - params: "8 1 1.0"
      table: table_d1.xvg
      quadruples:
        - [DOD:A1, DOD:B1, DOD:B2, DOD:B3]
      cg_bonded: true
  pairs:
    - params: "1 0.35 0.3"
      atom_pairs:
        - [DOD:C1, DOD:C4]
```

Entries without `cg_bonded: true` are AT-level cross interactions (weighted by `λ²`). Entries with `cg_bonded: true` are CG-level cross interactions (weighted by `1-λ²`).

The `dihedrals` sub-section SHALL support:
- **Function type 3** (Ryckaert-Bellemans): `params: "3 C0 C1 C2 C3 C4 C5"` — six RB coefficients in kJ/mol (GROMACS units). The generator SHALL transform coefficients to LAMMPS convention: `C_lammps[n] = (-1)^n × C_gromacs[n]` and convert units from kJ/mol to kcal/mol.
- **Function type 8** (tabulated): `params: "8 table_index scale_factor"` with `table:` field specifying the GROMACS `.xvg` table file.

The `pairs` sub-section SHALL support 1-4 non-bonded pairs with LJ parameters. These are relevant when `exclusion_nrexcl < 3`.

#### Scenario: PE cross bonds with dihedrals
- **WHEN** the YAML defines cross bonds, cross angles, and cross dihedrals for polyethylene
- **THEN** the generator SHALL assign distinct LAMMPS types: intra-CG bonds get static `harmonic`, AT cross bonds get `backmap/harmonic at`, CG cross bonds get `backmap/table cg`, AT cross dihedrals get `backmap/ryckaert at`, CG cross dihedrals get `backmap/table cg`

#### Scenario: RB coefficient conversion
- **WHEN** GROMACS RB dihedral params are `"3 9.2048 12.552 -13.120 -33.472 6.2726 33.464"` (kJ/mol)
- **THEN** the generator SHALL output LAMMPS coefficients with sign flips on odd indices: `C0=2.2000, C1=-3.0000, C2=-3.1362, C3=8.0000, C4=1.4998, C5=-8.0000` (kcal/mol, approximate)

#### Scenario: No cross interactions (water)
- **WHEN** the `cross_interactions` section is empty or absent
- **THEN** the generator SHALL use only static bond/angle styles (no `backmap/*` styles needed)

#### Scenario: Dihedral table conversion
- **WHEN** a cross dihedral specifies `table: table_d1.xvg` with `cg_bonded: true`
- **THEN** the generator SHALL convert the `.xvg` file to LAMMPS dihedral table format with angles in degrees and energies in kcal/mol

### Requirement: LAMMPS data file generation

The generator SHALL produce a LAMMPS data file (`atom_style full`) containing:
- Header with atom/bond/angle/dihedral counts and type counts
- Box dimensions (converted from nm to Angstrom)
- Masses section (one entry per atom type)
- Atoms section with: atom_id, molecule_id, type, charge, x, y, z
- Bonds section with bond topology
- Angles section with angle topology
- Dihedrals section (if cross-CG dihedrals or intra-CG dihedrals exist)

Dihedral type assignment SHALL follow the same three-category scheme as bonds and angles:
- **Intra-CG dihedrals** (from AT source `[ dihedrals ]` where all four atoms belong to the same CG bead): assigned to static dihedral types (e.g., `multi/harmonic`)
- **Cross-CG AT dihedrals** (from `cross_interactions.dihedrals` without `cg_bonded`): assigned to `backmap/ryckaert at` or `backmap/table at` types
- **Cross-CG CG dihedrals** (from `cross_interactions.dihedrals` with `cg_bonded: true`): assigned to `backmap/ryckaert cg` or `backmap/table cg` types

#### Scenario: Dodecane data file with dihedrals
- **WHEN** generating data for dodecane molecules with cross dihedrals defined
- **THEN** the Dihedrals section SHALL contain entries for each cross-CG dihedral per molecule, with correct atom IDs and dihedral type

#### Scenario: Data file header includes dihedral counts
- **WHEN** the system has N_d total dihedrals and T_d dihedral types
- **THEN** the data file header SHALL include `N_d dihedrals` and `T_d dihedral types`

### Requirement: LAMMPS input script generation

The generator SHALL produce a LAMMPS input script (`.in`) configured for backmapping. The script SHALL include:
- `units real` and `atom_style full`
- `read_data` for the generated data file
- `pair_style backmap` with AT and CG sub-styles
- `pair_coeff` for all type pairs
- Bond/angle/dihedral styles:
  - If no cross interactions: static styles only
  - If cross interactions exist: `bond_style hybrid ...`, `angle_style hybrid ...`, `dihedral_style hybrid ...`
- Dihedral style routing: `dihedral_style hybrid backmap/ryckaert backmap/table` (only include sub-styles that are actually used)
- `dihedral_coeff` for each dihedral type, routing to the correct sub-style
- Group definitions for AT and CG atoms
- `fix backmap` with parameters from `simulation` section, including `phase` if two-phase mode is enabled
- `fix nve` and `fix langevin` applied to AT group only
- For two-phase mode: the run sequence SHALL include Phase 1 run, `fix_modify phase 2`, and Phase 2 run
- Thermo and dump configuration

#### Scenario: Polymer system with cross dihedrals
- **WHEN** the settings describe a polymer with cross bonds, cross angles, and cross dihedrals (both RB and tabulated)
- **THEN** the `.in` file SHALL contain `dihedral_style hybrid backmap/ryckaert backmap/table` and correct `dihedral_coeff` lines for each type

#### Scenario: Two-phase backmapping script
- **WHEN** `simulation.two_phase: true` is set
- **THEN** the `.in` file SHALL contain `fix backmap ... phase 1`, a Phase 1 run block, `fix_modify bm phase 2`, and a Phase 2 run block

#### Scenario: System with dihedrals but no tabulated dihedrals
- **WHEN** the settings define only RB cross dihedrals (no tabulated)
- **THEN** the `.in` file SHALL use `dihedral_style hybrid backmap/ryckaert` (without `backmap/table`)

### Requirement: Simulation parameters section

The `simulation` section SHALL configure backmapping parameters. The following parameters SHALL be added to support Phase 2 features:

```yaml
simulation:
  two_phase: false              # enable two-phase backmapping protocol
  alpha2: null                  # lambda increment for phase 2 (defaults to alpha)
  cap_force: null               # force capping value in kJ/(mol·nm)
  cap_force_ramp: null          # gradual force cap ramp (steps)
  em_steps: 0                   # energy minimization steps between phases
  em_ftol: 10.0                 # EM force tolerance
```

When `two_phase: true`, the generator SHALL produce an input script with Phase 1 and Phase 2 run blocks. If `alpha2` is specified, Phase 2 SHALL use a different lambda ramp rate.

#### Scenario: Two-phase parameters
- **WHEN** `two_phase: true` and `alpha2: 0.0002` are set
- **THEN** the Phase 1 run SHALL use `alpha` for the ramp rate and Phase 2 SHALL use `fix_modify bm alpha 0.0002` (or equivalent) for its ramp rate

#### Scenario: Force capping
- **WHEN** `cap_force: 1000.0` is set (kJ/(mol·nm))
- **THEN** the input script SHALL include `fix_modify bm cap_force 23.9006` (converted to kcal/(mol·Å))

### Requirement: Molecules section — CG-AT mapping

The `molecules` section SHALL define one or more CG molecule types and their atomistic-to-coarse-grained mapping.

**Phase 3 additions** — Reactive network support:

```yaml
molecules:
  - name: EPO
    ident: EPO
    source:
      coordinates:
        - file: epon-828.gro
          molecule_degree: 0
        - file: epon-828_deg1_A1.gro
          molecule_degree: 1
          when: A1
      topology:
        - file: epon-828.itp
          molecule_degree: 0
        - file: epon-828_deg1_A1.itp
          molecule_degree: 1
          when: A1
    beads:
      - name: A1
        type: A
        atoms_by_degree:
          - degree: 1
            molecule_degree: "0"
            atoms: [1:EPO:C1, 1:EPO:O1]
          - degree: 2
            molecule_degree: "1,2"
            atoms: [1:EPO:C1, 1:EPO:O1, 1:EPO:H25]
            active_site: "EPO:C1:4"
        remove:
          - active_site: "MOL:ATOM"
            atoms: [1:EPO:H8]
    charge_management:
      equilibrate: true
      transfers:
        - when: "IPD:N1:2"
          from_atom: "IPD:H8"
          to_atoms: "EPO:C1#H25,EPO:C21#H26"
    charge_map:
      A1: [0.1, -0.2]
    type_map:
      A1: [opls_135, opls_140]
```

The generator SHALL handle:
- **Degree-dependent source files**: Different AT coordinates and topologies for each bonding degree of the CG molecule
- **Degree-dependent bead atoms**: `atoms_by_degree` lists specifying different AT atom sets depending on the bead's bonding degree
- **Active sites**: AT atoms that can form new bonds, identified by `active_site: "molecule:atom:max_degree"`
- **Atom removal on bond formation**: `remove` entries specifying atoms to delete when a bond forms at an active site
- **Charge management**: `equilibrate` flag and `transfers` rules for redistributing charges upon bond formation
- **Charge and type maps**: Per-bead overrides for AT atom charges and types based on degree

#### Scenario: Degree-dependent bead with two degrees
- **WHEN** a bead defines `atoms_by_degree` with degree 1 (3 atoms) and degree 2 (4 atoms)
- **THEN** the generator SHALL create different atom type sets for each degree variant and assign CG beads to the correct variant based on the CG topology's connectivity

#### Scenario: Active site detection
- **WHEN** a bead has `active_site: "EPO:C1:4"` (atom C1 in molecule EPO, max degree 4)
- **THEN** the generator SHALL identify C1 as a potential bond-formation site and configure the appropriate extra bond types

#### Scenario: Multiple molecule types
- **WHEN** the YAML defines EPO, HDD, and IPD molecules with inter-molecule cross bonds
- **THEN** the generator SHALL handle cross-molecule bonds and assign unique LAMMPS types across all molecule types

### Requirement: Source file parsing

The generator SHALL parse AT and CG source files to extract topology and coordinates.

**Phase 4 additions** — Non-GROMACS format support:

The parser module SHALL support:
- **PDB files**: Atom positions, residue information, chain identifiers. Topology inferred from CONECT records or a separate topology file.
- **LAMMPS data files**: Direct reading of existing LAMMPS data files as CG input. No coordinate conversion needed.
- **XYZ files**: Atom positions only. Requires a separate topology source.

Format detection SHALL use the `format` field in the YAML settings or auto-detect from file extension:
- `.gro` → GROMACS coordinates
- `.top`, `.itp` → GROMACS topology
- `.pdb` → PDB coordinates
- `.data` → LAMMPS data file
- `.xyz` → XYZ coordinates

**LAMMPS `.table` passthrough**: When a table file has `.table` extension (native LAMMPS format), the generator SHALL use it directly without conversion.

#### Scenario: PDB coordinate input
- **WHEN** `cg_system.format: pdb` and `cg_system.coordinates: system.pdb`
- **THEN** the generator SHALL read atom positions from the PDB file and convert to LAMMPS units

#### Scenario: LAMMPS data file as CG input
- **WHEN** `cg_system.format: lammps` and `cg_system.coordinates: cg_system.data`
- **THEN** the generator SHALL read CG positions and topology directly from the LAMMPS data file without unit conversion

#### Scenario: Native LAMMPS table passthrough
- **WHEN** a cross-interaction table file has `.table` extension
- **THEN** the generator SHALL copy the file as-is without converting units or format

#### Scenario: Unsupported format
- **WHEN** the `format` field specifies an unrecognized format
- **THEN** the generator SHALL abort with: "Unsupported format 'X'. Supported formats: gromacs, pdb, lammps, xyz"

### Requirement: Feature phasing

The generator SHALL be implemented in phases. The YAML schema supports all features from Phase 1 through Phase 4, but deferred features SHALL produce a clear "not yet implemented" error if encountered in the settings file.

**Phase 2** — Extended bonded interactions:
- Cross dihedrals and 1-4 pairs in `cross_interactions`
- Dihedral section in LAMMPS data file
- `dihedral_style hybrid backmap/*` in LAMMPS input script
- Ryckaert-Bellemans coefficient conversion (GROMACS → LAMMPS convention)
- Dihedral table conversion (`.xvg` → `.table`)
- Two-phase backmapping (`two_phase: true`)
- Force capping (`cap_force`)

**Phase 3** — Reactive networks:
- Degree-dependent bead definitions (`atoms_by_degree`)
- Active sites and bond formation logic
- Charge management (equilibration, transfers, charge_map, type_map)
- Atom removal on bond formation (`remove`)
- Predefined active sites file
- Restricted cross-bond patterns
- Multiple molecule types in one system

**Phase 4** — Format flexibility:
- PDB, LAMMPS data, XYZ coordinate parsers
- Non-GROMACS CG input formats
- Native LAMMPS `.table` format as input (skip conversion)

#### Scenario: Phase 2 dihedral feature used
- **WHEN** the YAML contains `cross_interactions.dihedrals` entries
- **THEN** the generator SHALL process them and produce dihedral sections in the output files

#### Scenario: Phase 3 feature used before implementation
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) before Phase 3 is implemented
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 4 format before implementation
- **WHEN** `cg_system.format: pdb` is specified before Phase 4 is implemented
- **THEN** the generator SHALL abort with: "Format 'pdb' is not yet supported (planned for Phase 4)"
