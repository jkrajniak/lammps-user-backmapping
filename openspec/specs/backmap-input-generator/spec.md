## ADDED Requirements

### Requirement: YAML settings file as primary input

The generator SHALL read a YAML settings file as its primary input. This file describes the complete system configuration for backmapping: CG-AT mapping, source files, cross-CG interactions, simulation parameters, and output options. The format is a modernized version of the bakery XML settings format (VOTCA-inspired).

The YAML file has five top-level sections:

```yaml
molecules:          # CG molecule definitions and AT mapping
cg_system:          # CG configuration files
cross_interactions: # Cross-CG bonds/angles/dihedrals
simulation:         # Backmapping simulation parameters
output:             # Output file configuration
```

The generator SHALL validate the YAML against a Pydantic schema and report clear errors for missing or invalid fields.

#### Scenario: Valid settings file
- **WHEN** the generator reads a valid `settings.yaml` with all required sections
- **THEN** it SHALL parse without errors and proceed to file generation

#### Scenario: Missing required field
- **WHEN** the YAML is missing a required field (e.g., `molecules[0].name`)
- **THEN** the generator SHALL abort with a validation error naming the missing field and its expected type

#### Scenario: Unknown fields
- **WHEN** the YAML contains unknown fields
- **THEN** the generator SHALL issue a warning but continue (forward compatibility)

### Requirement: Molecules section — CG-AT mapping

The `molecules` section SHALL define one or more CG molecule types and their atomistic-to-coarse-grained mapping. Each molecule has:

```yaml
molecules:
  - name: DOD                       # CG molecule name
    ident: DOD                      # Name in AT topology (may differ from name)
    source:                         # AT source files
      coordinates: dodecane_single.gro
      topology: topol_aa.top
    beads:                          # CG bead definitions
      - name: A1
        type: A
        atoms:                      # AT atoms mapped to this bead
          - 1:DOD:C1
          - 1:DOD:C2
      - name: B1
        type: B
        atoms:
          - 1:DOD:C3
          - 1:DOD:C4
```

Atom references use the bakery format: `chain_idx:molecule_name:atom_name`.

The bead order in the YAML defines the CG bead ordering. All AT atoms listed in a bead's `atoms` array belong to that CG bead (intra-CG). AT atoms in different beads that are bonded to each other form cross-CG connections.

**Phase 1 (MVP):** Single molecule type, simple `atoms` list per bead, single source coordinate/topology file.

**Deferred features** (spec'd for completeness, implemented in Phase 3):

```yaml
molecules:
  - name: EPO
    ident: EPO
    source:
      # Degree-dependent source files for reactive networks
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
        # Degree-dependent atom lists
        atoms_by_degree:
          - degree: 1
            molecule_degree: "0"
            atoms: [1:EPO:C1, 1:EPO:O1, ...]
          - degree: 2
            molecule_degree: "1,2"
            atoms: [1:EPO:C1, 1:EPO:O1, ..., 1:EPO:H25]
            active_site: "EPO:C1:4"
        # Atom removal on bond formation
        remove:
          - active_site: "MOL:ATOM"
            atoms: [1:EPO:H8]
    # Charge management
    charge_management:
      equilibrate: true
      transfers:
        - when: "IPD:N1:2"
          from_atom: "IPD:H8"
          to_atoms: "EPO:C1#H25,EPO:C21#H26"
    # Per-bead charge/type overrides
    charge_map:
      A1: [0.1, -0.2, ...]    # one value per atom, or "*" for original
    type_map:
      A1: [opls_135, ...]     # one value per atom, or "*" for original
```

The generator SHALL handle (Phase 3):
- **Degree-dependent source files**: Different AT coordinates and topologies for each bonding degree of the CG molecule
- **Degree-dependent bead atoms**: `atoms_by_degree` lists specifying different AT atom sets depending on the bead's bonding degree
- **Active sites**: AT atoms that can form new bonds, identified by `active_site: "molecule:atom:max_degree"`
- **Atom removal on bond formation**: `remove` entries specifying atoms to delete when a bond forms at an active site
- **Charge management**: `equilibrate` flag and `transfers` rules for redistributing charges upon bond formation
- **Charge and type maps**: Per-bead overrides for AT atom charges and types based on degree

#### Scenario: Dodecane molecule definition
- **WHEN** the YAML defines a DOD molecule with 6 beads (A1, B1..B4, A2), each mapping 2 carbon atoms
- **THEN** the generator SHALL create 6 CG beads + 12 AT atoms per molecule, with correct CG-AT mapping

#### Scenario: Multiple molecule types (deferred)
- **WHEN** the YAML defines EPO, HDD, and IPD molecules for a reactive network
- **THEN** the generator SHALL handle all molecule types and their inter-molecule cross bonds

#### Scenario: Degree-dependent bead with two degrees
- **WHEN** a bead defines `atoms_by_degree` with degree 1 (3 atoms) and degree 2 (4 atoms)
- **THEN** the generator SHALL create different atom type sets for each degree variant and assign CG beads to the correct variant based on the CG topology's connectivity

#### Scenario: Active site detection
- **WHEN** a bead has `active_site: "EPO:C1:4"` (atom C1 in molecule EPO, max degree 4)
- **THEN** the generator SHALL identify C1 as a potential bond-formation site and configure the appropriate extra bond types

### Requirement: CG system section

The `cg_system` section SHALL specify the CG configuration files:

```yaml
cg_system:
  coordinates: cg_conf.gro
  topology: cg_topol.top
  format: gromacs                   # currently only gromacs supported
  # Deferred: predefined active sites for reactive networks
  # predefined_active_sites: active_sites.txt
```

The generator SHALL read CG coordinates and topology to determine molecule count, CG positions, and CG bonded topology.

#### Scenario: GROMACS CG system
- **WHEN** the CG system specifies `format: gromacs` with `.gro` and `.top` files
- **THEN** the generator SHALL read CG positions, box dimensions, and CG bond topology from these files

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

### Requirement: Simulation parameters section

The `simulation` section SHALL configure backmapping parameters:

```yaml
simulation:
  # Dynamic resolution
  alpha: 0.001                       # lambda increment per step
  initial_resolution: 0.0            # initial lambda value
  nonuniform_lambda: false           # staggered initial lambda

  # Time stepping
  timestep: 0.001                    # ps (converted to fs for LAMMPS)
  timestep_backmapping: 0.001        # ps, during backmapping phase

  # Run lengths
  equilibration_steps: 10000         # CG equilibration
  production_steps: 10000            # AT production after backmapping

  # Temperature and thermostat
  temperature: 423.0                 # K
  thermostat: langevin               # langevin | velocity_rescaling
  thermostat_gamma: 0.5              # coupling constant (1/ps)
  thermostat_target: atomistic       # atomistic | all | cg_only

  # Cutoffs
  lj_cutoff: 1.2                     # nm (converted to Angstrom)
  cg_cutoff: 1.4                     # nm
  coulomb_cutoff: 0.9                # nm

  # CG table groups
  table_groups: [WCG]                # CG atom types with tabulated potentials

  # Exclusions
  exclusion_nrexcl: 3                # nrexcl for special_bonds

  # Output
  energy_interval: 1000              # thermo output every N steps
  trajectory_interval: 1000          # dump every N steps

  # Random seed
  rng_seed: -1                       # -1 for random
```

**Phase 2 parameters** (two-phase backmapping, force capping):

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

**Other deferred parameters** (spec'd, implemented later):

```yaml
simulation:
  second_phase_em: false             # energy minimization in phase 2
  em_gamma: 0.0001                   # EM gamma parameter
  disable_angles: false              # disable angle interactions
  disable_dihedrals: false           # disable dihedral interactions
  coulomb_epsilon1: 1.0              # reaction field epsilon1
  coulomb_epsilon2: 78.0             # reaction field epsilon2
```

All parameters in GROMACS units. The generator converts to LAMMPS `real` units in the output.

#### Scenario: Default parameters
- **WHEN** the `simulation` section omits optional fields
- **THEN** the generator SHALL use the documented default values

#### Scenario: Parameter validation
- **WHEN** `alpha` is negative or `temperature` is zero
- **THEN** the generator SHALL abort with a validation error

#### Scenario: Two-phase parameters
- **WHEN** `two_phase: true` and `alpha2: 0.0002` are set
- **THEN** the Phase 1 run SHALL use `alpha` for the ramp rate and Phase 2 SHALL use `fix_modify bm alpha 0.0002` (or equivalent) for its ramp rate

#### Scenario: Force capping
- **WHEN** `cap_force: 1000.0` is set (kJ/(mol·nm))
- **THEN** the input script SHALL include `fix_modify bm cap_force 23.9006` (converted to kcal/(mol·Å))

### Requirement: Output section

The `output` section SHALL configure generator output:

```yaml
output:
  prefix: system                     # output file prefix
  format: lammps                     # target MD engine (only lammps for now)
  units: real                        # LAMMPS unit system
```

#### Scenario: Default output
- **WHEN** the `output` section is omitted
- **THEN** the generator SHALL use prefix `system`, format `lammps`, units `real`

### Requirement: Source file parsing

The generator SHALL parse AT and CG source files to extract topology and coordinates. In Phase 1 (MVP), the supported format is GROMACS:

- **`.gro` files**: atom positions and box dimensions
- **`.top`/`.itp` files**: atom types, charges, masses, bond/angle/dihedral topology, `#include` directives, combination rules

The parsing logic SHALL handle:
- `[ atomtypes ]` with `particletype` field (`A` = atomistic, `V` = virtual/CG)
- `[ atoms ]` with charge, mass, and type assignment
- `[ bonds ]`, `[ angles ]`, `[ dihedrals ]` with interaction parameters
- `[ virtual_sites3 ]` definitions
- `#include` directives and `#define` preprocessor directives
- `[ molecules ]` section for replication count
- Combination rules for LJ parameters (rules 1, 2, 3)

The parser SHALL be modular so that additional source formats (PDB, LAMMPS data, XYZ) can be added in Phase 4.

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

#### Scenario: GROMACS AT topology
- **WHEN** the source topology is a GROMACS `.top` file with `#include` directives
- **THEN** the generator SHALL resolve includes and parse the full topology tree

#### Scenario: GROMACS coordinates
- **WHEN** the source coordinates are a `.gro` file
- **THEN** the generator SHALL extract positions (nm) and box dimensions (nm)

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

Bond type assignment SHALL follow the three-category scheme from the design:
- **Intra-CG bonds** (from AT source `[ bonds ]`): assigned to static bond types
- **Cross-CG AT bonds** (from `cross_interactions.bonds` without `cg_bonded`): assigned to `backmap/*` AT bond types
- **Cross-CG CG bonds** (from `cross_interactions.bonds` with `cg_bonded: true`): assigned to `backmap/*` CG bond types

The same three-category distinction applies to angles, dihedrals, and 1-4 pairs.

Particle ordering within each molecule SHALL be: CG atom first, then AT atoms.

#### Scenario: Dodecane data file
- **WHEN** generating data for 100 dodecane molecules (6 beads, 12 AT atoms each)
- **THEN** each molecule SHALL contain 18 atoms (6 CG + 12 AT) with the same molecule_id, CG atoms listed first

#### Scenario: Dodecane data file with dihedrals
- **WHEN** generating data for dodecane molecules with cross dihedrals defined
- **THEN** the Dihedrals section SHALL contain entries for each cross-CG dihedral per molecule, with correct atom IDs and dihedral type

#### Scenario: Data file header includes dihedral counts
- **WHEN** the system has N_d total dihedrals and T_d dihedral types
- **THEN** the data file header SHALL include `N_d dihedrals` and `T_d dihedral types`

#### Scenario: Bond type mapping output
- **WHEN** the data file is generated
- **THEN** the generator SHALL print mapping tables:
  - Atom type mapping (e.g., "Type 1 = A, Type 2 = B, Type 3 = C1, ...")
  - Bond type mapping with category (e.g., "Bond type 1 = harmonic (intra-CG), Bond type 2 = backmap/harmonic at (cross-CG AT), ...")

### Requirement: Unit conversion

All quantities SHALL be converted from GROMACS units to LAMMPS `real` units:

| Quantity     | GROMACS          | LAMMPS real       | Factor     |
|-------------|------------------|-------------------|------------|
| Distance    | nm               | Angstrom          | x10        |
| Energy      | kJ/mol           | kcal/mol          | x0.239006  |
| Force       | kJ/(mol*nm)      | kcal/(mol*Angstrom) | x0.0239006 |
| Time        | ps               | fs                | x1000      |
| Charge      | e                | e                 | x1         |
| Mass        | g/mol            | g/mol             | x1         |
| Spring const (bond) | kJ/(mol*nm^2) | kcal/(mol*Angstrom^2) | x0.00239006 |
| Spring const (angle) | kJ/(mol*rad^2) | kcal/(mol*rad^2) | x0.239006 |

#### Scenario: Position conversion
- **WHEN** a GROMACS coordinate is 3.75 nm
- **THEN** the LAMMPS coordinate SHALL be 37.5 Angstrom

#### Scenario: LJ parameter conversion
- **WHEN** GROMACS sigma = 0.3166 nm and epsilon = 0.6502 kJ/mol
- **THEN** LAMMPS sigma SHALL be 3.166 Angstrom and epsilon SHALL be 0.15539 kcal/mol

#### Scenario: Bond parameter conversion
- **WHEN** GROMACS harmonic bond has k = 345000 kJ/(mol*nm^2) and r0 = 0.1 nm
- **THEN** LAMMPS harmonic bond SHALL have k = 82.76 kcal/(mol*Angstrom^2) and r0 = 1.0 Angstrom

### Requirement: Tabulated potential conversion

The generator SHALL convert tabulated potential files to LAMMPS table format. This applies to:
- **Non-bonded CG tables** (e.g., `table_WCG_WCG.xvg`) -> LAMMPS `pair_style table` format
- **Bonded CG tables** (e.g., `table_b1.xvg`) -> LAMMPS `bond_style backmap/table` format
- **Dihedral CG tables** (e.g., `table_d1.xvg`) -> LAMMPS `dihedral_style backmap/table` format

Input table formats supported in Phase 1: GROMACS `.xvg` (columns: r, V, F).
Phase 4 adds native LAMMPS `.table` format as input (no conversion needed).

Table file paths are specified in `cross_interactions` entries with `table:` field, or auto-detected from CG topology based on `table_groups`.

#### Scenario: CG pair table conversion
- **WHEN** the generator converts `table_WCG_WCG.xvg`
- **THEN** it SHALL produce `table_WCG_WCG.table` with distances in Angstrom, energies in kcal/mol, forces in kcal/(mol*Angstrom)

#### Scenario: CG bond table conversion
- **WHEN** the generator converts `table_b1.xvg` (a cross-CG CG bond)
- **THEN** it SHALL produce `table_b1.table` in LAMMPS bond table format

#### Scenario: LAMMPS table passthrough (Phase 4)
- **WHEN** a table file already has `.table` extension
- **THEN** the generator SHALL use it directly without conversion

### Requirement: LAMMPS input script generation

The generator SHALL produce a LAMMPS input script (`.in`) configured for backmapping. The script SHALL include:
- `units real` and `atom_style full`
- `read_data` for the generated data file
- `pair_style backmap` with AT and CG sub-styles, configured from `simulation.lj_cutoff` and `simulation.cg_cutoff`
- `pair_coeff` for all type pairs
- Bond/angle/dihedral styles:
  - If no cross interactions: static styles only (e.g., `bond_style harmonic`)
  - If cross interactions exist: `bond_style hybrid ...`, `angle_style hybrid ...`, `dihedral_style hybrid ...`
- Dihedral style routing: `dihedral_style hybrid backmap/ryckaert backmap/table` (only include sub-styles that are actually used)
- `bond_coeff` / `angle_coeff` / `dihedral_coeff` for each type, routing to correct sub-style based on category (intra-CG static vs cross-CG AT vs cross-CG CG)
- Group definitions for AT and CG atoms
- `fix backmap` with parameters from `simulation` section, including `phase` if two-phase mode is enabled
- `fix nve` and `fix langevin` applied to AT group only (or as configured by `thermostat_target`)
- For two-phase mode: the run sequence SHALL include Phase 1 run, `fix_modify phase 2`, and Phase 2 run
- `special_bonds` from `simulation.exclusion_nrexcl`
- Thermo output every `simulation.energy_interval` steps
- Dump configuration every `simulation.trajectory_interval` steps
- Three-phase run sequence: CG equilibration -> backmapping -> AT production

#### Scenario: Water2 system (no cross bonds)
- **WHEN** the settings describe a single-bead water system with no cross interactions
- **THEN** the `.in` file SHALL use static `bond_style harmonic` and `angle_style harmonic`, `pair_style backmap` for non-bonded, and `fix backmap` for lambda ramp

#### Scenario: Polymer system (with cross bonds)
- **WHEN** the settings describe a polymer with cross bonds, cross angles, and cross dihedrals
- **THEN** the `.in` file SHALL use `bond_style hybrid harmonic backmap/harmonic backmap/table`, routing each bond type to the correct sub-style

#### Scenario: Polymer system with cross dihedrals
- **WHEN** the settings describe a polymer with cross bonds, cross angles, and cross dihedrals (both RB and tabulated)
- **THEN** the `.in` file SHALL contain `dihedral_style hybrid backmap/ryckaert backmap/table` and correct `dihedral_coeff` lines for each type

#### Scenario: Two-phase backmapping script
- **WHEN** `simulation.two_phase: true` is set
- **THEN** the `.in` file SHALL contain `fix backmap ... phase 1`, a Phase 1 run block, `fix_modify bm phase 2`, and a Phase 2 run block

#### Scenario: System with dihedrals but no tabulated dihedrals
- **WHEN** the settings define only RB cross dihedrals (no tabulated)
- **THEN** the `.in` file SHALL use `dihedral_style hybrid backmap/ryckaert` (without `backmap/table`)

#### Scenario: Simulation parameters applied
- **WHEN** settings specify `timestep: 0.001` (ps), `temperature: 423.0`, `alpha: 0.0005`
- **THEN** the `.in` file SHALL contain `timestep 1.0` (fs), `fix langevin ... 423.0 423.0 ...`, `fix backmap ... alpha 0.0005`

### Requirement: Command-line interface

The generator SHALL be invoked as:
```
backmap-prep settings.yaml [--output-prefix PREFIX]
```

The settings YAML file is the single required argument. All configuration is in the YAML file. The optional `--output-prefix` overrides `output.prefix` from the YAML.

For backward compatibility with bakery workflows, the generator SHALL also accept:
```
backmap-prep --settings settings.yaml
```

#### Scenario: Minimal invocation
- **WHEN** the user runs `backmap-prep settings.yaml`
- **THEN** the generator SHALL produce `system.data`, `in.backmap`, and any converted table files (using the prefix from `output.prefix`)

#### Scenario: Custom prefix
- **WHEN** the user runs `backmap-prep settings.yaml --output-prefix my_system`
- **THEN** output files SHALL be named `my_system.data` and `in.my_system`

### Requirement: Python 3 compatibility

All Python code SHALL be Python 3.10+ compatible. The generator SHALL use:
- `PyYAML` for YAML parsing
- `Pydantic v2` for settings validation
- Standard library for file I/O

No dependency on Python 2 code or AdResSLab modules. Source file parsers (GROMACS `.gro`/`.top`) SHALL be new, clean Python 3 implementations.

#### Scenario: Python 3 execution
- **WHEN** the generator is run with Python 3.10+
- **THEN** it SHALL execute without errors

### Requirement: Exclusion list handling

The generator SHALL handle non-bonded exclusions based on `simulation.exclusion_nrexcl`:
- `nrexcl = 1`: exclude 1-2 pairs
- `nrexcl = 2`: exclude 1-2 and 1-3 pairs
- `nrexcl = 3`: exclude 1-2, 1-3, and 1-4 pairs (most common)

Exclusions SHALL be written as `special_bonds` configuration in the LAMMPS input script.

#### Scenario: nrexcl = 3
- **WHEN** `exclusion_nrexcl: 3`
- **THEN** the input script SHALL contain `special_bonds lj 0.0 0.0 0.0 coul 0.0 0.0 0.0` (exclude all 1-2, 1-3, 1-4 interactions)

#### Scenario: nrexcl = 2 with 1-4 scaling
- **WHEN** `exclusion_nrexcl: 2` and the topology defines fudgeLJ/fudgeQQ factors
- **THEN** the input script SHALL contain `special_bonds lj 0.0 0.0 FUDGE_LJ coul 0.0 0.0 FUDGE_QQ`

### Requirement: Feature phasing

The generator SHALL be implemented in phases. The YAML schema supports all features from Phase 1 through Phase 4, but deferred features SHALL produce a clear "not yet implemented" error if encountered in the settings file.

**Phase 1 (MVP)** — Linear molecules:
- YAML parsing with all five sections
- Single molecule type, simple `atoms` list per bead
- Cross bonds and angles from `cross_interactions`
- GROMACS source files (`.gro`/`.top`/`.itp`)
- Tabulated potential conversion (`.xvg` to `.table`)
- Unit conversion (GROMACS to LAMMPS real)
- Full LAMMPS `.data` and `.in` generation with `backmap/*` styles
- Validation: dodecane, PE (polyethylene), PE4, PE-10, PE-AA, melamine

**Phase 2** — Extended bonded interactions:
- Cross dihedrals and 1-4 pairs
- Tabulated bond/dihedral potentials
- Two-phase backmapping (`two_phase: true`)
- Energy minimization
- Force capping

**Phase 3** — Reactive networks:
- Degree-dependent bead definitions (`atoms_by_degree`)
- Active sites and bond formation logic
- Charge management (equilibration, transfers, charge_map, type_map)
- Atom removal on bond formation (`remove`)
- Predefined active sites file
- Restricted cross-bond patterns
- Multiple molecule types in one system
- Validation: epoxy networks (RIM135)

**Phase 4** — Format flexibility:
- Non-GROMACS source formats (PDB, LAMMPS data, XYZ)
- Non-GROMACS CG input formats
- Native LAMMPS `.table` format as input (skip conversion)

#### Scenario: Deferred feature used
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) in Phase 1
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 2 dihedral feature used
- **WHEN** the YAML contains `cross_interactions.dihedrals` entries
- **THEN** the generator SHALL process them and produce dihedral sections in the output files

#### Scenario: Phase 3 feature used before implementation
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) before Phase 3 is implemented
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 4 format before implementation
- **WHEN** `cg_system.format: pdb` is specified before Phase 4 is implemented
- **THEN** the generator SHALL abort with: "Format 'pdb' is not yet supported (planned for Phase 4)"

#### Scenario: PE validation in Phase 1
- **WHEN** `backmap-prep` processes the PE example `settings.yaml` with 50 beads/chain, 2 atoms/bead, tabulated CG bonds and angles
- **THEN** the output SHALL be a valid LAMMPS data file and input script that runs without errors

#### Scenario: Melamine validation in Phase 1
- **WHEN** `backmap-prep` processes the melamine example `settings.yaml` with 3 beads/molecule in a triangular topology
- **THEN** the output SHALL be a valid LAMMPS data file with correct triangular bonding and input script
