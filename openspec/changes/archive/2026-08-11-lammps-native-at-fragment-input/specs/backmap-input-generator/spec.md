## MODIFIED Requirements

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

or, for a LAMMPS-native AT fragment:

```yaml
molecules:
  - name: DOD
    ident: DOD
    source:
      format: lammps
      data: dodecane_at.data          # LAMMPS data file: box, per-atom type/charge/mass,
                                       # Bonds/Angles/Dihedrals connectivity by type ID
      input_script: dodecane_at_ff.in # bond_coeff/angle_coeff/dihedral_coeff/pair_coeff
    beads:
      - name: A1
        type: A
        atoms:
          - 1:DOD:1                   # LAMMPS numeric atom ID, as a string
          - 1:DOD:2
```

Atom references use the bakery format: `chain_idx:molecule_name:atom_name`. When the owning molecule's `source.format: lammps`, `atom_name` SHALL be the LAMMPS numeric atom ID (as a string) from the `data` file's `Atoms` section, since a `data` file has no symbolic atom-name field — this mirrors the numeric-type-ID convention `cg_system.format: lammps` already uses for `beads[].type`.

The bead order in the YAML defines the CG bead ordering. All AT atoms listed in a bead's `atoms` array belong to that CG bead (intra-CG). AT atoms in different beads that are bonded to each other form cross-CG connections.

**Phase 1 (MVP):** Single molecule type, simple `atoms` list per bead, single source coordinate/topology file.

**LAMMPS-native AT fragment source** — implemented for single-file (non-degree-dependent) sources on the linear engine (see "Source file parsing" requirement below). Not supported for degree-dependent sources (`atoms_by_degree`, below) or the network engine.

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

#### Scenario: LAMMPS-native AT fragment
- **WHEN** `molecules[0].source.format: lammps` with a `data` file and `input_script`, and `beads[].atoms` entries reference numeric LAMMPS atom IDs
- **THEN** the generator SHALL resolve bead membership by atom ID exactly as the GROMACS path resolves it by atom name

#### Scenario: Multiple molecule types (deferred)
- **WHEN** the YAML defines EPO, HDD, and IPD molecules for a reactive network
- **THEN** the generator SHALL handle all molecule types and their inter-molecule cross bonds

#### Scenario: Degree-dependent bead with two degrees
- **WHEN** a bead defines `atoms_by_degree` with degree 1 (3 atoms) and degree 2 (4 atoms)
- **THEN** the generator SHALL create different atom type sets for each degree variant and assign CG beads to the correct variant based on the CG topology's connectivity

#### Scenario: Active site detection
- **WHEN** a bead has `active_site: "EPO:C1:4"` (atom C1 in molecule EPO, max degree 4)
- **THEN** the generator SHALL identify C1 as a potential bond-formation site and configure the appropriate extra bond types

#### Scenario: LAMMPS-native source rejected for degree-dependent mapping
- **WHEN** `molecules[0].source.format: lammps` and `source.coordinates`/`source.topology` are given as degree-dependent lists (Phase 3 shape) instead of `data`/`input_script`
- **THEN** the generator SHALL abort with an error naming the unsupported combination

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

The parser SHALL be modular so that additional source formats (PDB, LAMMPS data, XYZ) can be added later, and so that additional CG-side and AT-fragment-side formats can be added independently of each other.

**LAMMPS-native CG system** — implemented (see "CG system section" requirement): a LAMMPS `data`-file parser producing the same internal shape (`GroFile`/`Topology`/`AtomType`/`TopAtom` equivalents) `builder.py` already consumes from the GROMACS `cg_gro`/`cg_top` path.

**LAMMPS-native AT fragment** — implemented for single-file, linear-engine sources (see "Molecules section — CG-AT mapping" requirement): a LAMMPS `data`-file reader that, unlike the CG reader, DOES parse `Bonds`/`Angles`/`Dihedrals` sections (these are load-bearing for AT fragments — `at_mol.bonds`/`.angles`/`.dihedrals` drive the actual intra-bead bond/angle/dihedral coefficients, unlike the CG side's `cg_mol.bonds`/`.angles`, which are parsed but never consumed). Paired with a bounded input-script reader covering exactly `bond_coeff` (for `bond_style harmonic`), `angle_coeff` (`angle_style harmonic`), `dihedral_coeff` (`dihedral_style ryckaert`), and `pair_coeff i i` diagonal (self) entries — the only forms `builder.py` consumes for AT fragments. Cross-type (`i != j`) `pair_coeff` lines, if present, SHALL be ignored (cross-type LJ parameters are always computed in Python via mixing, never read from source). Any bond/angle func or dihedral func the GROMACS path itself would skip (non-harmonic bonds/angles, non-RB dihedrals) SHALL be skipped identically on the LAMMPS-native path — no new functional-form support is added by this requirement.

**Still deferred** — PDB/XYZ CG or AT input, and degree-dependent (Phase 3) AT-fragment sources of any format:
- **PDB files**: Atom positions, residue information, chain identifiers. Topology inferred from CONECT records or a separate topology file.
- **XYZ files**: Atom positions only. Requires a separate topology source.
- **GROMACS `[virtual_sites3]`-equivalent AT fragments**: not supported for `format: lammps` — no LAMMPS `data`-file construct maps onto a GROMACS virtual/dummy interaction site; such fragments must stay GROMACS-format.

Format detection SHALL use the `format` field in the YAML settings:
- `cg_system.format: gromacs` (default) → GROMACS `.gro`/`.top` CG system
- `cg_system.format: lammps` → LAMMPS `data` CG system
- `molecules[].source.format: gromacs` (default) → GROMACS `.gro`/`.top` AT fragment
- `molecules[].source.format: lammps` → LAMMPS `data`/`input_script` AT fragment
- `.pdb` → PDB coordinates (deferred)
- `.xyz` → XYZ coordinates (deferred)

**LAMMPS `.table` passthrough**: unchanged from the CG-side requirement; AT fragments in this codebase are always analytic (never tabulated), so no AT-fragment table passthrough applies.

#### Scenario: GROMACS AT topology
- **WHEN** the source topology is a GROMACS `.top` file with `#include` directives
- **THEN** the generator SHALL resolve includes and parse the full topology tree

#### Scenario: GROMACS coordinates
- **WHEN** the source coordinates are a `.gro` file
- **THEN** the generator SHALL extract positions (nm) and box dimensions (nm)

#### Scenario: LAMMPS data file as CG input
- **WHEN** `cg_system.format: lammps` and `cg_system.data: cg_system.data`
- **THEN** the generator SHALL read CG positions and per-atom type/mass/charge directly from the LAMMPS data file without unit conversion

#### Scenario: LAMMPS-native AT fragment, bonded coefficients
- **WHEN** `molecules[0].source.format: lammps` with a `data` file whose `Bonds` section references bond type 1, and the paired `input_script` declares `bond_style harmonic` with `bond_coeff 1 <k> <r0>`
- **THEN** the generator SHALL use k and r0 unchanged (no unit conversion) for that bond's intra-bead coefficient

#### Scenario: LAMMPS-native AT fragment, dihedral coefficients
- **WHEN** the `input_script` declares `dihedral_style ryckaert` with `dihedral_coeff <id> <C0> <C1> <C2> <C3> <C4> <C5>`
- **THEN** the generator SHALL use the six coefficients unchanged as the intra-bead Ryckaert–Bellemans dihedral parameters, with no GROMACS RB→LAMMPS conversion applied

#### Scenario: LAMMPS-native AT fragment, ignored cross pair-coeff
- **WHEN** the `input_script` includes a `pair_coeff i j` line with `i != j`
- **THEN** the generator SHALL ignore it and compute the corresponding cross-type LJ parameters via the existing arithmetic/geometric mixing, exactly as it does for the GROMACS path

#### Scenario: Missing coefficient for a referenced type
- **WHEN** the `data` file's `Bonds`/`Angles`/`Dihedrals` section references a type ID with no corresponding `bond_coeff`/`angle_coeff`/`dihedral_coeff` in the `input_script`
- **THEN** the generator SHALL abort with an error naming the missing type ID

#### Scenario: Native LAMMPS table passthrough for CG pair tables
- **WHEN** `table_A_B.table` exists alongside the settings file and type names `A`/`B` (or, for a LAMMPS-native CG system, their numeric type-ID string equivalents) are listed in `simulation.table_groups`
- **THEN** the generator SHALL use `table_A_B.table` directly without converting units or format, in preference to `table_A_B.xvg`

#### Scenario: Unsupported format
- **WHEN** the `cg_system.format` or `molecules[].source.format` field specifies an unrecognized format
- **THEN** the generator SHALL abort with: "Unsupported format 'X'. Supported formats: gromacs, lammps"

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
- LAMMPS-native CG system (`data` file: box, `Masses`, `Atoms # full`; `.table` pair passthrough) — **implemented**
- LAMMPS-native AT fragment (`data` file + bounded `input_script`; numeric atom-ID bead references), single-file/linear-engine sources only — **implemented**
- PDB/XYZ CG or AT input formats — deferred
- Degree-dependent (Phase 3) AT-fragment sources in LAMMPS format — deferred
- GROMACS `[virtual_sites3]`-equivalent AT fragments in LAMMPS format — deferred (no direct LAMMPS construct)

#### Scenario: Deferred feature used
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) in Phase 1
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 2 dihedral feature used
- **WHEN** the YAML contains `cross_interactions.dihedrals` entries
- **THEN** the generator SHALL process them and produce dihedral sections in the output files

#### Scenario: Phase 3 feature used before implementation
- **WHEN** the YAML contains `atoms_by_degree` (Phase 3 feature) before Phase 3 is implemented
- **THEN** the generator SHALL abort with: "Feature 'atoms_by_degree' is not yet implemented (planned for Phase 3)"

#### Scenario: Phase 4 CG format available
- **WHEN** `cg_system.format: lammps` is specified
- **THEN** the generator SHALL process it as described in the "CG system section" requirement (no longer a deferred-feature error)

#### Scenario: Phase 4 AT-fragment format available
- **WHEN** `molecules[0].source.format: lammps` is specified with `data`/`input_script` (single-file, linear engine)
- **THEN** the generator SHALL process it as described in the "Molecules section — CG-AT mapping" and "Source file parsing" requirements (no longer a deferred-feature error)

#### Scenario: Phase 4 AT-fragment format still deferred for network/degree-dependent use
- **WHEN** `molecules[0].source.format: lammps` is combined with degree-dependent source files, or `prep.engine: network`
- **THEN** the generator SHALL abort with: "cg_system/source format 'lammps' is not yet supported for degree-dependent or network-engine sources (planned for a future phase)"

#### Scenario: PE validation in Phase 1
- **WHEN** `backmap-prep` processes the PE example `settings.yaml` with 50 beads/chain, 2 atoms/bead, tabulated CG bonds and angles
- **THEN** the output SHALL be a valid LAMMPS data file and input script that runs without errors

#### Scenario: Melamine validation in Phase 1
- **WHEN** `backmap-prep` processes the melamine example `settings.yaml` with 3 beads/molecule in a triangular topology
- **THEN** the output SHALL be a valid LAMMPS data file with correct triangular bonding and input script
