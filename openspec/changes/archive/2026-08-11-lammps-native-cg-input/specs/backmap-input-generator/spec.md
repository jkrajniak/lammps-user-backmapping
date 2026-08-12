## MODIFIED Requirements

### Requirement: CG system section

The `cg_system` section SHALL specify the CG configuration files:

```yaml
cg_system:
  coordinates: cg_conf.gro
  topology: cg_topol.top
  format: gromacs                   # gromacs (default) or lammps
  # Deferred: predefined active sites for reactive networks
  # predefined_active_sites: active_sites.txt
```

or, for a LAMMPS-native CG system:

```yaml
cg_system:
  format: lammps
  data: cg_system.data              # LAMMPS data file: box, per-atom type/charge/mass
                                     # (Masses + Atoms # full sections)
```

The generator SHALL read CG coordinates and topology to determine molecule count, CG positions, and per-atom type/mass/charge, regardless of `format`.

For `format: lammps`, the generator SHALL read CG atom positions, box dimensions, and per-atom type/charge/mass directly from the `data` file's `Masses` and `Atoms # full` sections, matched positionally per molecule in the same way as the `format: gromacs` path (no atom-name matching is used on the CG side, and no `data`-file `Bonds`/`Angles` section is read — CG-CG bonded connectivity is configured via `cross_interactions`, identically for both formats; see "Bonded and nonbonded CG-CG terms are format-independent" below). The `data` file SHALL be assumed to already be in LAMMPS `real` units (e.g. as produced by `write_data` under `units real`); no unit conversion is applied to CG positions or box dimensions read this way. Each CG atom's LAMMPS numeric type ID (from the `Atoms` section) SHALL be used, as a string, as that atom's type label — `molecules[].beads[].type` in the YAML settings MUST reference this same numeric ID (as a string, e.g. `type: "1"`) when `cg_system.format: lammps`, since a `data` file carries no symbolic type name.

AT-fragment sources (`molecules[].source`) remain GROMACS-only; LAMMPS-native AT fragment input is not supported by this requirement (tracked separately, blocked on an atom-identity naming design, since AT-fragment bead-mapping is matched by free-text atom name which a `data` file cannot supply).

**Bonded and nonbonded CG-CG terms are format-independent.** Regardless of `cg_system.format`, CG-CG bonded interactions (bonds/angles between beads) are configured via `cross_interactions` in the YAML settings, and CG-CG nonbonded (pair) interactions are configured via `simulation.table_groups` plus `table_*` files present alongside the settings file. Neither is read from `cg_system.topology`'s or `cg_system.data`'s own bonded/pair sections in either format — this was already true for `format: gromacs` before this change (the CG `.top` file's own `[bonds]`/`[angles]` sections are parsed but not otherwise consumed).

#### Scenario: GROMACS CG system
- **WHEN** the CG system specifies `format: gromacs` with `.gro` and `.top` files
- **THEN** the generator SHALL read CG positions, box dimensions, and per-atom type/mass/charge from these files

#### Scenario: LAMMPS-native CG system
- **WHEN** `cg_system.format: lammps` and `cg_system.data: cg_system.data` is a valid LAMMPS data file with `Masses` and `Atoms # full` sections
- **THEN** the generator SHALL read CG positions, box dimensions, and per-atom type/mass/charge from the `data` file, without any unit conversion, and use each atom's numeric LAMMPS type ID (as a string) as its type label

#### Scenario: Bead type references a numeric LAMMPS type ID
- **WHEN** `cg_system.format: lammps` and `molecules[].beads[].type: "1"` matches type ID 1 in the `data` file's `Masses` section
- **THEN** the generator SHALL resolve that bead's CG atom type to the mass and atom count declared for type 1

#### Scenario: Missing required section in LAMMPS data file
- **WHEN** `cg_system.format: lammps` and the `data` file has no `Masses` section, or no `Atoms` section
- **THEN** the generator SHALL abort with an error naming the missing section

#### Scenario: CG-CG bonded terms unaffected by CG format
- **WHEN** `cg_system.format: lammps` and `cross_interactions.bonds`/`.angles` declare CG-CG cross terms (`cg_bonded: true`)
- **THEN** the generator SHALL process them identically to the `format: gromacs` case, since they are read entirely from `cross_interactions`, not from `cg_system`

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

The parser SHALL be modular so that additional source formats (PDB, LAMMPS data, XYZ) can be added later, and so that additional CG-side formats can be added independently of AT-fragment formats.

**LAMMPS-native CG system** — implemented (see "CG system section" requirement above): a LAMMPS `data`-file parser producing the same internal shape (`GroFile`/`Topology`/`AtomType`/`TopAtom` equivalents) `builder.py` already consumes from the GROMACS `cg_gro`/`cg_top` path, so downstream molecule-replication, bead-mapping, and hybrid-build logic requires no changes.

**Still deferred** — non-GROMACS AT-fragment sources, and PDB/XYZ CG or AT input:
- **PDB files**: Atom positions, residue information, chain identifiers. Topology inferred from CONECT records or a separate topology file.
- **LAMMPS data files as AT-fragment input**: blocked on an atom-identity naming design, since AT-fragment bead-mapping is matched by free-text atom name (`atoms: [C1, H1, ...]`), which a LAMMPS `data` file has no field for.
- **XYZ files**: Atom positions only. Requires a separate topology source.

Format detection SHALL use the `format` field in the YAML settings:
- `cg_system.format: gromacs` (default) → GROMACS `.gro`/`.top` CG system
- `cg_system.format: lammps` → LAMMPS `data` CG system
- `.pdb` → PDB coordinates (deferred)
- `.xyz` → XYZ coordinates (deferred)

**LAMMPS `.table` passthrough**: When a table file has `.table` extension (native LAMMPS format), the generator SHALL use it directly without conversion. The `table_groups` CG-CG pair-table lookup SHALL check for a `.table`-extension file (`table_<a>_<b>.table`) before falling back to the `.xvg`-extension naming convention (`table_<a>_<b>.xvg`, converted). This applies identically for both `cg_system.format` values, since `table_groups` is unrelated to CG-system format.

#### Scenario: GROMACS AT topology
- **WHEN** the source topology is a GROMACS `.top` file with `#include` directives
- **THEN** the generator SHALL resolve includes and parse the full topology tree

#### Scenario: GROMACS coordinates
- **WHEN** the source coordinates are a `.gro` file
- **THEN** the generator SHALL extract positions (nm) and box dimensions (nm)

#### Scenario: LAMMPS data file as CG input
- **WHEN** `cg_system.format: lammps` and `cg_system.data: cg_system.data`
- **THEN** the generator SHALL read CG positions and per-atom type/mass/charge directly from the LAMMPS data file without unit conversion

#### Scenario: Native LAMMPS table passthrough for CG pair tables
- **WHEN** `table_A_B.table` exists alongside the settings file and type names `A`/`B` (or, for a LAMMPS-native CG system, their numeric type-ID string equivalents) are listed in `simulation.table_groups`
- **THEN** the generator SHALL use `table_A_B.table` directly without converting units or format, in preference to `table_A_B.xvg`

#### Scenario: Unsupported format
- **WHEN** the `cg_system.format` field specifies an unrecognized format
- **THEN** the generator SHALL abort with: "Unsupported format 'X'. Supported formats: gromacs, lammps"

### Requirement: Unit conversion

All quantities sourced from GROMACS input (`format: gromacs`) SHALL be converted from GROMACS units to LAMMPS `real` units:

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

CG positions and box dimensions sourced from a LAMMPS-native CG system (`cg_system.format: lammps`) SHALL NOT be converted: the `data` file is assumed to already be in `real` units (see "CG system section" requirement), so values are used as-is. (Mass and charge already have an identity conversion factor for the GROMACS path too, so this distinction only matters for distance-valued quantities.)

#### Scenario: Position conversion
- **WHEN** a GROMACS coordinate is 3.75 nm
- **THEN** the LAMMPS coordinate SHALL be 37.5 Angstrom

#### Scenario: LJ parameter conversion
- **WHEN** GROMACS sigma = 0.3166 nm and epsilon = 0.6502 kJ/mol
- **THEN** LAMMPS sigma SHALL be 3.166 Angstrom and epsilon SHALL be 0.15539 kcal/mol

#### Scenario: Bond parameter conversion
- **WHEN** GROMACS harmonic bond has k = 345000 kJ/(mol*nm^2) and r0 = 0.1 nm
- **THEN** LAMMPS harmonic bond SHALL have k = 82.76 kcal/(mol*Angstrom^2) and r0 = 1.0 Angstrom

#### Scenario: No conversion for LAMMPS-native CG position
- **WHEN** `cg_system.format: lammps` and the `data` file's `Atoms` section gives a CG atom position of `x = 37.5`
- **THEN** the generator SHALL use x = 37.5 Angstrom unchanged (no nm→Angstrom conversion)

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
- Non-GROMACS AT-fragment source formats (LAMMPS, PDB, XYZ) — deferred; LAMMPS AT-fragment input additionally blocked on an atom-identity naming design (AT-fragment bead-mapping is name-keyed, unlike the CG side)
- PDB/XYZ CG input formats — deferred

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

#### Scenario: Phase 4 AT-fragment format before implementation
- **WHEN** a `molecules[].source` entry specifies a non-GROMACS format
- **THEN** the generator SHALL abort with: "Format 'X' is not yet supported for AT-fragment sources (planned for a future phase)"

#### Scenario: PE validation in Phase 1
- **WHEN** `backmap-prep` processes the PE example `settings.yaml` with 50 beads/chain, 2 atoms/bead, tabulated CG bonds and angles
- **THEN** the output SHALL be a valid LAMMPS data file and input script that runs without errors

#### Scenario: Melamine validation in Phase 1
- **WHEN** `backmap-prep` processes the melamine example `settings.yaml` with 3 beads/molecule in a triangular topology
- **THEN** the output SHALL be a valid LAMMPS data file with correct triangular bonding and input script
