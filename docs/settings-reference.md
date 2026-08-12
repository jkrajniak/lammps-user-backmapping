# Settings Reference

The `backmap-prep` tool reads a YAML settings file that describes the CG-to-AT
mapping, cross interactions, simulation parameters, and output configuration.
This page documents every field.

## Top-Level Structure

```yaml
molecules:       # List of CG molecule definitions with AT mapping
cg_system:       # CG system coordinate and topology files
cross_interactions:  # Cross-CG bonded interactions
simulation:      # Backmapping simulation parameters
output:          # Output file configuration
```

---

## `molecules`

A list of molecule definitions. Each molecule describes how CG beads map to AT
atoms.

**Type**: `list[MoleculeDef]`

### `molecules[].name`

Name of the molecule type (used in cross-interaction references).

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `DOD` |

### `molecules[].ident`

Optional identifier for the molecule. Defaults to `name` if omitted. Used to
match residue names in topology files.

| | |
|---|---|
| **Type** | `string` or `null` |
| **Default** | `null` (uses `name`) |
| **Example** | `DOD` |

### `molecules[].source`

Paths to the atomistic coordinate and topology files for this molecule type.

#### `molecules[].source.coordinates`

Path to the atomistic coordinate file (GROMACS `.gro` format). Relative paths
are resolved from the settings file directory.

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `dodecane_single.gro` |

#### `molecules[].source.topology`

Path to the atomistic topology file (GROMACS `.top` format).

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `topol_aa.top` |

### `molecules[].beads`

List of CG bead definitions within this molecule.

**Type**: `list[BeadDef]`

#### `molecules[].beads[].name`

Name of this CG bead within the molecule. Must be unique within the molecule.
Used in cross-interaction references as `MOLECULE:BEAD` (e.g., `DOD:A1`).

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `A1` |

#### `molecules[].beads[].type`

CG bead type identifier. Beads with the same type share non-bonded interaction
parameters. When `cg_system.format: lammps`, this MUST be the LAMMPS numeric
atom type ID from the CG `data` file's `Masses` section, as a string (e.g.
`"1"`) — a `data` file has no symbolic type name, so the numeric ID is the
identifier. `simulation.table_groups` must reference the same strings.

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `A` (GROMACS CG system) or `"1"` (LAMMPS-native CG system) |

#### `molecules[].beads[].atoms`

List of AT atoms belonging to this bead. Each entry is a string in
`RESID:RESNAME:ATOMNAME` format matching the topology file.

| | |
|---|---|
| **Type** | `list[string]` |
| **Required** | yes |
| **Example** | `["1:DOD:C1", "1:DOD:C2"]` |

### Example

```yaml
molecules:
  - name: DOD
    ident: DOD
    source:
      coordinates: dodecane_single.gro
      topology: topol_aa.top
    beads:
      - name: A1
        type: A
        atoms:
          - 1:DOD:C1
          - 1:DOD:C2
      - name: B1
        type: B
        atoms:
          - 1:DOD:C3
          - 1:DOD:C4
```

---

## `cg_system`

Specifies the coarse-grained system files. Two formats are supported,
selected by `cg_system.format`: `gromacs` (default) reads coordinates and
topology from GROMACS `.gro`/`.top`; `lammps` reads a self-contained LAMMPS
`data` file instead.

CG-CG bonded interactions (bonds/angles between beads) and nonbonded
tabulated interactions are configured via `cross_interactions` and
`simulation.table_groups` respectively, identically for both `cg_system`
formats — neither is read from `cg_system` itself in either format.

### `cg_system.coordinates`

Path to the CG coordinate file (GROMACS `.gro` format). Required when
`format: gromacs`.

| | |
|---|---|
| **Type** | `string` |
| **Required** | when `format: gromacs` |
| **Example** | `cg_conf.gro` |

### `cg_system.topology`

Path to the CG topology file (GROMACS `.top` format). Required when
`format: gromacs`.

| | |
|---|---|
| **Type** | `string` |
| **Required** | when `format: gromacs` |
| **Example** | `topol_cg.top` |

### `cg_system.data`

Path to a LAMMPS `data` file supplying the CG system: box bounds, a
`Masses` section, and an `Atoms # full` section (atom ID, molecule ID,
type ID, charge, x, y, z). Required when `format: lammps`.

Values are read as-is, with **no unit conversion** — the file is assumed to
already be in LAMMPS `real` units, e.g. as written by `write_data` under
`units real`. Any `Bonds`/`Angles`/other sections in the file are tolerated
but not read (see the note above `cg_system.coordinates`).

Atoms must be ordered in contiguous per-molecule blocks (all of molecule
1's atoms, then all of molecule 2's, ...) — the same ordering assumption
`cg_system.coordinates` (`.gro`) already relies on.

| | |
|---|---|
| **Type** | `string` |
| **Required** | when `format: lammps` |
| **Example** | `cg_system.data` |

### `cg_system.format`

Input file format for the CG system.

| | |
|---|---|
| **Type** | `"gromacs"` \| `"lammps"` |
| **Default** | `"gromacs"` |

### Example

```yaml
cg_system:
  coordinates: cg_conf.gro
  topology: topol_cg.top
  format: gromacs
```

or, for a LAMMPS-native CG system:

```yaml
cg_system:
  format: lammps
  data: cg_system.data
```

See `examples/dodecane-lammps-cg/` for a complete worked example (a
LAMMPS-native CG system paired with the same GROMACS AT fragment and
`cross_interactions` as `examples/dodecane/`, verified to build to a
numerically identical hybrid system).

---

## `cross_interactions`

Defines bonded interactions that span CG bead boundaries. These are weighted
by lambda during the backmapping simulation.

### `cross_interactions.bonds`

List of cross-CG bond definitions.

**Type**: `list[CrossBond]`

#### `cross_interactions.bonds[].params`

Bond parameters as a space-separated string. The interpretation depends on
the bond style:

- For **harmonic** (`cg_bonded: false`): `"style_id K r0"` where `K` is the
  force constant and `r0` is the equilibrium distance
- For **tabulated** (`cg_bonded: true`): `"N style_id scale"` where `N` is
  the number of table entries

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `"1 0.153 334720.0"` (harmonic) or `"8 1 1.0"` (tabulated) |

#### `cross_interactions.bonds[].pairs`

List of atom pairs defining which bonds to create. Each pair is a two-element
list of `MOLECULE:BEAD` or `MOLECULE:ATOM` identifiers.

For CG cross bonds, use bead names: `[DOD:A1, DOD:B1]`.
For AT cross bonds, use atom names: `[DOD:C2, DOD:C3]`.

| | |
|---|---|
| **Type** | `list[list[string]]` |
| **Required** | yes |

#### `cross_interactions.bonds[].table`

Path to the tabulated potential file (for CG tabulated bonds). Omit for
harmonic bonds.

| | |
|---|---|
| **Type** | `string` or `null` |
| **Default** | `null` |
| **Example** | `table_b1.xvg` |

#### `cross_interactions.bonds[].cg_bonded`

Whether this bond is a CG-CG interaction (weighted by 1 - &lambda;<sub>i</sub>&lambda;<sub>j</sub>)
or an AT-AT interaction (weighted by &lambda;<sub>i</sub>&lambda;<sub>j</sub>).

| | |
|---|---|
| **Type** | `bool` |
| **Default** | `false` |

### `cross_interactions.angles`

List of cross-CG angle definitions.

**Type**: `list[CrossAngle]`

#### `cross_interactions.angles[].params`

Angle parameters as a space-separated string: `"style_id theta0 K"` where
`theta0` is the equilibrium angle in degrees and `K` is the force constant.

| | |
|---|---|
| **Type** | `string` |
| **Required** | yes |
| **Example** | `"1 111.0 530.0"` |

#### `cross_interactions.angles[].triples`

List of atom triples defining which angles to create. Each triple is a
three-element list of `MOLECULE:ATOM` identifiers.

| | |
|---|---|
| **Type** | `list[list[string]]` |
| **Required** | yes |

#### `cross_interactions.angles[].table`

Path to a tabulated angle potential file (if using tabulated angles).

| | |
|---|---|
| **Type** | `string` or `null` |
| **Default** | `null` |

#### `cross_interactions.angles[].cg_bonded`

Whether this angle uses CG weighting (1 - &lambda;&lambda;) or AT weighting
(&lambda;&lambda;).

| | |
|---|---|
| **Type** | `bool` |
| **Default** | `false` |

### `cross_interactions.dihedrals`

List of cross-CG dihedral definitions. Same structure as angles but with
`quadruples` instead of `triples`.

**Type**: `list[CrossDihedral]`

| Field | Type | Description |
|-------|------|-------------|
| `params` | `string` | Dihedral parameters |
| `quadruples` | `list[list[string]]` | Atom quadruples |
| `table` | `string` or `null` | Tabulated potential file |
| `cg_bonded` | `bool` | CG weighting flag (default: `false`) |

### Example

```yaml
cross_interactions:
  bonds:
    # CG tabulated cross bonds
    - params: "8 1 1.0"
      table: table_b1.xvg
      pairs:
        - [DOD:A1, DOD:B1]
        - [DOD:B1, DOD:B2]
      cg_bonded: true
    # AT harmonic cross bonds
    - params: "1 0.153 334720.0"
      pairs:
        - [DOD:C2, DOD:C3]
        - [DOD:C4, DOD:C5]

  angles:
    - params: "1 111.0 530.0"
      triples:
        - [DOD:C1, DOD:C2, DOD:C3]
        - [DOD:C2, DOD:C3, DOD:C4]
```

---

## `simulation`

Backmapping simulation parameters. All length and energy values are in
GROMACS units (nm, kJ/mol) in the settings file; `backmap-prep` converts
them to LAMMPS units (angstrom, kcal/mol for `real` units) in the generated
input files.

### Lambda Control

#### `simulation.alpha`

Lambda ramp rate. Lambda increases by this amount at each timestep.

| | |
|---|---|
| **Type** | `float` |
| **Default** | `0.001` |
| **Constraint** | must be positive |

#### `simulation.initial_resolution`

Initial lambda value for all atoms.

| | |
|---|---|
| **Type** | `float` |
| **Default** | `0.0` |

#### `simulation.nonuniform_lambda`

Whether to use staggered (non-uniform) initial lambda values. When `true`,
each atom gets a random offset applied to its starting lambda.

| | |
|---|---|
| **Type** | `bool` |
| **Default** | `false` |

### Time Integration

#### `simulation.timestep`

Timestep for the **optional post-backmap production** segment (Phase 3), when
`production_steps` > 0 (ps). It is **not** used for the λ-frozen equilibration
segment (Phase 1) or the λ ramp (Phase 2); those use `timestep_backmapping`.

| | |
|---|---|
| **Type** | `float` |
| **Default** | `0.001` |

#### `simulation.timestep_backmapping`

Timestep for the **first dynamics segment** in the generated input (ps): the
λ-frozen equilibration when `equilibration_steps` > 0, and the λ ramp. Set
smaller than `timestep` if the hybrid is stiff or dense. The writer also emits
`neigh_modify delay 0 every 1 check yes` after the AT/CG groups.

| | |
|---|---|
| **Type** | `float` |
| **Default** | `0.001` |

### Run Length

#### `simulation.equilibration_steps`

Optional in-hybrid relaxation with λ frozen (`fix_modify bm active no`): hybrid
CG–AT coupling (COM tracking, force distribution) continues; only the λ ramp
is paused. CG sites are not integrated directly — apply thermostats to AT
groups as usual. Runs **before** the λ ramp, using `timestep_backmapping`.
Use `0` when the CG melt was equilibrated before building the hybrid and the
hybrid starts stable; use a few thousand steps for dense melts if you see
missing-bond errors early in the ramp.

| | |
|---|---|
| **Type** | `int` |
| **Default** | `0` |

#### `simulation.production_steps`

Optional extra MD steps after λ reaches 1 **in the same generated input file**.
Use `0` (default) for **backmapping only**: the writer emits `write_data …_hybrid.data`
after the ramp; run pure-atomistic production (e.g. RDF) in a **separate** LAMMPS
input on the extracted AT system. Set to a positive value only if you want a
post-backmap segment in the same script (Phase 3 when `equilibration_steps > 0`).

| | |
|---|---|
| **Type** | `int` |
| **Default** | `0` |

### Thermostat

#### `simulation.temperature`

Target temperature (K).

| | |
|---|---|
| **Type** | `float` |
| **Default** | `300.0` |
| **Constraint** | must be positive |

#### `simulation.thermostat`

Thermostat type.

| | |
|---|---|
| **Type** | `"langevin"` or `"velocity_rescaling"` |
| **Default** | `"langevin"` |

#### `simulation.thermostat_gamma`

Thermostat damping parameter (1/ps for Langevin).

| | |
|---|---|
| **Type** | `float` |
| **Default** | `0.5` |

#### `simulation.thermostat_target`

Which atoms the thermostat is applied to.

| | |
|---|---|
| **Type** | `"atomistic"`, `"all"`, or `"cg_only"` |
| **Default** | `"atomistic"` |

### Cutoffs

#### `simulation.lj_cutoff`

Lennard-Jones (AT non-bonded) interaction cutoff (nm).

| | |
|---|---|
| **Type** | `float` |
| **Default** | `1.2` |

#### `simulation.cg_cutoff`

CG non-bonded interaction cutoff (nm).

| | |
|---|---|
| **Type** | `float` |
| **Default** | `1.4` |

#### `simulation.coulomb_cutoff`

Coulomb interaction cutoff (nm).

| | |
|---|---|
| **Type** | `float` |
| **Default** | `0.9` |

### Neighbor List and Exclusions

#### `simulation.table_groups`

Pairs of CG bead types that interact via tabulated CG potentials. Used to
set up the `pair_coeff` entries. Table files are looked up as
`table_<a>_<b>.table` first, falling back to `table_<a>_<b>.xvg` (converted).
When `cg_system.format: lammps`, use the numeric type-ID strings (e.g.
`"1"`, `"2"`) here, matching `molecules[].beads[].type`.

| | |
|---|---|
| **Type** | `list[list[string]]` |
| **Default** | `[]` |
| **Example** | `[[A, B]]` or `[A, B]` (flat list treated as single group) |

#### `simulation.exclusion_nrexcl`

Number of bonds to exclude from non-bonded interactions (LAMMPS
`special_bonds` exclusion depth).

| | |
|---|---|
| **Type** | `int` |
| **Default** | `3` |

### Output Frequency

#### `simulation.energy_interval`

How often to output thermodynamic data (steps).

| | |
|---|---|
| **Type** | `int` |
| **Default** | `1000` |

#### `simulation.trajectory_interval`

How often to dump trajectory frames (steps).

| | |
|---|---|
| **Type** | `int` |
| **Default** | `1000` |

### Random Number Generator

#### `simulation.rng_seed`

Seed for the random number generator (used by thermostat and non-uniform
lambda initialization). Set to `-1` for automatic seeding.

| | |
|---|---|
| **Type** | `int` |
| **Default** | `-1` |

---

## `output`

Output file configuration.

### `output.prefix`

Prefix for generated files. Produces `<prefix>.data` and `in.<prefix>`.

| | |
|---|---|
| **Type** | `string` |
| **Default** | `"system"` |

### `output.format`

Output format. Currently only LAMMPS is supported.

| | |
|---|---|
| **Type** | `"lammps"` |
| **Default** | `"lammps"` |

### `output.units`

LAMMPS unit system for generated files.

| | |
|---|---|
| **Type** | `"real"` |
| **Default** | `"real"` |

### Example

```yaml
output:
  prefix: dodecane
  format: lammps
  units: real
```

---

## Complete Example

```yaml
molecules:
  - name: DOD
    ident: DOD
    source:
      coordinates: dodecane_single.gro
      topology: topol_aa.top
    beads:
      - name: A1
        type: A
        atoms: [1:DOD:C1, 1:DOD:C2]
      - name: B1
        type: B
        atoms: [1:DOD:C3, 1:DOD:C4]
      - name: B2
        type: B
        atoms: [1:DOD:C5, 1:DOD:C6]

cg_system:
  coordinates: cg_conf.gro
  topology: topol_cg.top

cross_interactions:
  bonds:
    - params: "8 1 1.0"
      table: table_b1.xvg
      pairs:
        - [DOD:A1, DOD:B1]
        - [DOD:B1, DOD:B2]
      cg_bonded: true
    - params: "1 0.153 334720.0"
      pairs:
        - [DOD:C2, DOD:C3]
        - [DOD:C4, DOD:C5]

  angles:
    - params: "1 111.0 530.0"
      triples:
        - [DOD:C1, DOD:C2, DOD:C3]

simulation:
  alpha: 0.0001
  timestep: 0.001
  equilibration_steps: 0
  production_steps: 0
  temperature: 298.0
  thermostat: langevin
  thermostat_gamma: 30.0
  lj_cutoff: 1.4
  cg_cutoff: 1.4
  rng_seed: 12345

output:
  prefix: dodecane
  format: lammps
  units: real
```
