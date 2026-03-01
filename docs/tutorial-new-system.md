# Tutorial: Setting Up a New System

This tutorial walks you through setting up backmapping for a new molecular
system from scratch. We use dodecane (C12H26) as a concrete example, but the
steps apply to any molecule.

## Prerequisites

Before starting, make sure you have:

- The LAMMPS Backmapping Package [installed](getting-started.md)
- The `backmap-prep` CLI available (`uv run backmap-prep --help`)
- A coarse-grained simulation of your system (coordinates + topology)
- An atomistic reference for your molecule (coordinates + topology)

Both coordinate and topology files must be in GROMACS format (`.gro` and
`.top` files).

## Step 1: Prepare Your Input Files

You need four input files:

| File | Description |
|------|-------------|
| AT coordinate file | Single-molecule atomistic structure (`.gro`) |
| AT topology file | Atomistic force field and topology (`.top`) |
| CG coordinate file | Full CG system coordinates (`.gro`) |
| CG topology file | CG system topology (`.top`) |

### Atomistic Reference

The AT coordinate file should contain a **single molecule** with all atoms
present. This serves as the template for placing AT atoms inside CG beads.

```
dodecane single molecule
   12
    1DOD     C1    1   1.000   1.000   1.000
    1DOD     C2    2   1.153   1.000   1.000
    1DOD     C3    3   1.306   1.000   1.000
    ...
    1DOD    C12   12   2.683   1.000   1.000
   5.0  5.0  5.0
```

The AT topology file defines the atomistic force field (atom types, masses,
charges, bonds, angles, dihedrals):

```
[ moleculetype ]
DOD    3

[ atoms ]
     1   CH3  1  DOD  C1   1  0.000  15.035
     2   CH2  1  DOD  C2   2  0.000  14.027
     ...
```

### CG System

The CG coordinate file is a snapshot from your CG simulation. It should
contain all CG beads in the simulation box.

The CG topology file defines bead types, masses, and CG bonded interactions.

!!! tip
    Make sure the residue names and atom names in your files match what
    you reference in the settings YAML. The `RESID:RESNAME:ATOMNAME` format
    in the `atoms` field must correspond to entries in the AT coordinate
    and topology files.

## Step 2: Define the CG-to-AT Mapping

Create a `settings.yaml` file. Start with the molecule definitions:

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

Each bead entry maps a CG bead name to the AT atoms it contains. The `type`
field groups beads that share the same non-bonded CG interactions.

**Key decisions:**

- How many AT atoms per bead? This should match your CG model's mapping.
- What are the bead types? Beads with different CG non-bonded potentials need
  different types.
- What are the atom names? They must match the AT topology exactly.

## Step 3: Specify the CG System Files

```yaml
cg_system:
  coordinates: cg_conf.gro
  topology: topol_cg.top
  format: gromacs
```

## Step 4: Define Cross Interactions

Cross interactions are the bonded terms that span CG bead boundaries. There
are two kinds:

### CG Cross Bonds (fade out during backmapping)

These maintain the CG structure during the transition. They use the CG bond
potential (typically tabulated from IBI or similar methods):

```yaml
cross_interactions:
  bonds:
    - params: "8 1 1.0"
      table: table_b1.xvg
      pairs:
        - [DOD:A1, DOD:B1]
        - [DOD:B1, DOD:B2]
        - [DOD:B2, DOD:B3]
        - [DOD:B3, DOD:B4]
        - [DOD:B4, DOD:A2]
      cg_bonded: true
```

!!! note
    The `table` field points to an XVG file with the tabulated potential.
    The `params` field provides metadata for the table interpolation.
    Set `cg_bonded: true` so the bond weight follows 1 - &lambda;&lambda;.

### AT Cross Bonds (fade in during backmapping)

These are the atomistic bonds between atoms in *different* CG beads. They
restore the correct AT connectivity:

```yaml
    - params: "1 0.153 334720.0"
      pairs:
        - [DOD:C2, DOD:C3]
        - [DOD:C4, DOD:C5]
        - [DOD:C6, DOD:C7]
        - [DOD:C8, DOD:C9]
        - [DOD:C10, DOD:C11]
```

Here `params` contains the harmonic bond parameters: `style_id r0 K`.

!!! warning
    Only list bonds that cross bead boundaries. Intra-bead bonds (e.g.,
    C1-C2 within bead A1) are handled automatically from the AT topology.

### AT Cross Angles

Similarly, define angles that span bead boundaries:

```yaml
  angles:
    - params: "1 111.0 530.0"
      triples:
        - [DOD:C1, DOD:C2, DOD:C3]
        - [DOD:C2, DOD:C3, DOD:C4]
        - [DOD:C3, DOD:C4, DOD:C5]
        ...
```

**How to identify cross interactions:**

1. List all AT bonds in your molecule
2. For each bond, check if both atoms belong to the same CG bead
3. If they belong to different beads, it's a cross bond
4. Repeat for angles (check if the triplet spans bead boundaries)

## Step 5: Configure Simulation Parameters

```yaml
simulation:
  alpha: 0.0001
  initial_resolution: 0.0
  nonuniform_lambda: false

  timestep: 0.001
  timestep_backmapping: 0.001

  equilibration_steps: 10000
  production_steps: 10000

  temperature: 298.0
  thermostat: langevin
  thermostat_gamma: 30.0
  thermostat_target: atomistic

  lj_cutoff: 1.4
  cg_cutoff: 1.4
  coulomb_cutoff: 0.9

  table_groups: [A, B]

  exclusion_nrexcl: 3
  energy_interval: 1000
  trajectory_interval: 1000
  rng_seed: 12345
```

**Parameter guidelines:**

| Parameter | Guidance |
|-----------|----------|
| `alpha` | Start with 10<sup>-4</sup>; decrease if the transition is unstable |
| `equilibration_steps` | Enough to equilibrate CG structure (monitor energy) |
| `temperature` | Match your CG simulation temperature |
| `thermostat_gamma` | Higher values = stronger damping (Langevin) |
| `table_groups` | List the CG bead types that have tabulated pair potentials |

## Step 6: Set Output Options

```yaml
output:
  prefix: mysystem
  format: lammps
  units: real
```

## Step 7: Generate LAMMPS Input Files

Run `backmap-prep`:

```bash
uv run backmap-prep settings.yaml
```

This generates:

- `mysystem.data` -- LAMMPS data file containing:
    - Both CG and AT atoms
    - All bond, angle, and dihedral topologies
    - Masses and atom types
- `in.mysystem` -- LAMMPS input script with:
    - `pair_style backmap` with AT and CG sub-styles
    - `bond_style` and `angle_style` definitions
    - `fix backmap` setup
    - Three-phase simulation protocol
- Table files -- converted from GROMACS XVG to LAMMPS format

You can override the output prefix:

```bash
uv run backmap-prep settings.yaml --output-prefix myrun
```

### Inspect the Generated Files

Before running, review the generated input script to verify:

- [ ] Atom types and masses look correct
- [ ] Pair coefficients are assigned to the right type pairs
- [ ] Bond and angle coefficients match your force field
- [ ] The three simulation phases have appropriate step counts

## Step 8: Run the Simulation

```bash
lmp -in in.mysystem
```

Or in parallel:

```bash
mpirun -np 4 lmp -in in.mysystem
```

### Monitor the Run

Watch the thermodynamic output for:

- **Energy convergence** during CG equilibration (Phase 1)
- **Smooth energy transition** during backmapping (Phase 2) -- expect some
  fluctuation as AT interactions ramp up
- **Stable energy** during AT production (Phase 3)

!!! warning "Common issues"
    - **Atoms lost**: reduce `alpha` for a gentler transition
    - **High energy spikes**: check cross-interaction definitions; missing
      cross bonds cause severe distortions
    - **Temperature drift**: adjust thermostat damping (`thermostat_gamma`)

### Visualize Results

The trajectory dump includes per-atom lambda values in the `f_bm` column.
You can use this to color atoms by resolution in VMD, OVITO, or similar
visualization tools.

## Checklist

Before running a production backmapping simulation, verify:

- [ ] AT reference molecule has correct geometry and topology
- [ ] CG-to-AT mapping covers all AT atoms (no atom left out of a bead)
- [ ] All cross-boundary bonds and angles are defined
- [ ] CG tabulated potentials cover the expected distance range
- [ ] Alpha is small enough for a stable transition
- [ ] Thermostat temperature matches the CG simulation
