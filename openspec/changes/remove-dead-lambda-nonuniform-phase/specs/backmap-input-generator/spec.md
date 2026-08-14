## MODIFIED Requirements

### Requirement: Simulation parameters section

The `simulation` section SHALL configure backmapping parameters:

```yaml
simulation:
  # Dynamic resolution
  alpha: 0.001                       # lambda increment per step
  initial_resolution: 0.0            # initial lambda value

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

`two_phase: true` is currently rejected at validation time (deferred feature, not yet implemented) rather than producing Phase 1/Phase 2 run blocks — see the "Two-phase backmapping is rejected" scenario below. `alpha2`, `cap_force`, and `cap_force_ramp` are unreachable while `two_phase` is rejected.

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

#### Scenario: Two-phase parameters (unreachable while two_phase is rejected)
- **WHEN** `two_phase: true` and `alpha2: 0.0002` are set
- **THEN** validation SHALL abort on `two_phase` before `alpha2` is ever consulted (see "Two-phase backmapping is rejected" above)

#### Scenario: Force capping (unreachable while two_phase is rejected)
- **WHEN** `cap_force: 1000.0` is set (kJ/(mol·nm)) alongside `two_phase: true`
- **THEN** validation SHALL abort on `two_phase` before `cap_force` is ever consulted


### Requirement: LAMMPS input script generation

The generator SHALL produce a LAMMPS input script (`.in`) configured for backmapping. The `fix backmap` command SHALL list ALL CG atom type IDs after the `cg_type` keyword:

```
fix bm all backmap cg_type T1 T2 ... alpha A lambda0 L0
```

Where `T1 T2 ...` are all atom type IDs marked as CG in the system's atom type list, sorted ascending.

The script SHALL include:
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
- `fix backmap` with ALL CG type IDs (not just one) plus parameters from `simulation` section
- `fix nve` and `fix langevin` applied to AT group only (or as configured by `thermostat_target`)
- Two-phase backmapping (`simulation.two_phase`) remains deferred (see "Deferred features" below) and is rejected at validation time; `fix backmap`'s command syntax has no `phase` keyword to target for it
- `special_bonds` from `simulation.exclusion_nrexcl`
- Thermo output every `simulation.energy_interval` steps
- Dump configuration every `simulation.trajectory_interval` steps
- Three-phase run sequence: CG equilibration -> backmapping -> AT production

#### Scenario: Single CG type system (water)
- **WHEN** the system has one CG type (WCG, type 1)
- **THEN** the fix command SHALL be `fix bm all backmap cg_type 1 alpha ...`

#### Scenario: Multiple CG type system (dodecane)
- **WHEN** the system has CG types A (type 1) and B (type 2)
- **THEN** the fix command SHALL be `fix bm all backmap cg_type 1 2 alpha ...`

#### Scenario: CG types listed in ascending order
- **WHEN** CG types are 3, 1, 5
- **THEN** the fix command SHALL list them as `cg_type 1 3 5` (sorted ascending)

#### Scenario: Water2 system (no cross bonds)
- **WHEN** the settings describe a single-bead water system with no cross interactions
- **THEN** the `.in` file SHALL use static `bond_style harmonic` and `angle_style harmonic`, `pair_style backmap` for non-bonded, and `fix backmap` for lambda ramp

#### Scenario: Polymer system (with cross bonds)
- **WHEN** the settings describe a polymer with cross bonds, cross angles, and cross dihedrals
- **THEN** the `.in` file SHALL use `bond_style hybrid harmonic backmap/harmonic backmap/table`, routing each bond type to the correct sub-style

#### Scenario: Polymer system with cross dihedrals
- **WHEN** the settings describe a polymer with cross bonds, cross angles, and cross dihedrals (both RB and tabulated)
- **THEN** the `.in` file SHALL contain `dihedral_style hybrid backmap/ryckaert backmap/table` and correct `dihedral_coeff` lines for each type

#### Scenario: Two-phase backmapping is rejected (deferred feature)
- **WHEN** `simulation.two_phase: true` is set
- **THEN** the generator SHALL abort at validation time with "Feature 'two_phase' backmapping is not yet implemented (planned for Phase 2)" — `fix backmap`'s command syntax has no `phase` keyword, so this deferred feature has no current target mechanism

#### Scenario: System with dihedrals but no tabulated dihedrals
- **WHEN** the settings define only RB cross dihedrals (no tabulated)
- **THEN** the `.in` file SHALL use `dihedral_style hybrid backmap/ryckaert` (without `backmap/table`)

#### Scenario: Simulation parameters applied
- **WHEN** settings specify `timestep: 0.001` (ps), `temperature: 423.0`, `alpha: 0.0005`
- **THEN** the `.in` file SHALL contain `timestep 1.0` (fs), `fix langevin ... 423.0 423.0 ...`, `fix backmap ... alpha 0.0005`
