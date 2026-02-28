## ADDED Requirements

### Requirement: Time-dependent lambda ramp

The fix SHALL increment the per-atom resolution parameter λ by a configurable rate each timestep: `λ(t+dt) = min(1.0, λ(t) + rate)`. The rate is specified as `alpha` in the fix command. All atoms within a molecule (CG + AT) SHALL receive the same λ value.

The fix SHALL support activating and deactivating the lambda ramp at runtime via `fix_modify`. When deactivated, λ values SHALL remain frozen at their current values.

Reference: Krajniak et al., "Generic Adaptive Resolution Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic Descriptions", JCTC 2016, DOI: 10.1021/acs.jctc.6b00595

#### Scenario: Uniform lambda ramp from 0 to 1
- **WHEN** the fix is configured with `alpha 0.001` and `lambda0 0.0`
- **THEN** after 1000 timesteps, all atoms SHALL have λ = 1.0

#### Scenario: Lambda is clamped at 1.0
- **WHEN** λ reaches 1.0 for a molecule
- **THEN** λ SHALL remain at 1.0 on subsequent timesteps and the fix SHALL not modify it further

#### Scenario: Lambda ramp deactivated
- **WHEN** the user issues `fix_modify backmap active no` at λ = 0.5
- **THEN** λ SHALL remain at 0.5 for all atoms until reactivated

#### Scenario: Non-uniform initial lambda
- **WHEN** the fix is configured with `nonuniform yes`
- **THEN** each molecule SHALL receive a random initial λ in the range `[-10000*alpha, 0]`, causing molecules to reach λ=0 at staggered times. Negative λ values SHALL be treated as λ=0 for force weighting.

### Requirement: CG-AT molecule mapping

The fix SHALL build an internal mapping between CG atoms and their AT atoms using LAMMPS molecule IDs (`atom->molecule`). Within each molecule, the CG atom is identified by the user-specified `cg_type` parameter. All other atoms in that molecule are treated as AT atoms.

The mapping SHALL be rebuilt whenever atoms migrate between processors (after domain decomposition).

#### Scenario: Valid molecule mapping
- **WHEN** a LAMMPS data file contains molecules where each molecule has exactly one atom of `cg_type`
- **THEN** the fix SHALL correctly identify CG-AT groups and report the number of mapped molecules during initialization

#### Scenario: Missing CG atom in molecule
- **WHEN** a molecule does not contain any atom of `cg_type`
- **THEN** the fix SHALL issue an error and abort

#### Scenario: CG mass consistency
- **WHEN** the fix initializes
- **THEN** it SHALL verify that each CG atom's mass equals the sum of its AT atoms' masses (within tolerance 1e-4) and issue a warning if not

### Requirement: CG position tracks COM of AT atoms

After each timestep, the fix SHALL update each CG atom's position to the center of mass of its AT atoms. The COM calculation SHALL use the minimum image convention to handle molecules near periodic boundaries.

Algorithm (matching ESPResSo++ `VelocityVerletHybrid::updateVS()`):
```
for each molecule:
    com_delta = (0, 0, 0)
    for each AT atom i:
        dr = minimum_image(r_AT_i - r_CG)
        com_delta += m_i * dr
    com_delta /= m_CG
    r_CG += com_delta
```

#### Scenario: COM update after AT atoms move
- **WHEN** AT atoms move during integration
- **THEN** the CG atom position SHALL equal the mass-weighted center of mass of its AT atoms after `end_of_step()`

#### Scenario: Molecule crossing periodic boundary
- **WHEN** AT atoms of a molecule span a periodic boundary
- **THEN** the COM calculation SHALL use minimum image vectors and produce the correct unwrapped COM position

### Requirement: CG force distribution to AT atoms

In `post_force()`, the fix SHALL distribute CG forces to AT atoms proportional to their mass ratio. This replicates ESPResSo++ `VelocityVerletHybrid::distributeVSforces()`.

The CG forces arriving at `post_force()` are already lambda-weighted by the interaction styles (`pair_style backmap`, `bond_style backmap/*`, etc.). The fix distributes these pre-weighted forces without additional scaling:

```
for each molecule:
    for each AT atom i:
        f[AT_i] += (m_i / m_CG) * f[CG]
    f[CG] = (0, 0, 0)
```

After distribution, the CG atom's forces SHALL be zeroed (CG atoms are not integrated).

The fix does NOT perform any lambda scaling of forces — that responsibility belongs entirely to the interaction styles (`pair_style backmap`, `bond_style backmap/*`, `angle_style backmap/*`, `dihedral_style backmap/*`).

#### Scenario: Equal mass AT atoms
- **WHEN** a molecule has 3 AT atoms each with mass 1.0 and CG mass 3.0
- **THEN** each AT atom SHALL receive exactly 1/3 of the CG force

#### Scenario: Unequal mass AT atoms
- **WHEN** a water molecule has OW (mass 16.0), HW1 (mass 1.0), HW2 (mass 1.0) and WCG (mass 18.0)
- **THEN** OW SHALL receive 16/18 of the CG force, each H SHALL receive 1/18

#### Scenario: CG forces zeroed after distribution
- **WHEN** CG forces have been distributed to AT atoms
- **THEN** all CG atom force components SHALL be zero

### Requirement: CG atoms excluded from integration

CG atoms SHALL NOT be integrated by standard LAMMPS integrators. The fix SHALL zero CG atom velocities and forces in `initial_integrate()` as a safety mechanism. The user MUST apply integration fixes (`fix nve`, `fix langevin`, etc.) only to AT atom groups.

#### Scenario: CG atoms remain stationary between COM updates
- **WHEN** `fix nve` is applied only to the AT group
- **THEN** CG atoms SHALL only move when the fix updates their position to COM in `end_of_step()`

#### Scenario: Safety zeroing if CG accidentally integrated
- **WHEN** a user accidentally includes CG atoms in an integration group
- **THEN** the fix SHALL zero CG velocities in `initial_integrate()` to prevent spurious CG motion

### Requirement: Simulation phases

The fix SHALL support the three-phase simulation workflow: CG equilibration → backmapping → AT production.

The backmapping phase SHALL support both single-phase and two-phase modes:
- **Single phase** (default, Phase 2): All lambda-weighted interactions (pair, cross-CG bonds, cross-CG angles, cross-CG dihedrals) use their respective custom styles with AT weight `λ²` and CG weight `1−λ²`. Intra-CG bonded forces remain at full strength (static). λ ramps from 0 to 1.
- **Two-phase mode**: Controlled via `fix_modify` with a state machine:
  - **Phase 1**: λ ramps 0→1. AT cross interactions are weighted by `λ²` (ramping in). CG cross interactions remain at full strength (weight = 1.0, ignoring lambda). This establishes atomistic connectivity under full CG restraint.
  - **Phase 2**: λ resets to 0.0 and ramps 0→1 again. AT cross interactions are weighted by `λ²` (ramping in). CG cross interactions are weighted by `1−λ²` (ramping out). This is the standard AdResS interpolation.

The timestep SHALL be changeable between phases via `timestep` command.

The current phase SHALL be accessible to interaction styles via `fix->extract("phase")`, returning a pointer to an integer (1 or 2). Single-phase mode uses phase=2 (the default MVP behavior).

`fix_modify` interface for phase control:
```
fix_modify <fix-id> phase 1       # enter Phase 1
fix_modify <fix-id> phase 2       # enter Phase 2, reset λ to 0
fix_modify <fix-id> active yes    # enable lambda ramp (existing)
fix_modify <fix-id> active no     # freeze lambda (existing)
```

#### Scenario: Single-phase backmapping (default)
- **WHEN** the fix is active with no explicit phase setting
- **THEN** the fix SHALL operate in Phase 2 mode: AT non-bonded pair forces weighted by `λ²`, CG non-bonded pair forces by `1−λ²`, cross-CG AT bonded forces by `λ²`, cross-CG CG bonded forces by `1−λ²`, and intra-CG bonded forces at full strength

#### Scenario: Two-phase backmapping — Phase 1 behavior
- **WHEN** `fix_modify bm phase 1` is issued
- **THEN** the fix SHALL set phase=1, and `fix->extract("phase")` SHALL return a pointer to integer value 1
- **AND** CG cross interactions SHALL receive weight 1.0 from `backmap_lambda.h` helpers regardless of lambda

#### Scenario: Two-phase backmapping — transition from Phase 1 to Phase 2
- **WHEN** Phase 1 completes (λ reaches 1.0) and the user issues `fix_modify bm phase 2`
- **THEN** λ SHALL be reset to 0.0 for all atoms, phase SHALL be set to 2, and the lambda ramp SHALL restart
- **AND** CG cross interactions SHALL begin using `1−λ²` weighting

#### Scenario: Phase 2 lambda reset
- **WHEN** `fix_modify bm phase 2` is issued while λ values are at various values
- **THEN** ALL per-atom λ values SHALL be reset to 0.0 regardless of their current values

#### Scenario: Phase queryable via extract
- **WHEN** an interaction style calls `fix->extract("phase", size)`
- **THEN** it SHALL receive a pointer to the current phase integer (1 or 2)
- **AND** `size` SHALL be set to 0 (scalar, not per-atom)

#### Scenario: Phase preserved in restart
- **WHEN** a simulation is restarted from a restart file that was written during Phase 1
- **THEN** the phase SHALL be restored to 1 and λ values SHALL be at their saved values

#### Scenario: Invalid phase value
- **WHEN** the user issues `fix_modify bm phase 3` or any value other than 1 or 2
- **THEN** the fix SHALL abort with: "fix backmap: phase must be 1 or 2"

### Requirement: Lambda accessible for output

The per-atom λ values SHALL be accessible via `fix->extract()` for use with `dump custom` or other LAMMPS analysis tools. The fix SHALL support writing λ to restart files so that simulations can be continued.

#### Scenario: Dump lambda values
- **WHEN** the user specifies `dump custom ... f_backmap[1]` (or equivalent accessor)
- **THEN** the output SHALL contain the current λ value for each atom

#### Scenario: Restart with lambda state
- **WHEN** a simulation is restarted from a restart file
- **THEN** λ values SHALL be restored to their values at the time the restart file was written

### Requirement: Fix command syntax

The fix SHALL be invoked with the following syntax:
```
fix ID group-ID backmap cg_type T alpha A lambda0 L0 [nonuniform yes/no] [phase P]
```

Where:
- `T` = atom type identifying CG particles (integer)
- `A` = lambda increment per timestep (float)
- `L0` = initial lambda value (float, default 0.0)
- `nonuniform` = staggered initial lambda (optional, default no)
- `P` = initial phase (optional, integer 1 or 2, default 2)

#### Scenario: Minimal invocation
- **WHEN** the user specifies `fix bm all backmap cg_type 3 alpha 0.001`
- **THEN** the fix SHALL initialize with lambda0=0.0, nonuniform=no, phase=2

#### Scenario: Full invocation with two-phase
- **WHEN** the user specifies `fix bm all backmap cg_type 3 alpha 0.0005 lambda0 0.0 phase 1`
- **THEN** the fix SHALL initialize in Phase 1 with lambda0=0.0

#### Scenario: Phase specified via fix_modify after init
- **WHEN** the user starts with default phase (2) and issues `fix_modify bm phase 1`
- **THEN** the fix SHALL switch to Phase 1 and reset λ to lambda0 for all atoms

#### Scenario: Full invocation with nonuniform
- **WHEN** the user specifies `fix bm all backmap cg_type 3 alpha 0.0005 lambda0 0.0 nonuniform yes`
- **THEN** the fix SHALL initialize with staggered lambda values for each molecule
