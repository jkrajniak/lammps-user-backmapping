## ADDED Requirements

### Requirement: Adaptive-resolution lambda ramp

The fix SHALL maintain a single global resolution parameter, `lambda_global`, incremented by a configurable rate each timestep: `lambda_global(t+dt) = min(1.0, lambda_global(t) + rate)`. The rate is specified as `alpha` in the fix command. `lambda_global` is the sole authoritative value every `backmap/*` interaction style weights by — there is no independent per-molecule or per-atom lambda state.

The fix SHALL support activating and deactivating the lambda ramp at runtime via `fix_modify`. When deactivated, `lambda_global` SHALL remain frozen at its current value.

Reference: Krajniak et al., "Generic Adaptive Resolution Method for Reverse Mapping of Polymers from Coarse-Grained to Atomistic Descriptions", JCTC 2016, DOI: 10.1021/acs.jctc.6b00595

#### Scenario: Uniform lambda ramp from 0 to 1
- **WHEN** the fix is configured with `alpha 0.001` and `lambda0 0.0`
- **THEN** after 1000 timesteps, `lambda_global` SHALL equal 1.0

#### Scenario: Lambda is clamped at 1.0
- **WHEN** `lambda_global` reaches 1.0
- **THEN** it SHALL remain at 1.0 on subsequent timesteps and the fix SHALL not modify it further

#### Scenario: Lambda ramp deactivated
- **WHEN** the user issues `fix_modify backmap active no` at `lambda_global` = 0.5
- **THEN** `lambda_global` SHALL remain at 0.5 until reactivated
- **AND** CG COM tracking and CG→AT force distribution SHALL continue on every timestep (only the increment in `end_of_step()` is frozen)

### Requirement: CG-AT molecule mapping

The fix SHALL build an internal mapping between CG beads and their AT atoms using LAMMPS molecule IDs (`atom->molecule`) and atom tag ordering. Within each molecule, CG atoms are identified by matching any of the user-specified CG types (one or more types via the `cg_type` parameter). All other atoms in that molecule are treated as AT atoms.

The mapping SHALL be per-CG-bead, not per-molecule. Each CG bead maps to a distinct subset of AT atoms. The mapping is inferred from atom tag ordering within each molecule:

1. Collect all CG atoms in the molecule, sorted by global tag.
2. Collect all AT atoms in the molecule, sorted by global tag.
3. Compute `atoms_per_bead = num_AT / num_CG`.
4. Assign AT atoms `[i * atoms_per_bead, (i+1) * atoms_per_bead)` to CG bead `i`.

The fix SHALL validate that `num_AT % num_CG == 0` for each molecule and abort with an error if not.

The mapping SHALL be rebuilt whenever atoms migrate between processors (after domain decomposition).

#### Scenario: Multi-bead molecule mapping (dodecane)
- **WHEN** a molecule has 6 CG beads (types 1,2,2,2,2,1) and 12 AT atoms (types 3,4,...), all in the same molecule ID
- **THEN** the fix SHALL create 6 bead mappings, each with 2 AT atoms, assigned by tag order

#### Scenario: Single-bead molecule mapping (backward compatible)
- **WHEN** a molecule has exactly 1 CG atom and N AT atoms
- **THEN** the fix SHALL create 1 bead mapping with all N AT atoms, identical to the previous behavior

#### Scenario: Missing CG atom in molecule
- **WHEN** a molecule does not contain any atom matching a CG type
- **THEN** the fix SHALL skip that molecule (it is not part of the backmapping system)

#### Scenario: Uneven AT/CG ratio
- **WHEN** a molecule has 6 CG beads and 13 AT atoms (not evenly divisible)
- **THEN** the fix SHALL abort with: "fix backmap: molecule M has N AT atoms and K CG beads, AT count must be divisible by CG count"

#### Scenario: CG mass consistency per bead
- **WHEN** the fix initializes
- **THEN** it SHALL verify that each CG bead's type mass equals the sum of its mapped AT atoms' type masses (within tolerance 1e-4) and issue a warning if not

### Requirement: CG position tracks COM of AT atoms

After each timestep, the fix SHALL update each CG bead's position to the center of mass of its own mapped AT atoms (not all AT atoms in the molecule). The COM calculation SHALL use the minimum image convention to handle atoms near periodic boundaries.

Algorithm (matching ESPResSo++ `VelocityVerletHybrid::updateVS()`):
```
for each bead mapping:
    com_delta = (0, 0, 0)
    m_total = 0
    for each AT atom i mapped to this CG bead:
        dr = minimum_image(r_AT_i - r_CG)
        com_delta += m_i * dr
        m_total += m_i
    com_delta /= m_total
    r_CG += com_delta
```

#### Scenario: COM update per bead
- **WHEN** a molecule has 6 CG beads, each with 2 AT atoms
- **THEN** each CG bead's position SHALL be the mass-weighted COM of only its own 2 AT atoms

#### Scenario: Molecule crossing periodic boundary
- **WHEN** AT atoms of a single CG bead span a periodic boundary
- **THEN** the COM calculation SHALL use minimum image vectors and produce the correct unwrapped COM position

### Requirement: CG force distribution to AT atoms

In `post_force()`, the fix SHALL distribute each CG bead's forces to only its own mapped AT atoms, proportional to their mass ratio. This replicates ESPResSo++ `VelocityVerletHybrid::distributeVSforces()`. The denominator SHALL be the per-bead sum of mapped AT type masses (`at_mass_sum`), not the CG bead's LAMMPS type mass, so that the distributed forces sum exactly to the CG force even when tabulated CG mass differs slightly from the atomistic fragment.

The CG forces arriving at `post_force()` are already lambda-weighted by the interaction styles (`pair_style backmap`, `bond_style backmap/*`, etc.). The fix distributes these pre-weighted forces without additional scaling:

```
for each bead mapping:
    for each AT atom i mapped to this CG bead:
        f[AT_i] += (m_i / at_mass_sum_bead) * f[CG_bead]
    f[CG_bead] = (0, 0, 0)
```

`at_mass_sum_bead` is computed at bead-map build time and communicated to ghost CG replicas via forward comm (same payload slot previously used for CG type mass).

After distribution, the CG bead's forces SHALL be zeroed (CG atoms are not integrated).

The fix does NOT perform any lambda scaling of forces — that responsibility belongs entirely to the interaction styles (`pair_style backmap`, `bond_style backmap/*`, `angle_style backmap/*`, `dihedral_style backmap/*`).

#### Scenario: Force distribution per bead (dodecane)
- **WHEN** CG bead A (mass 29.062) has forces (Fx, Fy, Fz), with AT atoms CH3 (mass 15.035) and CH2 (mass 14.027)
- **THEN** CH3 SHALL receive (15.035/29.062) of the force, CH2 SHALL receive (14.027/29.062), and bead B's forces are NOT distributed to bead A's AT atoms

#### Scenario: Unequal mass AT atoms (water, single-bead)
- **WHEN** a water molecule has OW (mass 16.0), HW1 (mass 1.0), HW2 (mass 1.0) and WCG (mass 18.0)
- **THEN** OW SHALL receive 16/18 of the CG force, each H SHALL receive 1/18

#### Scenario: CG forces zeroed after distribution
- **WHEN** CG forces have been distributed to AT atoms
- **THEN** all CG bead force components SHALL be zero

### Requirement: CG atoms excluded from integration

CG atoms SHALL NOT be integrated by standard LAMMPS integrators. The fix SHALL zero CG atom velocities and forces in `initial_integrate()` as a safety mechanism. The fix SHALL identify CG atoms by matching any of the configured CG types. The user MUST apply integration fixes (`fix nve`, `fix langevin`, etc.) only to AT atom groups.

#### Scenario: CG atoms remain stationary between COM updates
- **WHEN** `fix nve` is applied only to the AT group
- **THEN** CG atoms SHALL only move when the fix updates their position to COM in `end_of_step()`

#### Scenario: Safety zeroing with multiple CG types
- **WHEN** CG types are 1 and 2, and both type-1 and type-2 atoms exist
- **THEN** the fix SHALL zero velocities and forces for ALL atoms matching any CG type

#### Scenario: Safety zeroing if CG accidentally integrated
- **WHEN** a user accidentally includes CG atoms in an integration group
- **THEN** the fix SHALL zero CG velocities in `initial_integrate()` to prevent spurious CG motion

### Requirement: Lambda accessible for output

`lambda_global` SHALL be accessible to interaction styles via `fix->extract("lambda_global", dim)` (a pointer to the scalar, `dim=0`). It is additionally mirrored into a per-atom vector purely for output — every atom's mirrored value always equals `lambda_global` — accessible via `f_ID` in `dump custom` or thermo output. No interaction style reads the per-atom mirror; only `lambda_global` (plus `atom2cg` bead membership) determines any interaction weight.

The fix SHALL expose `lambda_global` as a global scalar via `compute_scalar()`, so it can be printed in thermo output as `f_bm` (for fix ID `bm`) with an optional `thermo_modify colname f_bm lambda` label.

`lambda_global` is NOT restart-persisted — the fix carries no state across `write_restart`/`read_restart`. A continuation run must reissue `fix backmap` with an explicit `lambda0` matching the desired resumption point.

#### Scenario: Thermo lambda output
- **WHEN** the input contains `thermo_style custom ... f_bm` and `fix bm` is active
- **THEN** each thermo line SHALL report the current value of `lambda_global`

#### Scenario: Dump lambda values
- **WHEN** the user specifies `dump custom ... f_bm` (or equivalent accessor)
- **THEN** the output SHALL contain `lambda_global`'s current value, repeated for every atom

### Requirement: Fix command syntax

The fix SHALL be invoked with the following syntax:
```
fix ID group-ID backmap cg_type T1 [T2 ...] alpha A lambda0 L0 [apb T1:N1 T2:N2 ...]
```

Where:
- `T1 [T2 ...]` = one or more atom types identifying CG particles (integers, read until next recognized keyword)
- `A` = lambda increment per timestep (float)
- `L0` = initial lambda value (float, default 0.0)

The `cg_type` keyword SHALL accept multiple type IDs. Parsing SHALL read integer values until a recognized keyword (`alpha`, `lambda0`, `apb`) or a non-integer token is encountered.

#### Scenario: Single CG type (backward compatible)
- **WHEN** the user specifies `fix bm all backmap cg_type 3 alpha 0.001`
- **THEN** the fix SHALL initialize with one CG type (3), lambda0=0.0

#### Scenario: Multiple CG types
- **WHEN** the user specifies `fix bm all backmap cg_type 1 2 alpha 0.0001 lambda0 0.0`
- **THEN** the fix SHALL initialize with CG types {1, 2}

#### Scenario: Invalid CG type
- **WHEN** the user specifies `cg_type 99` and atom->ntypes is 4
- **THEN** the fix SHALL abort with: "fix backmap cg_type 99 out of range [1,4]"

#### Scenario: No CG types provided
- **WHEN** the user specifies `fix bm all backmap cg_type alpha 0.001` (alpha parsed as a type)
- **THEN** the fix SHALL abort because "alpha" is not a valid integer type
