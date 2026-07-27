## MODIFIED Requirements

### Requirement: Lambda accessible for output

The per-atom λ values SHALL be accessible via `fix->extract()` for use with `dump custom` or other LAMMPS analysis tools. The fix SHALL support writing λ to restart files so that simulations can be continued.

The fix SHALL expose the group-averaged λ as a global scalar via `compute_scalar()`, so it can be printed in thermo output as `f_bm` (for fix ID `bm`) with an optional `thermo_modify colname f_bm lambda` label.

The fix SHALL additionally expose, via `fix->extract()`:
- `"atom2cg"`: a pointer to the per-atom CG-bead-membership map (`int*`,
  `dim=1`), giving the bead index that each local/ghost atom belongs to.
  This is the same mapping the fix already builds internally for COM
  tracking and force distribution.
- `"lambda_global"`: a pointer to a single authoritative global λ scalar
  (`double*`, `dim=0`), incremented by the same ramp rule as the per-atom
  `lambda[]` array (`lambda_global = min(1.0, lambda_global + alpha)` per
  timestep while the ramp is active), independent of any per-atom
  non-uniform variation.

Interaction styles (`pair_style backmap`, `bond/angle/dihedral_style
backmap/*`, `fix backmap/pairs`) SHALL read `"atom2cg"` and
`"lambda_global"` for lambda-weighting purposes instead of the per-atom
`"lambda"` array. The per-atom `"lambda"` array remains available and
continues to back thermo (`f_bm`), dump output, and restart.

#### Scenario: Thermo average lambda
- **WHEN** the input contains `thermo_style custom ... f_bm` and `fix bm` is active
- **THEN** each thermo line SHALL report the arithmetic mean of λ over atoms in the fix group

#### Scenario: Dump lambda values
- **WHEN** the user specifies `dump custom ... f_backmap[1]` (or equivalent accessor)
- **THEN** the output SHALL contain the current λ value for each atom

#### Scenario: Restart with lambda state
- **WHEN** a simulation is restarted from a restart file
- **THEN** λ values SHALL be restored to their values at the time the restart file was written

#### Scenario: atom2cg accessible via extract
- **WHEN** an interaction style calls `fix->extract("atom2cg", dim)`
- **THEN** it SHALL receive a pointer to the per-atom bead-index array and `dim` SHALL be set to 1
- **AND** the array SHALL be valid for both local and ghost atom indices (forward-communicated the same way per-atom `lambda[]` is)

#### Scenario: lambda_global accessible via extract
- **WHEN** an interaction style calls `fix->extract("lambda_global", dim)`
- **THEN** it SHALL receive a pointer to a single `double` holding the current global λ and `dim` SHALL be set to 0

#### Scenario: lambda_global tracks the same ramp as per-atom lambda
- **WHEN** the fix is configured with `alpha 0.001` and `lambda0 0.0`, with the ramp active
- **THEN** after N timesteps, `lambda_global` SHALL equal `min(1.0, N * 0.001)`, matching the per-atom `lambda[i]` value for a uniformly-ramped (non-`nonuniform`) system

#### Scenario: lambda_global unaffected by nonuniform per-atom staggering
- **WHEN** the fix is configured with `nonuniform yes`, causing per-atom `lambda[i]` values to be staggered
- **THEN** `lambda_global` SHALL still follow the single uniform ramp rule (not an average or other aggregate of the staggered per-atom values)
