## REMOVED Requirements

### Requirement: Simulation phases

**Reason**: Described a two-phase `fix_modify bm phase 1/2` state machine that was never implemented. `phase` was reserved in the argument parser's `KNOWN_KEYWORDS` with no handling branch (a user specifying it got a hard parse error), and `compute_weight3()`'s `phase` parameter was never passed a non-default value by any production call site. See `design.md` for the verification trail.

## MODIFIED Requirements

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
