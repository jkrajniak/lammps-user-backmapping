## MODIFIED Requirements

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
