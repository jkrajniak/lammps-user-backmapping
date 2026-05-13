## MODIFIED Requirements

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

Algorithm:
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

In `post_force()`, the fix SHALL distribute each CG bead's forces to only its own mapped AT atoms, proportional to their mass ratio. The mass ratio uses the CG bead's type mass as the denominator:

```
for each bead mapping:
    for each AT atom i mapped to this CG bead:
        f[AT_i] += (m_i / m_CG_bead) * f[CG_bead]
    f[CG_bead] = (0, 0, 0)
```

After distribution, the CG bead's forces SHALL be zeroed.

#### Scenario: Force distribution per bead (dodecane)
- **WHEN** CG bead A (mass 29.062) has forces (Fx, Fy, Fz), with AT atoms CH3 (mass 15.035) and CH2 (mass 14.027)
- **THEN** CH3 SHALL receive (15.035/29.062) of the force, CH2 SHALL receive (14.027/29.062), and bead B's forces are NOT distributed to bead A's AT atoms

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

### Requirement: Fix command syntax

The fix SHALL be invoked with the following syntax:
```
fix ID group-ID backmap cg_type T1 [T2 ...] alpha A lambda0 L0 [nonuniform yes/no] [phase P]
```

Where:
- `T1 [T2 ...]` = one or more atom types identifying CG particles (integers, read until next recognized keyword)
- `A` = lambda increment per timestep (float)
- `L0` = initial lambda value (float, default 0.0)
- `nonuniform` = staggered initial lambda (optional, default no)
- `P` = initial phase (optional, integer 1 or 2, default 2)

The `cg_type` keyword SHALL accept multiple type IDs. Parsing SHALL read integer values until a recognized keyword (`alpha`, `lambda0`, `nonuniform`, `phase`) or a non-integer token is encountered.

#### Scenario: Single CG type (backward compatible)
- **WHEN** the user specifies `fix bm all backmap cg_type 3 alpha 0.001`
- **THEN** the fix SHALL initialize with one CG type (3), lambda0=0.0, nonuniform=no, phase=2

#### Scenario: Multiple CG types
- **WHEN** the user specifies `fix bm all backmap cg_type 1 2 alpha 0.0001 lambda0 0.0`
- **THEN** the fix SHALL initialize with CG types {1, 2}

#### Scenario: Invalid CG type
- **WHEN** the user specifies `cg_type 99` and atom->ntypes is 4
- **THEN** the fix SHALL abort with: "fix backmap cg_type 99 out of range [1,4]"

#### Scenario: No CG types provided
- **WHEN** the user specifies `fix bm all backmap cg_type alpha 0.001` (alpha parsed as a type)
- **THEN** the fix SHALL abort because "alpha" is not a valid integer type
