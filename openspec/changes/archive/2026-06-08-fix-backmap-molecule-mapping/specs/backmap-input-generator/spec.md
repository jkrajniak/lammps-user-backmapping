## MODIFIED Requirements

### Requirement: LAMMPS input script generation

The generator SHALL produce a LAMMPS input script (`.in`) configured for backmapping. The `fix backmap` command SHALL list ALL CG atom type IDs after the `cg_type` keyword:

```
fix bm all backmap cg_type T1 T2 ... alpha A lambda0 L0 nonuniform yes/no
```

Where `T1 T2 ...` are all atom type IDs marked as CG in the system's atom type list.

The script SHALL include:
- `units real` and `atom_style full`
- `read_data` for the generated data file
- `pair_style backmap` with AT and CG sub-styles
- `pair_coeff` for all type pairs
- Bond/angle/dihedral styles with correct routing
- Group definitions for AT and CG atoms
- `fix backmap` with ALL CG type IDs, not just one
- `fix nve` and `fix langevin` applied to AT group only
- `special_bonds` from `simulation.exclusion_nrexcl`
- Three-phase run sequence

#### Scenario: Single CG type system (water)
- **WHEN** the system has one CG type (WCG, type 1)
- **THEN** the fix command SHALL be `fix bm all backmap cg_type 1 alpha ...`

#### Scenario: Multiple CG type system (dodecane)
- **WHEN** the system has CG types A (type 1) and B (type 2)
- **THEN** the fix command SHALL be `fix bm all backmap cg_type 1 2 alpha ...`

#### Scenario: CG types listed in ascending order
- **WHEN** CG types are 3, 1, 5
- **THEN** the fix command SHALL list them as `cg_type 1 3 5` (sorted ascending)

### Requirement: LAMMPS data file generation

The generator SHALL produce a LAMMPS data file (`atom_style full`) where atom coordinates are wrapped inside the simulation box `[0, L)` for each dimension. Atoms with coordinates outside the box SHALL be wrapped using periodic boundary conditions before writing.

Particle ordering within each molecule SHALL be: all CG atoms first (in bead order), then all AT atoms (in bead order, with atoms within each bead in their source file order).

#### Scenario: Atoms outside box are wrapped
- **WHEN** an AT atom has z-coordinate -1.57 and box is [0, 57.54]
- **THEN** the data file SHALL contain z = 55.97 (wrapped coordinate)

#### Scenario: Atoms inside box are unchanged
- **WHEN** an AT atom has z-coordinate 30.0 and box is [0, 57.54]
- **THEN** the data file SHALL contain z = 30.0

#### Scenario: Atom ordering convention
- **WHEN** a dodecane molecule has 6 CG beads and 12 AT atoms
- **THEN** the data file SHALL list atoms 1-6 as CG (in bead order) and atoms 7-18 as AT (grouped by bead, in bead order)
