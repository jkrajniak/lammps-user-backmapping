## Tasks

Tasks are grouped by component and ordered by dependency. Within each group, tasks build on each other sequentially. Groups can be worked in parallel where noted.

### Legend

- **[MVP]** = Phase 1, required for initial working system
- **[P2]** = Phase 2, extended bonded interactions
- **[P3]** = Phase 3, reactive networks
- **[P4]** = Phase 4, format flexibility

---

### Group A: Shared Lambda Infrastructure (C++)

These are foundational — all LAMMPS styles depend on them.

- [x] **A1 [MVP]**: Create `backmap_lambda.h` shared header
  - `find_fix_backmap()` — locate fix by scanning fix list
  - `get_lambda(int i)` — read per-atom λ from fix's array
  - `compute_weight(double lambda_i, double lambda_j, bool is_cg)` — returns `λ_i×λ_j` or `1−λ_i×λ_j`
  - `is_almost_zero(double w)` — threshold check (~1e-10)
  - Spec: `bonded-backmap-styles/spec.md` → "Shared lambda-access helper"

- [x] **A2 [MVP]**: Create LAMMPS package scaffolding (`USER-BACKMAP`)
  - `Install.sh` for the LAMMPS package system
  - `README` with package description
  - CMake integration (`CMakeLists.txt`) for modern LAMMPS builds

---

### Group B: `fix backmap` (C++)

The central fix — manages lambda, CG-AT mapping, COM tracking, force distribution.

- [x] **B1 [MVP]**: Implement `fix_backmap.cpp/.h` — constructor and initialization
  - Parse command syntax: `fix ID group-ID backmap cg_type T alpha A lambda0 L0 [nonuniform yes/no]`
  - Build CG-AT molecule mapping from `atom->molecule` and `cg_type`
  - Allocate per-atom λ array, initialize to `lambda0` (or staggered for nonuniform)
  - Validate CG mass = sum of AT masses (warning if mismatch)
  - Spec: `fix-backmap-resolution/spec.md` → "CG-AT molecule mapping", "Fix command syntax"

- [x] **B2 [MVP]**: Implement `initial_integrate()` hook
  - Zero CG atom velocities and forces (safety net)
  - Spec: `fix-backmap-resolution/spec.md` → "CG atoms excluded from integration"

- [x] **B3 [MVP]**: Implement `post_force()` hook — CG force distribution
  - For each molecule: distribute CG forces to AT atoms proportional to mass ratio
  - Zero CG forces after distribution
  - Spec: `fix-backmap-resolution/spec.md` → "CG force distribution to AT atoms"

- [x] **B4 [MVP]**: Implement `end_of_step()` hook — COM update and λ ramp
  - Update CG position to COM of AT atoms (minimum image convention)
  - Increment λ by alpha, clamp at 1.0
  - Spec: `fix-backmap-resolution/spec.md` → "CG position tracks COM of AT atoms", "Time-dependent lambda ramp"

- [x] **B5 [MVP]**: Implement λ access for other styles
  - `extract()` method returning per-atom λ array
  - `fix_modify` support for activating/deactivating ramp
  - Spec: `fix-backmap-resolution/spec.md` → "Lambda accessible for output"

- [x] **B6 [MVP]**: Implement MPI support
  - Rebuild molecule map on `post_exchange` callback (atom migration)
  - `comm_forward` / `comm_reverse` for ghost atom λ values
  - Spec: design.md → Risk "Molecule migration across MPI boundaries"

- [x] **B7 [MVP]**: Implement restart support
  - Write/read per-atom λ to/from restart files
  - Spec: `fix-backmap-resolution/spec.md` → "Restart with lambda state"

- [x] **B8 [P2]**: ~~Implement two-phase backmapping support~~ → Moved to `backmap-phase2-and-testing`

---

### Group C: `pair_style backmap` (C++)

Can be developed in parallel with Group B (needs A1 only).

- [x] **C1 [MVP]**: Implement `pair_backmap.cpp/.h` — constructor and settings
  - Parse command syntax: `pair_style backmap cutoff_at at_style at_args ... cutoff_cg cg_style cg_args ...`
  - Instantiate AT and CG sub-styles
  - Spec: `pair-backmap/spec.md` → "Pair style command syntax"

- [x] **C2 [MVP]**: Implement `pair_coeff` parsing
  - Route type pairs to `atomistic`, `cg`, or `none`
  - Auto-detect or accept fix ID for lambda access
  - Spec: `pair-backmap/spec.md` → "Pair style command syntax", "Lambda access from fix"

- [x] **C3 [MVP]**: Implement `compute()` — lambda-weighted pair forces
  - For each pair: compute sub-style force, multiply by `w_AT = λ_i×λ_j` or `w_CG = 1−λ_i×λ_j`
  - Skip computation when weight is negligible
  - Handle negative λ (clamp to 0)
  - Spec: `pair-backmap/spec.md` → "Lambda-weighted pair force computation"

- [x] **C4 [MVP]**: Implement energy and virial computation
  - Correctly weight energy/virial by lambda for thermo output
  - Spec: `pair-backmap/spec.md` → "Energy computation"

- [x] **C5 [MVP]**: Implement `single()` method
  - For `compute` command and neighbor list diagnostics
  - Spec: `pair-backmap/spec.md` → "Sub-style delegation"

---

### Group D: Bonded `backmap/*` Styles (C++)

Can be developed in parallel with Groups B and C (needs A1 only). Ordered by likely usage frequency.

- [x] **D1 [MVP]**: Implement `bond_backmap_harmonic.cpp/.h`
  - Parse: `bond_coeff N at/cg K r0`
  - Compute: `F = -w × k × (r - r0)`, `E = w × 0.5 × k × (r - r0)²`
  - Read λ from fix via `backmap_lambda.h`
  - Restart file support
  - Spec: `bonded-backmap-styles/spec.md` → "Bond style backmap/harmonic"

- [x] **D2 [MVP]**: Implement `bond_backmap_table.cpp/.h`
  - Parse: `bond_coeff M at/cg filename keyword`
  - Compute: `F = w × F_table(r)`, `E = w × E_table(r)`
  - Follow LAMMPS `bond_style table` conventions for table format
  - Spec: `bonded-backmap-styles/spec.md` → "Bond style backmap/table"

- [x] **D3 [MVP]**: Implement `angle_backmap_harmonic.cpp/.h`
  - Parse: `angle_coeff N at/cg K theta0`
  - Weight from first and last atoms of the angle triple
  - Spec: `bonded-backmap-styles/spec.md` → "Angle style backmap/harmonic"

- [x] **D4 [P2]**: ~~Implement `dihedral_backmap_ryckaert.cpp/.h`~~ → Moved to `backmap-phase2-and-testing`

- [x] **D5 [P2]**: ~~Implement `dihedral_backmap_table.cpp/.h`~~ → Moved to `backmap-phase2-and-testing`

---

### Group E: Python Input Generator (`backmap-prep`)

Can be developed in parallel with C++ groups. Depends on YAML schema design (completed in specs).

- [x] **E1 [MVP]**: Set up Python project scaffolding
  - `pyproject.toml` with dependencies: PyYAML, Pydantic v2
  - CLI entry point: `backmap-prep`
  - Spec: `backmap-input-generator/spec.md` → "Command-line interface", "Python 3 compatibility"

- [x] **E2 [MVP]**: Implement Pydantic YAML schema and validation
  - Models for all five top-level sections (molecules, cg_system, cross_interactions, simulation, output)
  - Defaults for optional fields
  - Validation errors for missing/invalid fields
  - "Not yet implemented" errors for Phase 2-4 features
  - Spec: `backmap-input-generator/spec.md` → "YAML settings file as primary input", "Feature phasing"

- [x] **E3 [MVP]**: Implement GROMACS source file parsers
  - `parse_gro()` — read atom positions and box dimensions from `.gro` files
  - `parse_top()` — read atom types, charges, masses, bonds, angles from `.top`/`.itp` files (with `#include` resolution)
  - Modular design for future format additions
  - Spec: `backmap-input-generator/spec.md` → "Source file parsing"

- [x] **E4 [MVP]**: Implement unit conversion module
  - Functions for distance, energy, force, time, spring constant conversions
  - GROMACS → LAMMPS real
  - Spec: `backmap-input-generator/spec.md` → "Unit conversion"

- [x] **E5 [MVP]**: Implement CG-AT mapping builder
  - Read CG system (coordinates + topology)
  - Read AT source files per molecule type
  - Build hybrid system: place AT atoms around CG positions using bead definitions from YAML
  - Classify bonds into three categories (intra-CG static, cross-CG AT, cross-CG CG)
  - Assign LAMMPS atom types and bond types

- [x] **E6 [MVP]**: Implement LAMMPS data file writer
  - Write `atom_style full` data file
  - Header, box, masses, atoms (CG first per molecule), bonds, angles
  - Distinct bond/angle types for static vs `backmap/*` styles
  - Print type mapping tables
  - Spec: `backmap-input-generator/spec.md` → "LAMMPS data file generation"

- [x] **E7 [MVP]**: Implement tabulated potential converter
  - Read GROMACS `.xvg` tables
  - Convert to LAMMPS `.table` format with unit conversion
  - Support pair, bond, and dihedral table types
  - Spec: `backmap-input-generator/spec.md` → "Tabulated potential conversion"

- [x] **E8 [MVP]**: Implement LAMMPS input script writer
  - Generate complete `.in` file from simulation parameters
  - `pair_style backmap` with correct sub-styles
  - `bond_style hybrid` routing (static vs `backmap/*`)
  - Group definitions, fix backmap, fix nve, fix langevin
  - `special_bonds` from exclusion settings
  - Three-phase run sequence
  - Spec: `backmap-input-generator/spec.md` → "LAMMPS input script generation"

- [x] **E9 [P2]**: ~~Add cross-dihedral and 1-4 pair support~~ → Moved to `backmap-phase2-and-testing`

- [x] **E10 [P3]**: ~~Add reactive network support~~ → Moved to `backmap-phase2-and-testing`

- [x] **E11 [P4]**: ~~Add non-GROMACS source format support~~ → Moved to `backmap-phase2-and-testing`

---

### Group F: Integration Testing and Validation

Depends on all MVP tasks from Groups A-E.

- [x] **F1 [MVP]**: ~~Water2 system — single-bead water, no cross bonds~~ → Moved to `backmap-phase2-and-testing`

- [x] **F2 [MVP]**: Dodecane system — linear molecule with cross bonds
  - Generate YAML settings for dodecane (6 CG beads, 12 AT atoms per molecule)
  - Verify: cross-CG AT bonds weighted by λ², CG bonds by (1−λ²), intra-CG bonds static
  - Verify: correct bond type routing in data file and input script

- [x] **F3 [MVP]**: ~~PE (polyethylene) system — long chain with tabulated CG potentials~~ → Moved to `backmap-phase2-and-testing`

- [x] **F4 [MVP]**: ~~MPI parallel correctness~~ → Moved to `backmap-phase2-and-testing`

- [x] **F5 [P3]**: ~~Epoxy network (RIM135) — reactive multi-molecule system~~ → Moved to `backmap-phase2-and-testing`

---

### Dependency Graph

```
A1 ─────┬──── B1 → B2 → B3 → B4 → B5 → B6 → B7
        │
        ├──── C1 → C2 → C3 → C4 → C5
        │
        └──── D1, D2, D3 (parallel)
                              │
A2 ──────────────────────────┘

E1 → E2 → E3 → E4 → E5 → E6 → E7 → E8

All MVP ──→ F1 → F2 → F3 → F4
```

Groups B, C, D can be developed in parallel (all depend only on A1).
Group E (Python) can be developed in parallel with all C++ groups.
Group F requires all MVP tasks complete.

---

### Estimated Scope

| Group | MVP Tasks | Est. Lines (C++) | Est. Lines (Python) |
|-------|-----------|-------------------|---------------------|
| A     | 2         | ~80               | —                   |
| B     | 7         | ~600              | —                   |
| C     | 5         | ~400              | —                   |
| D     | 3 (MVP)   | ~450              | —                   |
| E     | 8 (MVP)   | —                 | ~1200               |
| F     | 4         | —                 | ~200 (test scripts) |
| **Total MVP** | **29** | **~1530**    | **~1400**           |
