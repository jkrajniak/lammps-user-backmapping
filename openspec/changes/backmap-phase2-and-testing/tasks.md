## Tasks

Tasks carried over from `lammps-backmapping-package` (MVP change). These are the remaining items not completed in Phase 1.

### Legend

- **[MVP]** = Phase 1, required for initial working system (integration tests)
- **[P2]** = Phase 2, extended bonded interactions
- **[P3]** = Phase 3, reactive networks
- **[P4]** = Phase 4, format flexibility

---

### Group B: `fix backmap` Extensions (C++)

- [ ] **B8 [P2]**: Implement two-phase backmapping support
  - Phase 1: bonded + unweighted CG, phase 2: reset λ, add AT
  - `fix_modify` to switch phases
  - Spec: `fix-backmap-resolution/spec.md` → "Simulation phases" (two-phase)

---

### Group D: Dihedral `backmap/*` Styles (C++)

- [ ] **D4 [P2]**: Implement `dihedral_backmap_ryckaert.cpp/.h`
  - Parse: `dihedral_coeff N at/cg C0 C1 C2 C3 C4 C5`
  - Weight from first and last atoms of the dihedral quadruplet
  - Spec: `bonded-backmap-styles/spec.md` → "Dihedral style backmap/ryckaert"

- [ ] **D5 [P2]**: Implement `dihedral_backmap_table.cpp/.h`
  - Parse: `dihedral_coeff M at/cg filename keyword`
  - Follow LAMMPS `dihedral_style table` conventions
  - Spec: `bonded-backmap-styles/spec.md` → "Dihedral style backmap/table"

---

### Group E: Python Input Generator Extensions (`backmap-prep`)

- [ ] **E9 [P2]**: Add cross-dihedral and 1-4 pair support
  - Extend data file writer for dihedrals section
  - Extend input script writer for `dihedral_style hybrid backmap/*`
  - Spec: `backmap-input-generator/spec.md` → "Cross-interactions section" (dihedrals)

- [ ] **E10 [P3]**: Add reactive network support
  - Degree-dependent bead definitions
  - Active site logic and bond formation
  - Charge management (equilibrate, transfers)
  - Atom removal, type/charge maps
  - Multiple molecule types
  - Spec: `backmap-input-generator/spec.md` → "Molecules section" (deferred features)

- [ ] **E11 [P4]**: Add non-GROMACS source format support
  - PDB, LAMMPS data, XYZ parsers
  - LAMMPS `.table` passthrough
  - Spec: `backmap-input-generator/spec.md` → "Source file parsing" (Phase 4)

---

### Group F: Integration Testing and Validation

- [ ] **F1 [MVP]**: Water2 system — single-bead water, no cross bonds
  - Generate YAML settings for SPC water + WCG (3456 molecules)
  - Run `backmap-prep` to produce LAMMPS files
  - Run LAMMPS backmapping simulation
  - Verify: λ ramp, COM tracking, force distribution, final AT structure
  - Compare against ESPResSo++ water2 results

- [ ] **F3 [MVP]**: PE (polyethylene) system — long chain with tabulated CG potentials
  - Generate YAML settings for PE (50 CG beads, 100 AT atoms)
  - Verify: tabulated CG bond/pair potentials, harmonic AT cross bonds
  - Verify: full backmapping 0→1 with correct force weighting

- [ ] **F4 [MVP]**: MPI parallel correctness
  - Run water2 and PE on 4+ processors
  - Verify: molecule mapping rebuilt after migration, ghost atom λ communication, identical results to serial

- [ ] **F5 [P3]**: Epoxy network (RIM135) — reactive multi-molecule system
  - Validate degree-dependent beads, active sites, charge transfers
  - Compare against bakery/ESPResSo++ results

---

### Dependency Graph

```
B8 (standalone — extends existing fix backmap)

D4, D5 (parallel — depend on existing backmap_lambda.h)

E9 depends on D4/D5
E10 (standalone — Phase 3)
E11 (standalone — Phase 4)

F1, F3 (parallel — need compiled LAMMPS + test data)
F4 depends on F1/F3
F5 depends on E10
```
