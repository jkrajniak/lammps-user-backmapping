## Proposal: Backmapping Phase 2+ and Remaining Integration Tests

### Context

The MVP implementation of the `USER-BACKMAP` LAMMPS package and `backmap-prep` Python generator was completed in the `lammps-backmapping-package` change. This follow-up change captures the remaining work: Phase 2 extended bonded interactions, Phase 3/4 Python generator features, and MVP/Phase 3 integration tests.

### Problem

The MVP covers harmonic bonds, harmonic angles, tabulated bonds, and the pair style — but the full backmapping workflow for complex molecules also requires:
- Dihedral interactions (Ryckaert-Bellemans and tabulated)
- Two-phase backmapping protocol
- Reactive network support for epoxy-like systems
- Integration tests proving correctness on real systems (water2, PE, MPI, epoxy)

### Proposed Solution

Complete the remaining 10 tasks carried over from the original change, organized into three areas:

1. **Phase 2 C++ styles** — dihedral styles and two-phase backmapping in `fix backmap`
2. **Phase 2-4 Python generator** — dihedral support, reactive networks, non-GROMACS formats
3. **Integration testing** — water2, PE, MPI parallel, and epoxy validation

### Non-Goals

- Reimplementing any MVP features (those are done and tested)
- Changing the YAML settings schema (stable from MVP)
- Spatial AdResS regions (out of scope for backmapping)
