# Tasks: unify the hybrid engine

Design: [design.md](./design.md). Branch `feat/single-hybrid-engine`.

## 1. Reference and guard rails
- [ ] 1.1 GROMACS 2023 on the compute VM (user space, conda-forge)
- [ ] 1.2 `tests/reference/gmx_energy_check.py`: AT-only GROMACS rerun vs LAMMPS
      `run 0` at lambda = 1, per-term comparison (VM-only; documented)
- [ ] 1.3 Property tests: bead-on-fragment-COM, CG bonded completeness
- [ ] 1.4 Freeze current network outputs (rim135, PET, melamine_network) as
      byte-level baselines before touching the engine

## 2. Engine
- [x] 2.1 `v2_loader._prefix_atoms`: keep already-qualified names
- [ ] 2.2 Schema: remove `prep.engine` (deprecation warning), `hybrid` optional
      with defaults
- [ ] 2.3 United-atom LJ from `[ atomtypes ]` sigma/epsilon in lammps_builder
- [ ] 2.4 Communication cutoff: min-image based, not Cartesian extent
- [ ] 2.5 LAMMPS-format `cg_system` and `molecules[].source` -> temporary
      `.gro`/`.top`, exact unit round trip
- [ ] 2.6 CLI: `build`, `rebuild`, `cg-only` through the single engine
- [ ] 2.7 Delete `builder.build_system` and linear-only helpers; keep the
      data classes (`System`, `LammpsAtom`, ...) in a neutral module

## 3. Examples
- [ ] 3.1 Regenerate every example; COM and completeness checks pass
- [ ] 3.2 GROMACS energy check passes for dodecane, PE, pe4, pe_10, pe_aa,
      melamine, pe-lammps, dodecane-lammps-cg
- [ ] 3.3 Production inputs take their force-field block from the generated
      input (no restated coefficients)
- [ ] 3.4 Network baselines (1.4) unchanged

## 4. Tests, docs
- [ ] 4.1 Port `test_builder.py` cases to the single engine; delete
      linear-engine-only assertions that encoded the defects
- [ ] 4.2 Settings reference, tutorial, getting started, README, CHANGELOG
- [ ] 4.3 `make lint typecheck test docs`
