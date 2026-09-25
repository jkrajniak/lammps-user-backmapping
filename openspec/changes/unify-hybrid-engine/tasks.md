# Tasks: unify the hybrid engine

Design: [design.md](./design.md). Branch `feat/single-hybrid-engine`.

## 1. Reference and guard rails
- [x] 1.1 GROMACS 2023 on the compute VM (user space, conda-forge)
- [x] 1.2 `tests/reference/gmx_energy_check.py`: AT-only GROMACS rerun vs LAMMPS
      `run 0` at lambda = 1, per-term comparison (VM-only; documented)
- [x] 1.3 Property tests: bead-on-fragment-COM, CG bonded completeness
- [ ] 1.4 Freeze current network outputs (rim135, PET, melamine_network) as
      byte-level baselines before touching the engine

## 2. Engine
- [x] 2.1 `v2_loader._prefix_atoms`: keep already-qualified names
- [x] 2.2 Schema: remove `prep.engine` (deprecation warning), `hybrid` optional
      with defaults
- [x] 2.3 United-atom LJ from `[ atomtypes ]` sigma/epsilon in lammps_builder
- [x] 2.4 Communication cutoff: min-image based, not Cartesian extent
- [x] 2.5 LAMMPS-format `cg_system` and `molecules[].source` -> temporary
      `.gro`/`.top`, exact unit round trip
- [x] 2.6 CLI: `build`, `rebuild`, `cg-only` through the single engine
- [x] 2.8 Write outputs to the settings directory (or an explicit output
      directory), never into `prep.data_dir`: a PET build wrote into the
      published paper-data repository
- [x] 2.7 Delete `builder.build_system` and linear-only helpers; keep the
      data classes (`System`, `LammpsAtom`, ...) in a neutral module

## 2b. Found while porting (done)
- [x] Exact kJ->kcal (1/4.184, was 0.239006) and 10-digit coefficients
- [x] LJ mixing follows the topology combination rule (was always LB)
- [x] Clear error when a CG residue has no molecule definition
- [x] CG pair tables as LAMMPS `.table` files (passthrough)
- [x] Examples restored to the published CG models (dodecane, pe4, pe_aa,
      dodecane-lammps-cg); melamine CG residues per molecule

## 3. Examples
- [ ] 3.1 Regenerate every example; COM and completeness checks pass
- [x] 3.2 GROMACS energy check passes for dodecane, PE, pe4, pe_10, pe_aa,
      melamine, pe-lammps, dodecane-lammps-cg
- [ ] 3.3 Production inputs take their force-field block from the generated
      input (no restated coefficients)
- [ ] 3.4 Network baselines (1.4) unchanged

## 4. Tests, docs
- [ ] 4.1 Port `test_builder.py` cases to the single engine; delete
      linear-engine-only assertions that encoded the defects
- [ ] 4.2 Settings reference, tutorial, getting started, README, CHANGELOG
- [ ] 4.3 `make lint typecheck test docs`
