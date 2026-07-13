# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/).

## [Unreleased]

### Fixed

- **AT intra-bead lambda weighting**: `compute_weight3` incorrectly gated
  AT intra-bead (same-CG-bead) terms to `w = 0` at `lambda_global = 0`
  instead of always `w = 1`. Intra-molecular AT chemistry (bonds, angles,
  dihedrals, LJ within one fragment) exists independent of the CG/AT
  resolution ramp, matching ESPResSo++'s production driver
  (`start_backmapping.py`), which builds the full AT interaction list
  unconditionally before the ramp mechanism is ever activated. Fixed in
  `backmap_lambda_weights.h`; `pair_backmap`'s `LAMBDA_AT_ONSET` deferral
  now only applies to inter-bead (intermolecular) AT-AT pairs. Validated
  by 21 unit tests and a VM regression across dodecane/PE/rim135/melamine
  (all Tier B PASS after also adding a Phase 0b Langevin thermostat to
  each example, since AT strain now releases as real, un-damped kinetic
  energy during that stage).

- **Lambda-weighting formula**: all 8 `backmap/*` styles (`pair_backmap`,
  `bond_backmap_harmonic`/`table`, `angle_backmap_harmonic`/`table`,
  `dihedral_backmap_ryckaert`/`table`, `fix_backmap/pairs`) previously
  weighted interactions by the per-particle product `λ_i × λ_j` (AT) /
  `1 − λ_i × λ_j` (CG). Replaced with a single global λ scalar
  (`FixBackmap::lambda_global`, exposed via `extract("lambda_global", ...)`)
  combined with same-CG-bead membership (`FixBackmap::atom2cg`, exposed via
  `extract("atom2cg", ...)`): CG → `1 − λ_global`; AT atoms in the *same* CG
  bead → full strength once `λ_global > 0` (not scaled); AT atoms in
  *different* CG beads → `λ_global` (linear). See
  `backmap_lambda_weights.h::compute_weight3()`/`same_bead()`.

### Added

- **Fast C++ unit tests** (`tests/unit/`, GoogleTest via CMake `FetchContent`,
  no LAMMPS build required): 21 tests covering the lambda-weighting formula
  and same-bead detection. Wired into `.github/workflows/cpp-ci.yml`
  (`unit-tests` job) and `make test-cpp`.
- **`.clang-format`**: explicit `BasedOnStyle: Google` +
  `PointerAlignment: Right` config so the pre-commit/CI `--style=file` hook
  no longer relies on an undeclared fallback.
- **MPI-correct `fix backmap`**: COM tracking and CG force distribution now use
  atom-centric communication (`reverse_comm` for COM accumulation,
  `forward_comm` for CG force) following the `fix_rigid_small` pattern, so
  domain-decomposed runs no longer drop AT atoms on non-neighbour ranks or lose
  force writes to ghost atoms. Requires `comm_modify cutoff` >= maximum CG-AT
  distance within a bead (a few angstrom for backmapped fragments); the fix
  warns if a local AT atom has no CG partner in local+ghost range.
- **`compute_scalar`** on `fix backmap`: group-averaged lambda, printable as
  `f_bm` in `thermo_style` (with `thermo_modify colname f_bm lambda`).
- **MPI parity test**: `examples/dodecane/large/in.dodecane_mpi` +
  `test_mpi_serial_vs_4rank.sh` + `compare_mpi_data.py`.

### Fixed

- `fix backmap` under MPI: incomplete COM sum (AT atoms on non-ghost ranks
  missed) and lost CG-force contributions to ghost AT atoms.


## [0.1.0] - 2026-07-06

First tagged release: network backmapping engine (rim135 epoxy), Tier B dynamics
protocol, and structural validation vs JCC 2017/2018 paper reference.

### Added

- **Phase 3 network engine** (`backmap_prep.network`): Settings v2 YAML loader,
  hybrid GROMACS builder ported from bakery, LAMMPS unified `build` / `build-hybrid` /
  `finalize-cg` / `rebuild` CLI paths.
- **`examples/epoxy/`** rim135 example: `settings.v2.yaml`, `run_test.sh`, bundled
  OPLS-AA forcefield slice, Tier A parity tests.
- **`compare_rim135_structure.py`**: C–O / C–N RDF validation vs GROMACS `.xvg`
  references in `paper-reverse-mapping-polymer-networks`; pinned report 4/4 peak
  metrics PASS (Jul 2026).
- CG angle/dihedral/cross-pair table export; `fix backmap/pairs` C++ style;
  molecule-aware PBC export (`network/pbc.py`) with network image flags.
- PR4 Tier B bakery protocol: `cap_force`, Langevin `gamma`, velocity init,
  `comm_modify cutoff` auto from bonded extent for network hybrids.

### Fixed

- Network hybrid communication cutoff: large Cartesian bond extent requires
  `comm_modify cutoff` ≥ bonded ghost estimate (~115 Å for rim135); prevents
  missing bond atoms during λ ramp on VM.

### Changed

- `examples/epoxy/README.md`: documents supported `build` vs experimental `rebuild`
  paths.

## [Unreleased — pre-0.1.0 history retained below]

### Added

- Paper-grade RDF validation for the 500-molecule dodecane example:
  `examples/dodecane/large/in.dodecane_at_long` (post-backmap AT, multi-phase
  equilibration + 1 ns NVT production) and `in.dodecane_at_ref_long` (independent
  all-atom reference at the same density / temperature), both averaging g(r) in
  five 200 ps blocks via `fix ave/time 100 2000 200000`.
- `examples/dodecane/large/compare_rdf_blocks.py` post-processor that parses
  multi-block `fix ave/time` output, computes per-pair mean ± SEM across blocks,
  plots SEM-shaded bands, and reports PASS/FAIL on first-peak position, height,
  and L2(g_bm − g_ref) tolerances.
- Committed reference outputs `examples/dodecane/large/rdf_comparison_long.{png,txt}`
  documenting the validation: 9 / 9 metrics pass with first-peak positions
  matching to the bin width (≤ 0.07 Å), heights agreeing within 0.5 %, and
  L2(g_bm − g_ref) ≤ 0.003 across all three pairs.
- `apb` keyword in `fix backmap` for non-uniform atoms-per-bead mapping
  (`apb T1:N1 T2:N2 ...`). Required for systems where different CG bead
  types contain different numbers of AT atoms (e.g. all-atom PE with
  7-atom end beads and 6-atom interior beads).

### Fixed

- `fix backmap` multi-`run` segfault on the 500-molecule dodecane example
  (`Verlet::setup()` → `Domain::box_too_small_check()` reading a NULL `x[k]`
  on the second `run` after `unfix`/`fix` changes). Resolved by registering
  `atom->add_callback(Atom::GROW)` and `atom->add_callback(Atom::RESTART)`
  in `fix backmap`, which now triggers `grow_arrays()` whenever LAMMPS
  reallocates atom storage and prevents the heap-use-after-free of the
  per-atom `lambda` array. Verified end-to-end on
  `examples/dodecane/large/in.dodecane` (50 000 steps across all 3 phases,
  final T = 294 K, no crashes / mass mismatches / bond-missing errors).
- All PE example input scripts now use `cg_type 1 2` (both CG bead types)
  instead of `cg_type 1`, which caused incorrect bead-to-atom mapping.
- Melamine large example: `topol_cg.top` molecule count corrected from 50 to
  500; regenerated `melamine.data` (15 000 atoms).
- Melamine input scripts (`in.melamine`) use the robust multi-phase protocol
  instead of aggressive NVE + Langevin at dt = 1.0 fs.

- `fix backmap` now accepts multiple CG atom types via `cg_type T1 T2 ...`
  syntax, enabling correct bead-to-atom mapping in systems with more than one
  CG bead type.
- `fix backmap` COM tracking uses round-based PBC wrapping, correctly handling
  ghost atoms that may be more than one box length from the local CG bead.
- `fix backmap` bead map is rebuilt on every neighbor list rebuild (via
  `pre_force()` callback), preventing stale local indices after LAMMPS atom
  sorting.
- `pair_style backmap` defers AT pair interactions until both atoms reach
  &lambda; > 0.1, preventing LJ singularities from initial inter-molecular
  overlaps.
- Python writer emits `cg_type T1 T2 ...` with all CG type IDs and wraps
  atom coordinates into [0, L) in the LAMMPS data file.

### Changed

- Large-scale example input scripts use a robust multi-phase simulation
  protocol (minimise → `nve/limit` relaxation → `nve/limit` lambda ramp →
  gradual NVT equilibration) instead of aggressive single-phase NVE + Langevin,
  preventing "Bond atoms missing" errors in production-size systems.
- Documentation: README, `docs/index.md`, `docs/theory.md`, and `AGENTS.md` now
  tie the package motivation to **migrating from ESPResSo++** toward **LAMMPS**,
  and cite the **complex polymer / network** reverse-mapping paper with
  **December 2017** online publication date (DOI 10.1002/jcc.25129; print *J.
  Comput. Chem.* 2018).

- **backmap-prep** LAMMPS setup: emit `neigh_modify delay 0 every 1 check yes`
  after group definitions; use the **backmapping timestep** for the first
  dynamics segment (including λ-frozen equilibration when
  `equilibration_steps` > 0) instead of the larger production timestep, so
  hybrid relaxation does not start at an overly aggressive Δt.
- **examples/dodecane/n250** `settings.yaml`: `timestep_backmapping: 0.00025` ps
  and `equilibration_steps: 4000` so the default generated input completes the
  λ ramp without early missing-bond errors on the 250-molecule melt.
- **examples/dodecane/n250** RDF validation: `in.dodecane_n250_at`,
  `in.dodecane_n250_at_ref`, and `compare_rdf_n250.sh` (runs
  `../large/compare_rdf.py` strict + relaxed tolerances for small-N noise).

- **backmap-prep** default `simulation.production_steps` is now `0`: generated
  backmapping inputs end with `write_data …_hybrid.data` after the λ ramp only;
  atomistic production (e.g. RDF) is intended as a separate run on the extracted
  AT system. Set `production_steps` > 0 only to append an optional post-backmap
  segment in the same file.
- `scripts/run-backmap.sh`: supports single-segment restart layouts (only
  `*.phase1`, backmapping-only) in addition to two- and three-segment flows.
- Packaging and docs: Python project metadata (`pyproject.toml`, `uv.lock`) is now at repository root so `uv sync` from root installs `backmap-prep`; documentation and examples now use root-level `uv run backmap-prep ...` commands.
- README: repository layout section now lists all example directories and the `scripts/` validation script.

### Added

- **250-molecule dodecane** example layout: `examples/dodecane/n250/` with
  `build_cg_conf.py` (subset of `large/cg_conf.gro`), matching `topol_*_250.top`
  files, `settings.yaml`, and `prepare_inputs.sh`. Suited to laptops (hybrid
  ~4500 atoms; serial LAMMPS RSS typically low compared to 32 GB RAM).

- Small dodecane example: `examples/dodecane/in.dodecane_at` and
  `extract_at_frame.py` for the hybrid → AT → RDF workflow; README documents
  the `compare_rdf.py` sequence.

- Restart/checkpoint support for preemptible cloud instances: new
  `restart_interval` setting in `simulation` generates `write_restart`
  commands and per-phase input scripts. Includes `scripts/run-backmap.sh`
  entrypoint that detects restart files and resumes from the correct phase.
- Google Cloud Batch example (`examples/cloud-batch/`) with spot VM job
  template and setup instructions.
- "Running on Cloud / HPC" documentation page (`docs/cloud-hpc.md`) covering
  restart configuration, Cloud Batch, and Slurm/Apptainer workflows.
- Dockerfile for building LAMMPS with the backmapping package in a container
  (multi-stage build, configurable LAMMPS version via `LAMMPS_VERSION` build
  arg, includes `backmap-prep` CLI). Supports conversion to
  Singularity/Apptainer for HPC clusters.
- Docker documentation page (`docs/docker.md`) with build, run, MPI, and HPC
  conversion instructions.

- Large-scale example variants in `examples/<name>/large/` for dodecane, pe, pe4, pe_10, pe_aa, and melamine. Inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project; each `large/` contains a README and `settings.yaml` so `backmap-prep` can generate LAMMPS data and input files. See [Large-scale examples](docs/large-scale-examples.md) and the main README.

- C++ LAMMPS styles: `fix backmap`, `pair_style backmap`,
  `bond_style backmap/harmonic`, `bond_style backmap/table`,
  `angle_style backmap/harmonic` for time-dependent CG-to-AT backmapping.
- CMake build integration (`src/CMakeLists.txt`) and legacy `Install.sh`.
- Python CLI tool `backmap-prep` that generates LAMMPS data files, input
  scripts, and interaction tables from a YAML settings file.
- Pydantic v2 settings schema with validation and deferred-feature guards.
- GROMACS topology and coordinate parsers.
- LAMMPS data and input file writers.
- XVG-to-LAMMPS table converter.
- Dodecane example (`examples/dodecane/`) demonstrating full backmapping
  workflow with tabulated CG interactions.
- Makefile with convenience targets for install, lint, format, typecheck,
  test, and pre-commit.
- Pre-commit configuration with ruff and mypy hooks.
- Project constitution and OpenSpec change tracking (`openspec/`).
- Documentation requirements: `CHANGELOG.md` and `README.md` must be kept
  up to date with every change.
