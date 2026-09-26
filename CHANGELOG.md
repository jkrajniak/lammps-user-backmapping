# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Conventional Commits](https://www.conventionalcommits.org/).

## [Unreleased]

### Added

- **`fix backmap ... peratom full`**: per-atom array with lambda, the bead
  position (COM) and the CG force share each AT atom received (beads: their
  CG force), for static decomposition-parity tests. Default output unchanged.

- **`pair_style backmap ... cg_special w12 w13 w14`** and special-bond
  factors in general. The pair style evaluated every listed pair at factor 1,
  correct only because generated inputs exclude 1-2..1-4 fully; it now passes
  the `special_bonds` factors to its sub-styles, and `cg_special` sets
  separate weights for CG-CG pairs (bakery's `exclusion_cg`), e.g. a MARTINI
  CG model at nrexcl = 1 with an AT force field at nrexcl = 3. Results of the
  existing examples (`special_bonds` 0 0 0) are unchanged.

- **CHARMM-type bonded styles.** `angle_style backmap/charmm` (harmonic +
  Urey-Bradley), `dihedral_style backmap/fourier` (multi-term periodic with
  arbitrary phase) and `improper_style backmap/harmonic`, with the usual
  lambda weighting, restart and `write_data`. Kernels are the stock LAMMPS
  `angle charmm`, `dihedral fourier` and `improper harmonic`; tests check equality
  with them times the weight (energies and forces). Needed for CHARMM36/Slipids
  lipids (MARTINI 3 POPC example) and requested in review (impropers).

- **`cg_system.format: lammps`**: the CG side of `backmap-prep` can now be
  supplied as a native LAMMPS `data` file (box, `Masses`, `Atoms # full`)
  instead of GROMACS `.gro`/`.top`. No unit conversion is applied (the file
  is assumed to already be in `units real`); CG-CG bonded terms
  (`cross_interactions`) and nonbonded tables (`table_groups`) are
  configured identically to the GROMACS path. The CG pair-table lookup also
  now prefers a `.table` file over `.xvg` when both are present. See
  `examples/dodecane-lammps-cg/` and [Settings Reference: `cg_system`](https://jkrajniak.github.io/lammps-user-backmapping/settings-reference/#cg_system).
- **`molecules[].source.format: lammps`**: each molecule's AT fragment can
  now also be supplied as a native LAMMPS `data` file plus a bounded
  input-script fragment (`bond_coeff`/`angle_coeff`/`dihedral_coeff`/
  `pair_coeff`), instead of GROMACS `.gro`/`.top`. `beads[].atoms` then
  reference the LAMMPS numeric atom ID (as a string) rather than a
  symbolic atom name. Combined with `cg_system.format: lammps`, this makes
  `backmap-prep` fully GROMACS-free. See `examples/pe-lammps/` and
  [Settings Reference: `molecules[].source`](https://jkrajniak.github.io/lammps-user-backmapping/settings-reference/#moleculessource).
  Not supported for degree-dependent (Phase 3) sources, the network engine,
  or GROMACS-virtual-site AT fragments.

### Fixed

- **`fix backmap` COM update no longer depends on fix order.** It ran in
  `initial_integrate()`, so with `fix backmap` defined before the
  integration fix (as in the generated robust protocol, which redefines its
  integrators at every stage) the beads followed the AT positions of the
  previous step. It now runs in `post_integrate()`.

- **`fix backmap` distributes the CG forces at setup.** `setup()` did not
  call `post_force()`, so after the force evaluation that opens every run
  (and `run 0`) the CG forces stayed on the beads and the AT atoms got none;
  the first velocity half-kick of each `run` missed them. Affects only the
  first step of runs where CG forces are active (e.g. the start of a ramp).

- **`fix backmap/pairs` has no cutoff by default.** 1-4 pairs are bonded
  terms, but the fix skipped any pair beyond the pair-style cutoff (12 Å in
  the examples). In freshly built network frames 77 of 35612 RIM135 1-4
  pairs are longer than that, so their LJ and Coulomb terms were dropped
  (Coul-14 off by 5% against GROMACS, which applies no cutoff). `cut` is now
  opt-in. A pair whose partner is not present on the owning rank is an error
  instead of being skipped silently.
- **`fix backmap/pairs`: 1-4 Coulomb, energy output, setup, minimization and
  MPI.** The fix applied 1-4 LJ only, so with `special_bonds coul 0 0 0`
  the scaled 1-4 electrostatics (GROMACS `fudgeQQ`) were absent. It now
  takes an optional fifth pairs-file column, the 1-4 Coulomb scale. It
  reported no energy or virial; it now contributes both (scalar: total,
  vector: LJ-14, Coulomb-14), counted in `pe` and pressure by default. It had
  no `setup()` or `min_post_force()`, so its forces were missing from the
  first force evaluation of every run (including `run 0`) and from all
  minimizations. With `newton_pair on` it added forces to ghost atoms after
  the reverse communication, so the force on a partner owned by another
  rank was lost; every rank now evaluates each pair it owns an atom of and
  applies the force to its own atoms. Regression tests:
  `python/tests/test_lammps_pairs.py` (needs `BACKMAP_LMP`).

- **`write_data` left empty coefficient sections** for `bond_style
  backmap/harmonic`, `angle_style backmap/harmonic`, `dihedral_style
  backmap/harmonic`, `dihedral_style backmap/ryckaert` and `dihedral_style
  ryckaert`, which then made `read_data` fail ("Unexpected empty line in
  AngleCoeffs section"). These styles now implement `write_data()`.
  `backmap/table` bond and angle styles no longer write a section at all,
  as for stock `bond_style table`. Regression tests:
  `python/tests/test_lammps_write_data.py` (needs `BACKMAP_LMP`).
- **`pair_style backmap` reported twice the pair energy**: `PairBackmap::compute()`
  added each weighted pair energy to `eng_vdwl`/`eatom` by hand and then again
  through `ev_tally()`. Forces and the virial were always correct, so
  trajectories, structures and pressures are unaffected; every reported pair
  energy (`evdwl`, `pe`, `etotal`, `compute pe/atom`) was inflated by the pair
  term. Regression test: `python/tests/test_lammps_energy.py` (needs
  `BACKMAP_LMP`).

- **`build_network_lammps` under-sized `comm_modify cutoff` for crosslinked
  networks**: `build_system_from_hybrid` (the code path behind `backmap-prep
  build` for `engine: network` systems) already computed correct per-atom
  image flags via `prepare_network_coordinates()`, but never set
  `system.write_image_flags = True` -- the flag `writers.py::_compute_params`
  gates its network-aware `comm_modify cutoff` widening on. Every crosslinked
  network built via `build` (not `rebuild`/`finalize-cg`, which already set
  the flag) silently got a `comm_modify cutoff` sized only for the LJ/CG
  interaction cutoff (15 Å default) instead of the folded Cartesian extent of
  box-spanning crosslink bonds (82 Å for rim135), which is why every rim135
  and PET/Dacron production script in this repo has had to override
  `comm_modify cutoff` by hand. Fixed in `network/lammps_builder.py`
  (`build_system_from_hybrid`); regenerating `examples/epoxy/large/in.rim135`
  now writes `comm_modify cutoff 82.23` automatically, matching the
  hand-tuned value already documented in the manuscript's MPI section.
  `test_rim135_build_v2_lammps_smoke` (previously the suite's one failing
  test) now passes; full suite 207/207.

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

### Removed

- **Dead per-atom lambda, `nonuniform`, and `phase` machinery**, left over
  from `fix backmap`'s original broader AdResS-style design and made
  fully inert by the lambda-weighting fix above: `lambda_global` was
  already the single source of truth for every interaction style, so
  none of this had any effect on any validated example.
  - **BREAKING (LAMMPS input-script syntax)**: `fix backmap` no longer
    accepts `nonuniform yes/no`. Removed the per-atom `lambda[]` array's
    independent state (`RanMars`-based staggering, `pack_exchange`/
    `unpack_exchange`, restart persistence) -- `nonuniform` was `no` in
    every checked-in example, so the per-atom array was always
    numerically identical to `lambda_global`. The array survives, renamed
    `lambda_display`, as a pure broadcast mirror of `lambda_global`
    (refreshed unconditionally every `end_of_step()`) so `f_bm`
    dump/thermo output in existing examples (`in.pe`, `in.pe4`,
    `in.pe_robust`) is unaffected. `comm_forward` shrinks from 5 to 4
    doubles/ghost-atom (dropped the never-read ghost lambda slot).
    `compute_scalar()` simplifies to `return lambda_global;` --
    provably equal to the old per-atom mean for every real run, since
    `nonuniform=no` means every atom's value was always `lambda_global`
    already.
  - Also removed the `phase` fix-keyword: reserved in the argument
    parser but with no handling branch at all (a user specifying it got
    a hard parse error, so nothing that worked before still works
    differently now) and `compute_weight3()`'s `phase` parameter, never
    passed a non-default value by any of the 14 production call sites.
  - Corrected `docs/components/fix-backmap.md`'s Restart section, which
    claimed "seamless continuation" -- already false before this change
    (`lambda_global` was never restart-persisted; only the now-removed
    per-atom `lambda[]` was).
  - Rewrote the `pair-backmap` and `bonded-backmap-styles` openspec specs'
    weighting requirements, which still described the pre-lambda-weighting-fix
    `λ_i × λ_j` product model; removed `fix-backmap-resolution`'s
    "Simulation phases" requirement (never implemented).
  - See `openspec/changes/remove-dead-lambda-nonuniform-phase/` for the
    full verification trail.

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
