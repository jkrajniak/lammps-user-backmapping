## 1. Project Setup

- [x] 1.1 Create `docs/` directory structure matching the design layout
- [x] 1.2 Create `docs/requirements.txt` with `mkdocs` and `mkdocs-material` dependencies
- [x] 1.3 Create `mkdocs.yml` with site name, theme config, navigation, and strict mode enabled
- [x] 1.4 Add `docs` and `docs-serve` targets to the `Makefile`

## 2. Landing Page and Getting Started

- [x] 2.1 Write `docs/index.md` — project overview, feature highlights, and navigation pointers
- [x] 2.2 Write `docs/getting-started.md` — installation of C++ package and Python CLI, running the dodecane example

## 3. Theory and Background

- [x] 3.1 Write `docs/theory.md` — lambda parameter, alpha ramp rate, force weighting (bonded and non-bonded), cross-interaction concept

## 4. Settings Reference

- [x] 4.1 Write `docs/settings-reference.md` — document all `molecules` fields (name, ident, source, beads, atoms)
- [x] 4.2 Document `cg_system` fields (coordinates, topology, format)
- [x] 4.3 Document `cross_interactions` fields (bonds, angles, dihedrals — params, pairs/triples, table, cg_bonded)
- [x] 4.4 Document `simulation` fields (alpha, timestep, temperature, thermostat, cutoffs, etc.) with types and defaults
- [x] 4.5 Document `output` fields (prefix, format, units)

## 5. Tutorial

- [x] 5.1 Write `docs/tutorial-new-system.md` — overview and prerequisites
- [x] 5.2 Write tutorial section: preparing CG and AT input files (GROMACS format)
- [x] 5.3 Write tutorial section: writing the settings YAML step by step
- [x] 5.4 Write tutorial section: running `backmap-prep` and inspecting generated files
- [x] 5.5 Write tutorial section: executing the LAMMPS backmapping simulation and checking output

## 6. Component Documentation

- [x] 6.1 Write `docs/components/fix-backmap.md` — command syntax, arguments, usage examples
- [x] 6.2 Write `docs/components/pair-backmap.md` — pair_style syntax, coefficients, lambda weighting
- [x] 6.3 Write `docs/components/bond-styles.md` — backmap/harmonic and backmap/table syntax, coefficients
- [x] 6.4 Write `docs/components/angle-styles.md` — backmap/harmonic angle syntax, coefficients

## 7. CLI Reference

- [x] 7.1 Write `docs/cli/backmap-prep.md` — command syntax, flags, example invocations

## 8. CI/CD and Deployment

- [x] 8.1 Create `.github/workflows/docs.yml` — build and deploy docs on push to main
- [x] 8.2 Update `README.md` with link to published documentation site

## 9. Verification

- [x] 9.1 Run `mkdocs build --strict` locally and verify zero warnings/errors
- [x] 9.2 Run `mkdocs serve` and manually verify all pages render correctly with working navigation
