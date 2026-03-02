## Context

The `lammps-user-backmapping` project has a single working example (dodecane) demonstrating the YAML-driven `backmap-prep` workflow. The bakery project contains 7+ additional examples of varying complexity that were built for ESPResSo++. We need to migrate the subset that is compatible with the current Phase 1/2 capabilities of `backmap-prep` (linear polymers, tabulated bonds/angles/dihedrals).

The bakery examples use XML settings files, GROMACS-format input files (GRO, TOP, XVG tables), and the `prepare_files.py` + `start_backmapping.py` scripts. The new workflow replaces all of this with a single `settings.yaml` → `backmap-prep` → LAMMPS pipeline.

Current `backmap-prep` capabilities (Phase 1+2): single molecule type, simple atom lists per bead, cross bonds/angles/dihedrals, tabulated potentials, harmonic potentials, Ryckaert-Bellemans dihedrals. Network features (degree-dependent beads, active sites, charge management) are Phase 3 and out of scope for this change.

## Goals / Non-Goals

**Goals:**
- Migrate PE (2:1 UA mapping), PE4 (4:1 UA mapping), PE-10 (10:1 AA mapping), and PE-AA (2:1 AA mapping with explicit H) examples to the `backmap-prep` YAML workflow
- Migrate the melamine example (simple small molecule, no reactive network features needed for the basic case)
- Each example is self-contained with `settings.yaml`, source GROMACS files, generated LAMMPS files, and documentation
- Validate `backmap-prep` against diverse mapping ratios and force field types
- Provide regression test targets for CI

**Non-Goals:**
- Migrating network_backmapping examples (epoxy, pete, mf, ab2, abx) — these require Phase 3 features (degree-dependent beads, cross-links between molecule types, charge management)
- Implementing new `backmap-prep` features — we use existing Phase 1/2 capabilities only
- Performance benchmarking or production-quality simulation parameters
- Porting the espressopp simulation engine or `prepare_files.py` logic

## Decisions

### Decision 1: Migrate linear examples first, defer networks

**Choice:** Migrate pe, pe4, pe_10, pe_aa, and melamine. Defer all network_backmapping examples.

**Rationale:** The linear polymer examples exercise the core `backmap-prep` pipeline without requiring Phase 3 features. They cover the most important dimension of variation: mapping ratio (2:1 through 10:1) and force field type (UA vs AA). Network examples need `atoms_by_degree`, multi-molecule systems, and cross-link handling that isn't implemented yet.

**Alternatives considered:**
- Migrate everything now: Would require implementing Phase 3 features, scope creep
- Only migrate PE: Would miss the chance to validate different mapping ratios and force fields

### Decision 2: Directory structure mirrors dodecane pattern

**Choice:** Each example gets `examples/<name>/` with: `settings.yaml`, source GRO/TOP/XVG files, generated `<name>.data` and `in.<name>`, `README.md`, and optional verification scripts.

**Rationale:** Consistency with the existing dodecane example. Users can follow the same workflow pattern regardless of which system they're working with.

### Decision 3: Copy source files, don't symlink to bakery

**Choice:** Copy all necessary GROMACS source files (GRO, TOP, XVG tables) from bakery into each example directory. Rename to consistent conventions (e.g., `topol_aa.top`, `topol_cg.top`, `cg_conf.gro`).

**Rationale:** Examples must be self-contained within the lammps-user-backmapping repo. Users shouldn't need the bakery repo. Consistent naming helps documentation and discoverability.

**Alternatives considered:**
- Git submodule for bakery: Adds complexity, bakery has ESPResSo++ deps
- Reference bakery files by path: Not portable, breaks for other users

### Decision 4: Translate bakery XML → YAML settings manually

**Choice:** Manually write each `settings.yaml` using the bakery XML as a reference, not automated conversion.

**Rationale:** The XML and YAML schemas are structurally different (VOTCA-style XML vs flat YAML sections). Manual translation ensures correctness and lets us optimize the YAML for readability. The number of examples (5) is small enough that manual work is practical.

### Decision 5: Melamine as non-network molecule

**Choice:** Migrate melamine as a simple 3-bead molecule without reactive network features.

**Rationale:** The basic melamine molecule (3 CG beads, 9 AT atoms) can be backmapped without degree-dependent beads or cross-links. The network aspects (MF polymerization) are a separate concern handled in Phase 3. This gives users a non-linear molecule example.

### Decision 6: Use reduced system sizes for quick testing

**Choice:** Provide examples with small system sizes (10-25 molecules) suitable for quick validation. Document original bakery system sizes for reference.

**Rationale:** The dodecane example uses 10 chains. Keeping PE examples at ~10-25 chains allows fast testing and movie generation. Users can scale up by modifying the CG system files.

## Risks / Trade-offs

- **[Risk: Schema gaps]** The PE examples use CG dihedrals which require Phase 2 `cross_interactions.dihedrals` support. If this isn't fully working, PE examples won't run.
  → Mitigation: Test dihedral support during migration; file bugs if needed. Simpler PE variants without dihedrals can ship first.

- **[Risk: Atom name mismatches]** Bakery uses specific atom naming in GRO files that must exactly match `settings.yaml` atom references (`chain:mol:atom` format).
  → Mitigation: Validate generated LAMMPS data against bakery-generated hybrid files for correctness.

- **[Risk: Table format differences]** GROMACS XVG tables from bakery may have different column layouts or units than what `backmap-prep` expects.
  → Mitigation: Verify table conversion output against known-good LAMMPS tables. The dodecane tables already work, and PE tables use the same format.

- **[Risk: Large files in git]** PE systems with 75 chains have larger GRO files (~100KB+).
  → Mitigation: Use reduced system sizes (10-25 chains). Git handles 100KB text files fine.

- **[Trade-off: Incomplete example set]** Shipping only linear polymers + melamine may disappoint users who need network examples.
  → Mitigation: Document this clearly. Network examples are Phase 3 scope with their own migration task.

## Open Questions

1. Should PE-10 (all-atom, ~36 atoms per bead) use the full OPLS/AA force field with explicit charges, or a simplified UA variant for the initial example?
2. What reduced system size should we target for PE examples? 10 chains (like dodecane) or 25 chains for better statistics?
3. Should we include visualization scripts (`visualize_backmap_movie.py`) for each example, or share a single generic visualization tool?
