## Context

The MVP implementation of the `USER-BACKMAP` LAMMPS package is complete. It provides:

- **`fix backmap`** — single-phase lambda ramp (0→1), CG-AT mapping via molecule IDs, COM tracking, CG force distribution
- **`pair_style backmap`** — lambda-weighted non-bonded pair forces delegating to AT/CG sub-styles
- **`bond_style backmap/harmonic`**, **`bond_style backmap/table`** — lambda-weighted cross-CG bond forces
- **`angle_style backmap/harmonic`** — lambda-weighted cross-CG angle forces
- **`backmap_lambda.h`** — shared helper for lambda access and weight computation across all styles
- **`backmap-prep`** Python CLI — generates LAMMPS data files, input scripts, and interaction tables from GROMACS topology + YAML settings

All styles use the same weighting scheme: `w_AT = λ_i × λ_j`, `w_CG = 1 − λ_i × λ_j`, reading per-atom lambda from `fix backmap` via `fix->extract()`.

This design covers the remaining work: dihedral styles (D4/D5), two-phase backmapping (B8), Python generator extensions (E9–E11), and integration testing (F1/F3/F4/F5).

### Constraints

- Must follow the same architectural patterns established in the MVP (lambda access via `backmap_lambda.h`, `at`/`cg` keyword convention, LAMMPS coding style)
- Dihedral styles must integrate with `dihedral_style hybrid` alongside standard LAMMPS dihedral styles
- Two-phase protocol must be controllable from the LAMMPS input script without recompilation
- Integration tests need compiled LAMMPS with the package + Python environment for `backmap-prep`
- Reactive network support (E10) targets the ESPResSo++ bakery epoxy workflow

## Goals / Non-Goals

**Goals:**

- Implement `dihedral_style backmap/ryckaert` and `dihedral_style backmap/table` for lambda-weighted dihedral forces
- Add two-phase backmapping support to `fix backmap` via `fix_modify` phase switching
- Extend `backmap-prep` to generate dihedral sections in data files and input scripts
- Validate the full package against ESPResSo++ reference results on water2, PE, and MPI parallel runs
- Add reactive network support to the Python generator for epoxy-like systems

**Non-Goals:**

- Reimplementing any MVP features (bonds, angles, pair style, fix core mechanics)
- Changing the YAML settings schema structure (it is stable; we add fields, not restructure)
- Spatial AdResS regions (out of scope — that is `lammps-user-adress`)
- Non-linear lambda schedules (sigmoid, staged ramps) — defer to a future change
- GPU/Kokkos acceleration

## Decisions

### D1: Dihedral styles mirror existing bond/angle architecture

**Decision**: Implement `dihedral_backmap_ryckaert` and `dihedral_backmap_table` following the exact same patterns as the existing `bond_backmap_harmonic`, `bond_backmap_table`, and `angle_backmap_harmonic` styles. Specifically:

- Use `backmap_lambda.h` (`find_fix_backmap()`, `get_lambda()`, `compute_backmap_weight()`) for lambda access and weight computation
- Accept `at`/`cg` keyword as the first `dihedral_coeff` argument to set weighting direction
- Weight is computed from lambda values of the first and last atoms (i and l) of the dihedral quadruplet i-j-k-l, matching ESPResSo++'s convention
- Scale force, energy, and virial by the same weight factor
- Skip computation when weight < 1e-10 (the `is_almost_zero()` optimization)

**Rationale**: The MVP established a clean, consistent pattern for lambda-weighted styles. The dihedral styles are specified in the bonded-backmap-styles spec from the MVP change. Following the same pattern reduces cognitive overhead and ensures all styles behave consistently. The weight-from-endpoints convention (atoms i and l) matches ESPResSo++'s `FixedQuadrupleListAdressInteractionTemplate`.

**Input syntax**:
```
dihedral_style hybrid backmap/ryckaert backmap/table
dihedral_coeff 1 backmap/ryckaert at C0 C1 C2 C3 C4 C5
dihedral_coeff 2 backmap/table cg table_d1.table ENTRY
```

### D2: Ryckaert-Bellemans dihedral uses LAMMPS sign convention

**Decision**: The `backmap/ryckaert` style uses the standard LAMMPS Ryckaert-Bellemans formula:

```
E = w × Σ(n=0..5) Cn × cos^n(φ)
```

where φ is the dihedral angle and Cn are the six RB coefficients. The LAMMPS convention defines φ differently from GROMACS (LAMMPS uses the "polymer convention" where trans = 180°, GROMACS uses trans = 0°). The Python generator handles the coefficient transformation during table generation.

**Rationale**: LAMMPS's built-in `dihedral_style multi/harmonic` and older RB implementations follow this convention. Matching it avoids user confusion. The GROMACS-to-LAMMPS coefficient mapping is: `C_lammps[n] = (-1)^n × C_gromacs[n]` for the cosine power series. This transformation is a straightforward sign flip on odd-indexed coefficients and belongs in the Python generator, not the C++ style.

### D3: Tabulated dihedral follows LAMMPS `dihedral_style table` format

**Decision**: The `backmap/table` dihedral style uses the same table file format and interpolation approach as LAMMPS's built-in `dihedral_style table`:

- Table entries: angle (degrees), energy (kcal/mol), force (kcal/mol/radian)
- Linear or spline interpolation, selectable via style argument
- Periodic boundary handling for the dihedral angle (-180 to 180)

**Rationale**: Users familiar with LAMMPS's native `dihedral_style table` will find the format identical. The only addition is the lambda weighting layer.

### D4: Two-phase backmapping via `fix_modify` state machine

**Decision**: Extend `fix backmap` with a simple two-state machine controlled by `fix_modify`:

| Phase | Lambda behavior | CG forces | AT forces | Entered via |
|-------|----------------|-----------|-----------|-------------|
| **Phase 1** (default) | λ ramps 0→1 at `rate` | Full strength (static, no lambda weighting) | Weighted by λ² | `fix backmap ... phase 1` or default |
| **Phase 2** | λ resets to 0, ramps 0→1 at `rate` | Weighted by (1−λ²) | Weighted by λ² | `fix_modify bm phase 2` |

In Phase 1, cross-CG CG bonds are at full strength (ignoring lambda) while cross-CG AT bonds ramp in. When Phase 2 activates, lambda resets to 0 and CG forces begin ramping out as AT forces ramp in — the standard AdResS-style interpolation.

**`fix_modify` interface**:
```
fix_modify bm phase 2            # switch to phase 2, reset lambda to 0
fix_modify bm phase 1            # switch to phase 1 (rarely needed)
fix_modify bm active yes/no      # existing: enable/disable lambda ramp
```

**Implementation**: Add a `phase` integer (1 or 2) to the fix. During Phase 1, the fix communicates `phase=1` through `extract("phase")` so that `backmap/*` CG styles can apply full strength instead of `(1−λ²)`. In Phase 2 (the default behavior from MVP), styles behave as currently implemented.

**Rationale**: The two-phase protocol is described in Krajniak et al. 2016, Section 2.3. Phase 1 first establishes atomistic connectivity under full CG restraint. Phase 2 then gradually removes CG forces. This staged approach improves stability for complex molecules where turning on AT interactions and removing CG interactions simultaneously causes instability. Using `fix_modify` for phase switching is idiomatic LAMMPS — it allows scriptable transitions at specific timesteps without modifying the fix definition.

**Alternatives considered**:
- *Separate fixes for each phase*: Would require duplicating state management and transferring lambda arrays between fixes. Overly complex.
- *Phase encoded in lambda range (e.g., λ < 0 = phase 1)*: Confusing semantics, negative lambda is meaningless physically.
- *Automatic phase switching at λ=1*: Less flexible — users may want to equilibrate at λ=1 before switching phases, or use only phase 1 for simple systems.

### D5: Python generator dihedral support extends existing architecture

**Decision**: Extend `backmap-prep` for dihedrals by adding to each existing module:

| Module | Addition |
|--------|----------|
| `schema.py` | `CrossDihedral` model already exists; add `cross_dihedrals` list to `CrossInteractions` if not present |
| `builder.py` | Add `LammpsDihedral` dataclass, populate dihedral topology from cross-interaction definitions |
| `writers.py` | Add `Dihedrals` section to data file writer; add `dihedral_style hybrid backmap/*` to input script |
| `parsers/top_parser.py` | Already parses `[ dihedrals ]` section — verify it handles function type 3 (RB) and type 8 (tabulated) |

**Rationale**: The existing architecture handles bonds and angles with a clean pattern: schema defines the interaction → builder creates topology entries → writer outputs LAMMPS format. Dihedrals follow the same pattern with minimal new code.

### D6: Integration tests use shell scripts with numerical validation

**Decision**: Integration tests are structured as self-contained directories under `examples/`, each containing:

1. A YAML settings file for `backmap-prep`
2. Reference GROMACS source files (`.gro`, `.top`, `.xvg` tables)
3. A `run_test.sh` script that: runs `backmap-prep`, runs LAMMPS, extracts metrics
4. A Python analysis script that compares results against acceptance criteria

**Test matrix**:

| Test | System | Validates | Acceptance criteria |
|------|--------|-----------|-------------------|
| F1: water2 | 3456 SPC water + WCG | Lambda ramp, COM tracking, force distribution | Final AT structure matches ESPResSo++ (RMSD < 0.1 Å), density within 2% of pure AT |
| F3: PE | 50 CG → 100 AT | Tabulated CG potentials, harmonic AT cross bonds | Full 0→1 backmapping completes, bond lengths within 5% of equilibrium |
| F4: MPI | water2 + PE on 4 procs | Ghost atom communication, domain decomposition | Identical trajectory to serial (bitwise or within floating-point tolerance) |

**Rationale**: Shell-script-based tests are the convention in LAMMPS's own examples. Python analysis scripts provide quantitative validation beyond "it didn't crash." Each test is independently runnable for debugging.

**Alternatives considered**:
- *pytest-based integration tests*: Would require LAMMPS Python bindings or subprocess calls. Shell scripts are simpler and more portable for C++ package testing.
- *LAMMPS regression test framework*: Exists but is internal to LAMMPS CI. Our tests need to work standalone.

### D7: Reactive network support in Python generator (Phase 3)

**Decision**: Extend `backmap-prep` to handle epoxy-like reactive systems by adding to the YAML schema:

- **Degree-dependent bead definitions**: A CG bead type can map to different AT fragments depending on its bonding degree (number of reacted cross-links)
- **Active sites**: Specific AT atoms within a bead that can form new bonds during the backmapping process
- **Charge management**: Per-degree charge assignments and charge equilibration rules

This mirrors the bakery/AdResSLab reactive network functionality. The implementation is contained entirely in the Python generator — no C++ changes needed. The generated LAMMPS files use the same `fix backmap` + `backmap/*` styles but with additional atom types and bond types for each degree variant.

**Rationale**: The reactive network feature is the primary differentiator for epoxy (RIM135) systems. The bakery XML settings file already encodes degree-dependent mappings; the YAML schema mirrors this structure. Keeping the complexity in the generator means the C++ package stays simple and general.

### D8: Non-GROMACS format support via parser abstraction (Phase 4)

**Decision**: Add parser implementations for PDB, LAMMPS data, and XYZ formats behind a common interface. The builder calls `parse_coordinates()` and `parse_topology()` without knowing the source format. Format detection uses file extension or an explicit `format` field in the YAML settings.

**Rationale**: The current parser module (`parsers/top_parser.py`, `parsers/gro_parser.py`) already returns standardized data structures. Adding new format parsers is a matter of implementing the same interface for different file formats. This is lowest priority (Phase 4) because GROMACS format covers all current use cases.

## Risks / Trade-offs

**[Risk] Two-phase protocol changes fix-style interaction contract** → Phase 1 requires CG styles to behave differently (full strength instead of `(1−λ²)`). This means `backmap/*` styles need to check the phase flag via `fix->extract("phase")`, adding a runtime check per force computation. Mitigation: the check is a single integer comparison, negligible cost. The `backmap_lambda.h` helper encapsulates this logic so individual styles don't need to handle it directly.

**[Risk] Dihedral angle convention mismatch between GROMACS and LAMMPS** → GROMACS RB dihedrals use trans=0° convention; LAMMPS uses trans=180°. If the coefficient transformation is wrong, dihedral energies will be incorrect and systems may blow up. Mitigation: unit test the coefficient transformation in the Python generator against known GROMACS → LAMMPS conversions. Validate with the PE integration test which uses RB dihedrals.

**[Risk] MPI parallel correctness for large molecules** → Molecules spanning multiple processor domains require ghost atom communication for COM calculation and force distribution. The MVP handles this for small molecules (water, 4 atoms). Larger molecules (PE, 100+ atoms) may trigger edge cases in ghost communication range. Mitigation: ensure `comm_forward` and `comm_reverse` sizes are correct for the larger molecule case. Test explicitly with 4+ MPI ranks on the PE system.

**[Risk] Reactive network complexity in Python generator** → Degree-dependent beads, active sites, and charge management are the most complex features in the bakery workflow. They interact: changing a bead's degree changes its atom types, which changes its charges, which affects the force field. Mitigation: implement incrementally (degree-dependent beads first, then active sites, then charges). Each sub-feature gets its own unit tests. The epoxy integration test (F5) validates the full pipeline.

**[Trade-off] Integration tests require external LAMMPS binary** → Tests can't run in the Python CI alone; they need a compiled LAMMPS with the package installed. This makes CI setup more complex. Mitigation: document the build prerequisites clearly. Provide a `Makefile` target that builds LAMMPS with the package for testing. Consider containerized CI in the future.

**[Trade-off] Phase 4 non-GROMACS parsers add maintenance burden** → Each new file format parser needs to be maintained and tested. PDB and XYZ parsing in particular have many edge cases. Mitigation: use well-tested libraries (MDAnalysis, Biopython) where possible rather than custom parsers. Keep parsers minimal — extract only the fields `backmap-prep` needs.
