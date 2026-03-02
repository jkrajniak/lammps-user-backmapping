## ADDED Requirements

### Requirement: Water2 integration test

The package SHALL include a water2 integration test that validates the full backmapping workflow on a single-bead water system (SPC water mapped to WCG coarse-grained beads). This is the simplest backmapping system — no cross-CG bonds, only non-bonded lambda-weighted interactions.

The test SHALL:
1. Provide a YAML settings file for 3456 SPC water molecules mapped to WCG beads
2. Provide GROMACS source files: CG coordinates (`.gro`), AT topology (`.top`), CG pair table (`.xvg`)
3. Run `backmap-prep` to generate LAMMPS data file, input script, and tables
4. Run LAMMPS backmapping simulation (lambda 0→1)
5. Validate results against acceptance criteria

Acceptance criteria:
- Lambda ramp completes from 0.0 to 1.0 without simulation failure
- Final AT structure has correct water geometry: O-H bond length within 5% of 1.0 Å, H-O-H angle within 5% of 109.47°
- Density of the final AT system is within 2% of pure SPC/E water density (~1.0 g/cm³ at 300K)
- CG positions track COM of AT atoms throughout the simulation (max deviation < 0.01 Å)

#### Scenario: Water2 backmapping runs to completion
- **WHEN** the water2 test is executed with 3456 SPC water molecules and alpha=0.001
- **THEN** the simulation SHALL complete without errors and lambda SHALL reach 1.0 for all molecules

#### Scenario: Water2 structural validation
- **WHEN** the water2 backmapping completes
- **THEN** the final AT structure SHALL have mean O-H bond length within [0.95, 1.05] Å and mean H-O-H angle within [104.0, 115.0] degrees

#### Scenario: Water2 density validation
- **WHEN** the water2 backmapping completes
- **THEN** the system density SHALL be within 2% of the reference pure SPC/E density

#### Scenario: Water2 COM tracking
- **WHEN** the simulation is running
- **THEN** at every output frame, each CG atom position SHALL equal the COM of its AT atoms within 0.01 Å tolerance

### Requirement: PE integration test

The package SHALL include a PE (polyethylene) integration test that validates backmapping of a linear polymer with cross-CG bonds, cross-CG angles, and tabulated CG potentials. This is the first test requiring `bond_style backmap/*` and `angle_style backmap/*`.

The test SHALL:
1. Provide a YAML settings file for a PE system (50 CG beads → 100 AT united-atom carbons)
2. Provide GROMACS source files with AT topology including intra-CG and cross-CG bonds/angles
3. Provide tabulated CG bond and pair potentials
4. Run `backmap-prep` to generate LAMMPS files
5. Run LAMMPS backmapping simulation (lambda 0→1)
6. Validate results

Acceptance criteria:
- Full backmapping from λ=0 to λ=1 completes without simulation failure
- AT bond lengths are within 5% of equilibrium values
- AT angle values are within 10% of equilibrium values
- Cross-CG AT bond forces scale correctly with λ² (verified by thermo output)
- Cross-CG CG bond forces scale correctly with 1−λ² (verified by thermo output)
- No unphysical atom overlaps (minimum pairwise distance > 0.5 Å)

#### Scenario: PE backmapping runs to completion
- **WHEN** the PE test is executed with 50 CG beads mapped to 100 AT carbons
- **THEN** the simulation SHALL complete without errors and lambda SHALL reach 1.0

#### Scenario: PE bond length validation
- **WHEN** the PE backmapping completes
- **THEN** the mean C-C bond length SHALL be within 5% of the equilibrium value (1.53 Å for united-atom PE)

#### Scenario: PE angle validation
- **WHEN** the PE backmapping completes
- **THEN** the mean C-C-C angle SHALL be within 10% of the equilibrium value (~111°)

#### Scenario: PE tabulated CG potential
- **WHEN** the PE system uses tabulated CG bond potentials (`.table` format)
- **THEN** the `backmap/table cg` bond forces SHALL be read from the table and weighted by `1−λ²`

### Requirement: MPI parallel correctness test

The package SHALL include an MPI parallel correctness test that validates the backmapping produces identical (or numerically equivalent) results when run on multiple processors.

The test SHALL:
1. Run the water2 test system on 1 processor (serial reference)
2. Run the same system on 4 processors
3. Run the PE test system on 1 processor
4. Run the PE system on 4 processors
5. Compare trajectories and thermodynamic output between serial and parallel runs

Acceptance criteria:
- Thermodynamic quantities (total energy, temperature, pressure) SHALL match between serial and parallel runs within floating-point tolerance (relative error < 1e-10 or machine epsilon level)
- Final atom positions SHALL match within floating-point tolerance
- Molecule mapping SHALL be correctly rebuilt after atoms migrate between processor domains
- Ghost atom lambda values SHALL be correctly communicated across processor boundaries

#### Scenario: Water2 MPI consistency
- **WHEN** the water2 system is run on 1 processor and on 4 processors with the same random seed
- **THEN** final atom positions SHALL match within floating-point tolerance (< 1e-8 Å)

#### Scenario: PE MPI consistency
- **WHEN** the PE system is run on 1 processor and on 4 processors
- **THEN** final atom positions SHALL match within floating-point tolerance

#### Scenario: Molecule migration during simulation
- **WHEN** molecules migrate between processor domains during the PE simulation on 4 processors
- **THEN** the molecule map SHALL be correctly rebuilt and CG-AT associations SHALL remain valid

#### Scenario: Ghost atom lambda communication
- **WHEN** atoms of a molecule are split across processor boundaries
- **THEN** ghost atom lambda values SHALL be correctly communicated via `comm->forward_comm()` and the weight computation SHALL use the correct lambda values

### Requirement: Epoxy integration test

The package SHALL include an epoxy (RIM135) integration test that validates backmapping of a reactive multi-molecule system. This test validates the Phase 3 reactive network features of the Python generator.

The test SHALL:
1. Provide a YAML settings file for a cured epoxy network with EPO (EPON 828), IPD (isophorone diamine), and HDD (1,6-hexanediol diglycidyl ether) molecule types
2. Include degree-dependent bead definitions, active sites, and charge management rules
3. Run `backmap-prep` to generate LAMMPS files
4. Run LAMMPS backmapping simulation
5. Compare results against bakery/ESPResSo++ reference data

Acceptance criteria:
- Degree-dependent beads correctly identified and mapped based on CG topology connectivity
- Active sites correctly placed in the AT structure
- Charge transfers applied correctly upon bond formation
- Final AT structure is chemically valid (correct valence, no broken bonds, net charge conserved)
- Structural properties (RDF, density) within 5% of ESPResSo++ bakery reference

#### Scenario: Epoxy degree-dependent mapping
- **WHEN** a cured epoxy network has EPO beads with degrees 0, 1, and 2
- **THEN** each bead SHALL be mapped to the correct AT fragment for its degree, with the correct number of atoms

#### Scenario: Epoxy charge conservation
- **WHEN** charge transfers are applied during the backmapping setup
- **THEN** the total system charge SHALL be conserved (sum of all atomic charges unchanged within 1e-6 e)

#### Scenario: Epoxy comparison with ESPResSo++
- **WHEN** the epoxy backmapping completes
- **THEN** the radial distribution function of AT atoms SHALL match the ESPResSo++ reference within 5% at peak positions

### Requirement: Test infrastructure

Each integration test SHALL be a self-contained directory under `examples/` with a consistent structure:

```
examples/<test_name>/
  settings.yaml          # backmap-prep input
  source/                # GROMACS source files (.gro, .top, .xvg)
  reference/             # reference data for validation (optional)
  run_test.sh            # test runner script
  analyze.py             # Python analysis and validation script
  README.md              # test description and expected results
```

The `run_test.sh` script SHALL:
1. Check prerequisites (LAMMPS binary, `backmap-prep` command)
2. Run `backmap-prep settings.yaml`
3. Run LAMMPS with the generated input
4. Run `analyze.py` to validate results
5. Print PASS/FAIL status with summary metrics

The `analyze.py` script SHALL:
- Accept the LAMMPS output directory as argument
- Compute validation metrics (bond lengths, angles, density, RDF, etc.)
- Compare against acceptance criteria
- Print structured results (metric, value, threshold, pass/fail)
- Exit with code 0 on PASS, non-zero on FAIL

#### Scenario: Test runner prerequisites check
- **WHEN** `run_test.sh` is executed without LAMMPS installed
- **THEN** it SHALL print an error message indicating the missing prerequisite and exit with code 1

#### Scenario: Test runner full execution
- **WHEN** `run_test.sh` is executed with all prerequisites available
- **THEN** it SHALL run the full pipeline (generate → simulate → analyze) and report PASS/FAIL

#### Scenario: Analysis script structured output
- **WHEN** `analyze.py` completes
- **THEN** it SHALL print a summary table with columns: metric name, measured value, acceptance threshold, pass/fail status

### Requirement: Integration test coverage for migrated examples

The integration test suite SHALL validate the migrated bakery examples in addition to the existing dodecane test. Each migrated example SHALL have a test that:

1. Runs `backmap-prep settings.yaml` to generate LAMMPS data and input files
2. Validates the generated data file structure (atom counts, bond counts, molecule IDs)
3. Optionally runs a short LAMMPS simulation to verify the input script is functional

#### Scenario: PE example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe/settings.yaml`
- **THEN** the generated `pe.data` SHALL have the correct number of atoms, bonds, angles, and molecule IDs for the PE system

#### Scenario: PE4 example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe4/settings.yaml`
- **THEN** the generated `pe4.data` SHALL have the correct number of atoms per molecule (25 CG + 100 AT = 125) and correct bond topology

#### Scenario: PE-10 example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe_10/settings.yaml`
- **THEN** the generated `pe_10.data` SHALL have the correct atom counts reflecting the high AT-to-CG ratio

#### Scenario: PE-AA example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/pe_aa/settings.yaml`
- **THEN** the generated `pe_aa.data` SHALL include hydrogen atoms with correct charges and OPLS/AA types

#### Scenario: Melamine example generates valid output
- **WHEN** the CI pipeline runs `backmap-prep` on `examples/melamine/settings.yaml`
- **THEN** the generated `melamine.data` SHALL have triangular CG bonding (3 CG bonds per molecule) and correct atom counts (30 atoms per molecule)

#### Scenario: All examples produce consistent output
- **WHEN** `backmap-prep` is run on any migrated example twice with the same settings
- **THEN** the output files SHALL be byte-identical (deterministic generation)
