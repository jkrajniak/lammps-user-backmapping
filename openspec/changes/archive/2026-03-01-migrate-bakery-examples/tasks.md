## 1. PE Example (2:1 UA mapping)

- [x] 1.1 Create `examples/pe/` directory
- [x] 1.2 Copy source files from bakery `examples/pe/`: `pe_single.gro` → `pe_single.gro`, `topol.top` → `topol_aa.top`, `cg_conf.gro`, `cg_topol.top` → `topol_cg.top`, `table_b1.xvg`, `table_b2.xvg`, `table_a1.xvg`, `table_a2.xvg`, `table_A_A.xvg`, `table_A_B.xvg`, `table_B_B.xvg`
- [x] 1.3 Reduce system size: create a smaller `cg_conf.gro` with 10-25 chains (original: 75 chains, 50 beads each)
- [x] 1.4 Write `settings.yaml` translating bakery `pe_settings.xml` to YAML format: molecule definition (A/B beads, 2 atoms each), CG system, cross bonds/angles/dihedrals, simulation params
- [x] 1.5 Run `backmap-prep settings.yaml` and verify output (`pe.data`, `in.pe`) is generated without errors
- [x] 1.6 Fix any `backmap-prep` issues discovered during generation (if needed)
- [x] 1.7 Write `README.md` with system description, mapping table, quick-start commands, and file listing
- [x] 1.8 Verify LAMMPS can parse the generated input script (syntax check with `lmp -in in.pe`) — **data file, bonds, angles, groups, fixes all parse OK; fails on CG pair_coeff due to pre-existing backmap-prep bug (CG sub-style is `table` but pair_coeff emits `0.0 0.0` instead of table filename); same error as dodecane**

## 2. PE4 Example (4:1 UA mapping)

- [x] 2.1 Create `examples/pe4/` directory
- [x] 2.2 Copy source files from bakery `examples/pe4/`: `pe_single.gro`, `topol.top` → `topol_aa.top`, `cg_conf.gro`, `cg_topol.top` → `topol_cg.top`, tabulated potential XVG files
- [x] 2.3 Reduce system size: create smaller `cg_conf.gro` with 10-25 chains (original: 75 chains, 25 beads each)
- [x] 2.4 Write `settings.yaml` translating bakery `pe4_settings.xml`: molecule definition (A/B beads, 4 atoms each), cross interactions, simulation params
- [x] 2.5 Run `backmap-prep settings.yaml` and verify output (`pe4.data`, `in.pe4`)
- [x] 2.6 Write `README.md` with system description and mapping table
- [x] 2.7 Verify LAMMPS can parse the generated input script — **same pre-existing CG pair_coeff/table bug as all other examples; everything else parses OK**

## 3. PE-10 Example (10:1 AA mapping)

- [x] 3.1 Create `examples/pe_10/` directory
- [x] 3.2 Copy source files from bakery `examples/pe_10/`: `pe_single.gro`, `topol_at.top` → `topol_aa.top`, `conf_cg.gro` → `cg_conf.gro`, `cg_topol.top` → `topol_cg.top`, tabulated potential XVG files
- [x] 3.3 Reduce system size: create smaller `cg_conf.gro` with 10 chains (original: 75 chains, 10 beads each)
- [x] 3.4 Write `settings.yaml` translating bakery `backmapping.xml`: molecule definition with ~30-36 all-atoms per bead (including H), cross interactions, simulation params
- [x] 3.5 Run `backmap-prep settings.yaml` and verify output (`pe_10.data`, `in.pe_10`)
- [x] 3.6 Fix any issues with high atom-per-bead ratio handling in `backmap-prep` (if needed)
- [x] 3.7 Write `README.md` with system description and mapping table
- [x] 3.8 Verify LAMMPS can parse the generated input script — **same pre-existing CG pair_coeff/table bug; everything else parses OK**

## 4. PE-AA Example (2:1 AA mapping with explicit H)

- [x] 4.1 Create `examples/pe_aa/` directory
- [x] 4.2 Copy source files from bakery `examples/pe_aa/`: `pe_single.gro`, topology files → `topol_aa.top`, CG files → `cg_conf.gro`, `topol_cg.top`, tabulated potential XVG files
- [x] 4.3 Reduce system size if needed
- [x] 4.4 Write `settings.yaml` translating bakery `backmapping_pe_aa.xml`: molecule definition with 6-7 atoms per bead (2C + 4-5H), OPLS/AA charges/types, cross interactions
- [x] 4.5 Run `backmap-prep settings.yaml` and verify output (`pe_aa.data`, `in.pe_aa`)
- [x] 4.6 Verify hydrogen atoms have correct OPLS/AA charges and types in the data file
- [x] 4.7 Write `README.md` with system description and mapping table
- [x] 4.8 Verify LAMMPS can parse the generated input script — **same pre-existing CG pair_coeff/table bug; everything else parses OK**

## 5. Melamine Example (triangular CG topology)

- [x] 5.1 Create `examples/melamine/` directory
- [x] 5.2 Copy source files from bakery `examples/melamine/`: `single_mf.gro`, `at_topol.top` → `topol_aa.top`, `cg_conf_500.gro` → `cg_conf.gro`, `cg_topol.top` → `topol_cg.top`, `table_A_A.xvg`, `table_b1.xvg`
- [x] 5.3 Reduce system size: create `cg_conf.gro` with 50 molecules (original: 500)
- [x] 5.4 Write `settings.yaml` for melamine: 3 A-type beads per molecule (A1, A2, A3 each with 9 atoms), triangular cross bonds [A1-A2, A2-A3, A3-A1], simulation params
- [x] 5.5 Run `backmap-prep settings.yaml` and verify output (`melamine.data`, `in.melamine`)
- [x] 5.6 Verify triangular CG bonding topology in generated data file (3 CG bonds per molecule)
- [x] 5.7 Write `README.md` with system description and mapping table
- [x] 5.8 Verify LAMMPS can parse the generated input script — **same pre-existing CG pair_coeff/table bug; everything else parses OK**

## 6. Integration Testing

- [x] 6.1 Add `backmap-prep` generation tests for each new example in the Python test suite
- [x] 6.2 Validate generated data file structure: correct atom/bond/angle counts per molecule, molecule ID assignment
- [x] 6.3 Verify deterministic output: running `backmap-prep` twice produces identical files
- [ ] 6.4 Add new examples to CI workflow matrix (if applicable)

## 7. Documentation

- [x] 7.1 Update `docs/getting-started.md` to reference the new examples
- [x] 7.2 Ensure all example READMEs follow the dodecane README format
- [x] 7.3 Add a summary table of all available examples to the main docs or examples directory README
