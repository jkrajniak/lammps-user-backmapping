# Design — CG dihedrals

## Decisions

1. **Three C++ styles:** static `ryckaert` (intra-bead), `backmap/ryckaert` (cross AT),
   `backmap/table` (CG tabulated). Lambda weight on endpoints **i** and **l**.

2. **RB conversion:** GROMACS func-3 `C_lammps[n] = (-1)^n × energy(C_gromacs[n])`.

3. **Param resolution:** Hybrid TOP lines often have func only; lookup
   `(type_i, type_j, type_k, type_l)` in parsed `[ dihedraltypes ]` from forcefield includes.
   `prep.forcefield_dir` or `GMXDATA`/`GROMACS_DATA` for include paths.

4. **Routing:** All four atoms in same bead → static `ryckaert`; cross-bead → `backmap/ryckaert at`;
   func-8 → `backmap/table cg` via `table_d{tablenr}.xvg`.

5. **Validation gate:** rim135.data has ~33,421 dihedrals; step-0 LAMMPS parse with rebuilt binary.
