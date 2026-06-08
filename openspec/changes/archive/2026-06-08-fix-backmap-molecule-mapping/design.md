## Context

The `fix backmap` manages the CG→AT mapping, COM tracking, and force distribution during backmapping simulations. It currently assumes one CG atom per molecule (keyed by `molecule ID → single CG local index`). This breaks for multi-bead molecules like dodecane (6 CG beads per molecule), causing wrong mass sums, incorrect forces, and simulation crashes.

The pair, bond, and angle backmap styles are correct — they use per-atom lambda values from the fix and don't depend on the molecule map. Only the fix itself needs redesign.

## Goals / Non-Goals

**Goals:**

- Fix the CG→AT mapping to support multiple CG beads per molecule, each mapping to its own AT atom subset.
- Support multiple CG atom types (e.g., both A and B in dodecane).
- Maintain backward compatibility: single-CG-type inputs still work unchanged.
- Fix the ghost atom mass accounting bug.
- Fix the dodecane example to produce a stable backmapping simulation.

**Non-Goals:**

- Non-uniform AT-atoms-per-bead within a molecule (all beads in a molecule have the same AT atom count) — deferred until a real use case requires it.
- Changing pair/bond/angle styles — these are correct.
- Changing the YAML settings schema (the Python writer emits the fix command; the YAML already has the mapping info).
- Two-phase backmapping support (separate change).

## Decisions

### D1: Multiple CG types via extended `cg_type` keyword

**Decision:** Extend `cg_type` to accept one or more type IDs (read until the next recognized keyword). Single-type usage is unchanged.

```
# Backward compatible (single CG type):
fix bm all backmap cg_type 1 alpha 0.0001 ...

# New (multiple CG types):
fix bm all backmap cg_type 1 2 alpha 0.0001 ...
```

**Why over alternatives:**
- *Count-prefixed list* (`cg_types 2 1 2`) — non-standard for LAMMPS and confusing.
- *Renamed keyword* (`cg_types 1 2`) — breaking change with no benefit.
- *Multiple keyword instances* (`cg_type 1 cg_type 2`) — verbose, no LAMMPS precedent.

**Implementation:** In the constructor argument parser, after seeing `cg_type`, read integers with `strtol` until the next token is a recognized keyword (`alpha`, `lambda0`, `nonuniform`, `phase`) or isn't a valid integer. Store in `std::set<int> cg_types`.

### D2: Per-CG-bead mapping via atom ordering convention

**Decision:** Infer the CG bead → AT atom mapping from atom tag ordering within each molecule. The Python generator already writes atoms in this order:

```
[CG1, CG2, ..., CGn, AT1_1, AT1_2, AT2_1, AT2_2, ..., ATn_k]
```

Within each molecule, the fix:
1. Collects all CG atoms and sorts by global tag → gives CG bead order.
2. Collects all AT atoms and sorts by global tag → gives AT atom order.
3. Computes `atoms_per_bead = num_AT / num_CG` and validates `num_AT % num_CG == 0`.
4. Assigns AT atoms `[i*apb, (i+1)*apb)` to CG bead `i`.

**Why over alternatives:**
- *Explicit mapping file* — over-engineering for MVP; the convention is already enforced by the Python generator.
- *Bond connectivity inference* — requires traversing the bond topology during `init()`, which is complex and fragile. Better suited for a future enhancement.
- *Custom data file section* — non-standard LAMMPS data format; confuses users.

**Trade-off:** This breaks if users manually create data files with non-sequential ordering. Acceptable because the Python generator is the only supported data file creator.

### D3: Replace `MolMap` with `BeadMap` data structure

**Decision:** Replace the current `std::map<tagint, MolMap>` (molecule → single CG) with `std::vector<BeadMap>` where each entry maps one CG bead to its AT atoms:

```cpp
struct BeadMap {
  int cg_local;                // local index of the CG atom
  std::vector<int> at_local;   // local indices of AT atoms for this bead
  double at_mass_sum;          // sum of AT atom masses (for validation)
};
std::vector<BeadMap> bead_map;
```

The vector is rebuilt on every `init()` call (after domain decomposition). It only contains beads whose CG atom is **local** to this processor. AT atoms may be local or ghost.

### D4: Mass validation against CG bead mass (not molecule mass)

**Decision:** Validate `mass[cg_type] == sum(mass[at_type_i])` per CG bead, not per molecule. For dodecane:
- CG bead A (type 1, mass 29.062) vs AT atoms CH3+CH2 (15.035+14.027 = 29.062) ✓
- CG bead B (type 2, mass 28.054) vs AT atoms CH2+CH2 (14.027+14.027 = 28.054) ✓

**Ghost atom fix:** When computing `at_mass_sum`, use `atom->mass[atom->type[at]]` (per-type mass array), which is always available regardless of local/ghost status. This avoids the double-counting bug in the current code.

### D5: Force distribution and COM tracking per bead

**Decision:** Both `post_force()` and `end_of_step()` iterate over `bead_map` instead of `mol_map`. Each CG bead only interacts with its own AT atoms:

```
post_force: f[AT_i] += (m_AT_i / m_CG_bead) * f[CG_bead];  f[CG_bead] = 0
end_of_step: r[CG_bead] = COM(AT atoms of this bead)
```

The logic is identical to the current code, just scoped per-bead instead of per-molecule.

### D6: Python writer changes

**Decision:** Update `writers.py` to emit `cg_type T1 T2 ...` with all CG type IDs from `system.atom_types` where `is_cg=True`. Replace the single `system.cg_type_id` reference with a list. Keep the `cg_type` keyword (not `cg_types`) for LAMMPS convention consistency.

### D7: Fix dodecane data file coordinates

**Decision:** Ensure the Python generator wraps all atom coordinates inside the simulation box `[0, L)` before writing. Add a wrapping step in `write_lammps_data()` that applies `x = x % box_length` for each dimension. This resolves the "Inconsistent image flags" warning.

## Risks / Trade-offs

- **[Risk] Atom ordering convention breaks manual data files** → Mitigation: Document the convention clearly. The Python generator is the supported path; manual data files are advanced usage.
- **[Risk] Non-uniform AT-per-bead (e.g., terminal beads with 3 atoms, middle beads with 2)** → Mitigation: The current design detects this case (`num_AT % num_CG != 0`) and aborts with a clear error. Future enhancement can add a mapping file.
- **[Risk] Breaking change to `cg_type` syntax** → Mitigation: The change is backward compatible — `cg_type T` still works for single-type molecules.
- **[Risk] `build_molecule_map` performance with many molecules** → Mitigation: O(N) with N = total atoms. Acceptable for the typical system sizes in backmapping (thousands, not millions of atoms).
