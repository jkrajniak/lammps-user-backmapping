## Context

The current backmapping workflow uses ESPResSo++ (espressopp-bakery) with three custom components:

- **`DynamicResolution`** — integrator extension that increments per-particle λ by a fixed rate each timestep
- **`VelocityVerletHybrid`** — modified integrator that (a) skips CG virtual particles during integration, (b) updates CG position to COM of AT atoms, (c) distributes CG forces to AT atoms proportional to mass ratio
- **`FixedVSList`** — data structure mapping CG particle → list of AT particle IDs

The Python layer (AdResSLab) parses GROMACS topology, creates particles, sets up interactions, and drives the simulation through ESPResSo++'s Python bindings.

The target is LAMMPS, which has a fundamentally different architecture: input-script driven, fixes for extending integration, pair styles for force computation, and groups for particle selection.

### Constraints

- Must reproduce the same physics as the ESPResSo++ implementation
- Must work with standard LAMMPS releases (no core modifications)
- Must handle MPI parallelization (molecules split across processors)
- The package should be self-contained (no dependency on `lammps-user-adress`)

## Goals / Non-Goals

**Goals:**

- A LAMMPS user package (`USER-BACKMAP`) containing:
  - `fix backmap` — manages λ ramp, CG-AT mapping, COM tracking, CG force distribution
  - `pair_style backmap` — lambda-weighted non-bonded pair forces
  - `bond_style backmap/*`, `angle_style backmap/*`, `dihedral_style backmap/*` — lambda-weighted cross-CG bonded forces
- A Python tool that converts GROMACS topology + coordinates into LAMMPS-ready input files
- Validation against ESPResSo++ results on the water2 test system

**Non-Goals:**

- Spatial AdResS (position-dependent λ with explicit/hybrid/CG regions)
- Thermodynamic force correction
- H-AdResS (Hamiltonian adaptive resolution)
- Moving or spherical resolution regions
- Many-body potentials (only pairwise)
- Reverse direction (AT→CG coarse-graining)

## Decisions

### D1: Custom lambda-weighted styles for all interaction types

**Decision**: Implement custom LAMMPS styles that apply lambda weighting at force-computation time for both non-bonded (pair) and bonded (bond, angle, dihedral) interactions. This requires:

- **`pair_style backmap`** — wraps AT and CG pair sub-styles, applies lambda weighting to non-bonded forces
- **`bond_style backmap/harmonic`**, **`bond_style backmap/table`** — lambda-weighted bond styles for cross-CG bonds
- **`angle_style backmap/harmonic`** — lambda-weighted angle style for cross-CG angles
- **`dihedral_style backmap/ryckaert`**, **`dihedral_style backmap/table`** — lambda-weighted dihedral styles for cross-CG dihedrals

All custom styles read per-atom λ from `fix backmap`'s per-atom array and compute weights:
- **AT weight**: `w_AT = λ_i × λ_j` (for uniform λ: `λ²`)
- **CG weight**: `w_CG = 1 − λ_i × λ_j` (for uniform λ: `1 − λ²`)

This ensures `w_AT + w_CG = 1` always — a physically consistent interpolation where total effective force smoothly transitions from CG to AT.

**Rationale**: The ESPResSo++ implementation applies lambda weighting inside the force computation templates (`VerletListAdressInteractionTemplate` for pairs, `FixedPairListAdressInteractionTemplate` for bonds/angles/dihedrals). All use the same weighting scheme: `w12 = λ₁×λ₂` for AT and `1−w12` for CG. Doing this at the style level in LAMMPS is the natural equivalent — it keeps force/energy/virial computation consistent, avoids post-hoc corrections, and integrates with LAMMPS's reporting infrastructure.

**Example input script** (polymer backmapping):
```
# Non-bonded interactions
pair_style backmap 12.0 lj/cut/coul/cut 12.0 9.0  9.5 table linear 1000
pair_coeff 1 1 atomistic ...    # OW-OW (weighted by λ²)
pair_coeff 3 3 cg CG.table WCG  # WCG-WCG (weighted by 1−λ²)
pair_coeff 1 3 none              # no cross-type interaction

# Bonded interactions — uses hybrid to mix static and lambda-weighted
bond_style hybrid harmonic backmap/harmonic backmap/table
bond_coeff 1 harmonic 500.0 1.0                   # intra-CG bond (static)
bond_coeff 2 backmap/harmonic at 500.0 1.0         # cross-CG AT bond (×λ²)
bond_coeff 3 backmap/table cg table_b1.table ENTRY # cross-CG CG bond (×(1−λ²))

angle_style hybrid harmonic backmap/harmonic
angle_coeff 1 harmonic 50.0 109.47                 # intra-CG angle (static)
angle_coeff 2 backmap/harmonic at 50.0 120.0        # cross-CG AT angle (×λ²)

dihedral_style hybrid backmap/ryckaert backmap/table
dihedral_coeff 1 backmap/ryckaert at c0 c1 c2 ...  # cross-CG AT dihedral (×λ²)
dihedral_coeff 2 backmap/table cg table_d1.table ENTRY # cross-CG CG dihedral (×(1−λ²))
```

Each `backmap/*` style takes an `at` or `cg` keyword that determines the weighting direction.

**Alternatives considered**:
- *Standard `pair_style hybrid` with force scaling in `post_force()`*: Cannot distinguish bonded from non-bonded contributions in the accumulated force array. Intra-CG bonds must remain at full strength while cross-CG bonds are scaled — impossible to achieve with post-hoc scaling.
- *Fix computes all cross-bonded forces directly*: Would bypass LAMMPS's bond/angle/dihedral infrastructure, losing energy tallying, virial computation, and restart support. Also duplicates existing LAMMPS force computation code.
- *`fix adapt` to dynamically modify coefficients*: Works for harmonic bonds (scale k by λ²) but not for tabulated potentials. Fragile and limited.

### D1a: Three categories of bonded interactions

**Decision**: Bonded interactions (bonds, angles, dihedrals, 1-4 pairs) are classified into three categories, matching the ESPResSo++ `[ bonds ]` / `[ cross_bonds ]` distinction:

| Category | Topology section | Lambda weighting | LAMMPS style |
|---|---|---|---|
| **Intra-CG bead** | `[ bonds ]`, `[ angles ]`, etc. | None (static, always full strength) | Standard `bond_style harmonic`, etc. |
| **Cross-CG, AT level** | `[ cross_bonds ]` (AT atom entries) | `w = λ₁ × λ₂` | `bond_style backmap/harmonic at ...` |
| **Cross-CG, CG level** | `[ cross_bonds ]` (CG entries, marked `cg_bonded`) | `w = 1 − λ₁ × λ₂` | `bond_style backmap/table cg ...` |

**Rationale**: In ESPResSo++, the `setBondedInteractions()` function routes bonds based on the `cross_bonds` flag:
- `cross_bonds=False` → `is_cg=None` → static `FixedPairListHarmonic` (no lambda weighting)
- `cross_bonds=True` + AT atoms → `is_cg=False` → `FixedPairListAdressHarmonic` with AT weighting
- `cross_bonds=True` + CG atoms → `is_cg=True` → `FixedPairListAdressHarmonic` with CG weighting

Physical meaning:
- **Intra-CG bonds** hold atoms together *within* a CG bead (e.g., O-H bond in water where all 3 atoms map to one WCG bead). These must always be active to maintain the bead's internal structure.
- **Cross-CG AT bonds** connect atoms that belong to *different* CG beads (e.g., backbone bond in a polymer connecting monomers in adjacent CG beads). These represent the atomistic connectivity between beads and scale in with the AT resolution.
- **Cross-CG CG bonds** are the coarse-grained bonds between CG particles (e.g., CG-CG bond in a polymer chain, often tabulated). These scale out as atomistic resolution increases.

The same pattern extends to `[ cross_angles ]`, `[ cross_dihedrals ]`, and `[ cross_pairs ]`.

### D2: Single fix handles λ ramp, CG force distribution, and COM tracking

**Decision**: One fix (`fix backmap`) implements: λ ramp, CG force distribution to AT atoms, COM tracking, and CG zeroing. Lambda-weighted force computation is delegated to the custom pair/bond/angle/dihedral styles (D1).

**Rationale**: These operations are tightly coupled and must execute in a specific order within each timestep. The fix does NOT scale forces by λ — that responsibility belongs to the interaction styles, which apply weighting at force-computation time. The fix's role is managing the backmapping state (λ) and the CG virtual-site mechanics (COM tracking, force distribution).

**LAMMPS integration points**:

```
Timestep flow:
                                           fix backmap hooks:
  ┌─────────────────────────┐
  │ initial_integrate()     │◄──── zero CG forces/velocities
  │ (fix nve on AT group)   │      (prevent CG from being integrated
  │                         │       by any other fix)
  ├─────────────────────────┤
  │ neighbor list build     │
  │ (if needed)             │
  ├─────────────────────────┤
  │ pair_style backmap       │      AT pair forces weighted by λ²
  │ + bond_style hybrid      │      CG pair forces weighted by (1−λ²)
  │   (static + backmap/*)   │      intra-CG bonds: full strength
  │ + angle_style, etc.      │      cross-CG bonds: weighted by λ²/(1−λ²)
  │ → lambda-weighted        │
  │   forces on all atoms    │
  ├─────────────────────────┤
  │ post_force()            │◄──── 1. distribute CG forces to AT atoms:
  │                         │         f[AT_i] += (m_i/m_CG) · f[CG]
  │                         │      2. zero f[CG]
  ├─────────────────────────┤
  │ final_integrate()       │
  │ (fix nve on AT group)   │
  ├─────────────────────────┤
  │ end_of_step()           │◄──── 1. CG.pos = COM(AT atoms)
  │                         │      2. λ += rate (clamped to 1.0)
  └─────────────────────────┘
```

Note: CG forces arriving at `post_force()` are already weighted by `(1−λ²)` from the interaction styles. The fix distributes these pre-weighted CG forces to AT atoms proportional to mass ratio, then zeros CG forces so CG atoms are not integrated.

### D3: CG-AT mapping via molecule IDs

**Decision**: Use LAMMPS's built-in `atom->molecule` array to identify which AT atoms belong to which CG particle. The CG atom is identified by its type (`cg_type` parameter).

**Rationale**: LAMMPS data files already assign molecule IDs. For water2, all 4 atoms (OW, HW1, HW2, WCG) in one water molecule share the same molecule ID. The fix scans once at init to build the mapping.

**Data structure** (built in `init()`):
```cpp
struct MoleculeMap {
    int cg_atom;              // local index of CG atom
    std::vector<int> at_atoms; // local indices of AT atoms
    double cg_mass;           // sum of AT masses
};
std::map<tagint, MoleculeMap> mol_map;  // keyed by molecule ID
```

Rebuilt on neighbor list rebuilds (when atoms migrate between processors).

**Alternatives considered**:
- *Explicit tuple file*: Would mirror ESPResSo++'s `FixedVSList.addTuples()`. More flexible but adds an extra input file and parsing code. Molecule IDs are simpler and already standard in LAMMPS.
- *Group-based*: Define AT and CG groups. Insufficient — need per-molecule mapping, not just per-type.

### D4: CG atoms excluded from integration via group

**Decision**: The user applies integration fixes (`fix nve`, `fix langevin`) only to the AT atom group. CG atoms are not integrated by any standard fix. The backmap fix manages CG positions directly.

```
group at_atoms type 1 2          # OW, H
group cg_atoms type 3            # WCG

fix integrate at_atoms nve       # only AT atoms
fix thermo at_atoms langevin ... # only AT atoms
fix bm all backmap cg_type 3 rate 0.0001
```

**Rationale**: Cleanest separation of concerns. The fix doesn't need to fight against or undo another fix's integration of CG atoms. LAMMPS's group mechanism is designed exactly for this.

The fix's `initial_integrate()` hook zeros CG velocities as a safety net (in case the user accidentally includes CG atoms in an integration group).

### D5: Lambda stored as fix-owned per-atom array

**Decision**: Lambda is stored as a per-atom array owned by the fix, accessible via `fix->extract()` for diagnostics. Not exposed as a general per-atom property.

**Rationale**: Only the fix needs λ during simulation. Output for analysis uses `compute property/atom` or custom dump output via the fix's vector interface.

Lambda is written to restart files so simulations can be continued.

### D6: YAML settings file as primary input (modernized from bakery XML)

**Decision**: The Python input generator reads a YAML settings file as its single source of truth. This file (modernized from bakery's VOTCA-inspired XML format) describes the CG-AT mapping, cross-CG interactions, source file references, simulation parameters, and output configuration. GROMACS files are still used as source data for AT/CG topology and coordinates, but the generator does not depend on AdResSLab or its Python 2 parsers.

**Rationale**: The bakery XML settings file is the actual source of truth in the existing workflow — GROMACS hybrid topology files (`hyb_topol.top` with `[ cross_bonds ]`) are generated from it, not the other way around. Using the settings file directly eliminates the intermediate GROMACS hybrid generation step and removes the dependency on AdResSLab's Python 2 codebase. YAML is chosen over XML for readability and native Python support (PyYAML + Pydantic validation).

**Architecture**:
```
┌───────────────────────────────────────────────────────────────┐
│ backmap-prep                                                  │
│                                                               │
│  ┌──────────────────────┐                                     │
│  │ settings.yaml        │─── CG-AT mapping, cross bonds,     │
│  │ (primary input)      │    simulation params, source refs   │
│  └──────────┬───────────┘                                     │
│             │                                                 │
│             ▼                                                 │
│  ┌──────────────────────┐   ┌────────────────────────┐        │
│  │ Source file parsers   │   │ LAMMPS writers         │        │
│  │ (new, Python 3)       │   │                        │        │
│  │                       │   │ write_lammps_data()    │        │
│  │ parse_gro()           │──▶│ write_lammps_input()   │        │
│  │ parse_top()           │   │ convert_table()        │        │
│  │ (Phase 4: parse_pdb)  │   │ unit_conversion        │        │
│  └───────────────────────┘   └────────────────────────┘        │
│                                                               │
│  Dependencies: PyYAML, Pydantic v2                            │
└───────────────────────────────────────────────────────────────┘
```

**YAML top-level structure**:
```yaml
molecules:          # CG molecule definitions and AT-CG bead mapping
cg_system:          # CG configuration files (coordinates + topology)
cross_interactions:  # Cross-CG bonds/angles/dihedrals (lambda-weighted)
simulation:         # Backmapping parameters (alpha, temperature, cutoffs, ...)
output:             # Output file configuration
```

The settings format supports all features from the bakery XML (reactive networks, degree-dependent beads, charge management) but is implemented in phases: MVP handles linear molecules, later phases add reactive network support.

**Unit conversion** (GROMACS → LAMMPS real units):

| Quantity    | Factor          | Example                   |
|-------------|-----------------|---------------------------|
| Distance    | ×10             | 0.316 nm → 3.16 Å        |
| Energy      | ×0.239006       | 0.65 kJ/mol → 0.155 kcal/mol |
| Force       | ×0.0239006      | kJ/mol/nm → kcal/mol/Å   |
| Charge      | ×1 (same)       | elementary charges        |
| Mass        | ×1 (same)       | g/mol                     |
| Time        | ×1000           | 0.001 ps → 1.0 fs        |

### D7: COM calculation uses minimum image convention

**Decision**: COM position update uses LAMMPS's `domain->minimum_image()` to handle molecules near periodic boundaries, following the same approach as ESPResSo++'s `updateVS()`.

**Algorithm** (from ESPResSo++ `VelocityVerletHybrid::updateVS`):
```
for each molecule:
    com_delta = (0, 0, 0)
    for each AT atom i in molecule:
        dr = minimum_image(r_AT_i - r_CG)
        com_delta += m_i * dr
    com_delta /= m_CG
    r_CG += com_delta
```

The CG position is updated relative to its current position using minimum-image vectors from each AT atom. This avoids issues with atoms on opposite sides of the periodic box.

## Risks / Trade-offs

**[Risk] Molecule migration across MPI boundaries** → The molecule map must be rebuilt when atoms migrate between processors. The fix registers for the `post_exchange` callback and rebuilds the map when LAMMPS performs domain decomposition. Molecules split across processors require ghost atom communication for COM calculation — handled by setting `comm_forward` and `comm_reverse` in the fix.

**[Risk] Energy is not conserved during backmapping** → Expected and inherent to the method (the Hamiltonian changes over time). The Langevin thermostat handles temperature control. This should be documented clearly so users don't interpret energy drift as a bug.

**[Risk] CG mass consistency** → CG particle mass must equal the sum of its AT atom masses. The Python generator ensures this, but a runtime check in `fix backmap::init()` validates it and errors if inconsistent.

**[Risk] Multiple custom styles increase code surface** → Creating custom pair/bond/angle/dihedral styles means more C++ code to maintain than a single-fix approach. Mitigation: all `backmap/*` styles share a common lambda-access helper class that reads from `fix backmap`'s per-atom array. Each style adds only the potential-specific force computation + weighting logic on top of this shared base.

**[Trade-off] One custom style per potential type** → Unlike ESPResSo++'s template approach (where `FixedPairListAdressInteractionTemplate<Potential>` wraps any potential), LAMMPS's style architecture requires a concrete style per potential type (e.g., `backmap/harmonic`, `backmap/table`). In practice, only a small number of styles are needed per system (typically: harmonic bonds, harmonic angles, tabulated CG bonds, Ryckaert-Bellemans or tabulated dihedrals).

**[Trade-off] Python 2 → Python 3 migration for AdResSLab modules** → AdResSLab's `gromacs_topology.py` is Python 2 code (print statements, `.iteritems()`, etc.). The generator will need a ported version or a compatibility layer. This is straightforward (mostly syntax changes).

## Open Questions

1. **Rate vs. total steps**: Should the fix accept `rate` (λ increment per step) or `nsteps` (total steps for 0→1 transition, then rate = 1/nsteps)? The latter is more user-friendly. Could support both.

2. **Lambda schedule**: Linear ramp (`λ += rate`) is what ESPResSo++ uses. Should we support non-linear schedules (e.g., sigmoid, staged) in the initial version or defer?

3. ~~**Bonded force weighting**~~ → **Resolved**: Cross-CG bonded interactions (bonds, angles, dihedrals, 1-4 pairs) ARE weighted by lambda, matching the ESPResSo++ implementation exactly. Intra-CG bonded interactions remain at full strength. See D1a for the three categories. The weighting scheme is identical for bonded and non-bonded: `w_AT = λ₁×λ₂`, `w_CG = 1 − λ₁×λ₂`.

4. **Package location**: Should this be a new repo, or a subdirectory of `lammps-user-adress`? A new repo is cleaner but adds maintenance overhead.

5. **LAMMPS version compatibility**: Which minimum LAMMPS version to target? The stable 2Aug2023 release or newer?
