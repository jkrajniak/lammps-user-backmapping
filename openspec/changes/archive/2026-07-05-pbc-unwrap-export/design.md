# Design — molecule-aware PBC unwrap

## End-user model

Hybrid network systems export like ESPResSo++ unfolded coordinates: each backmap
fragment (one CG bead + its AT atoms, one `mol_id`) forms a contiguous Cartesian
cluster. Cross bonds between fragments are satisfied by **translating entire
molecules** by integer box shifts, not by moving individual atoms.

## Algorithm

1. **Intra-mol unwrap:** BFS from CG root on bonds where both endpoints share
   `mol_id`; `_place_near_reference` per edge (existing min-image logic).
2. **Inter-mol shortening:** For bonds spanning different `mol_id`, translate all
   atoms in the higher `mol_id` fragment so the bonded partner sits in the nearest
   periodic image of its anchor. Repeat up to 32 passes until convergence.
3. **Export:** Assign bond-tree image flags and write `(x, y, z, ix, iy, iz)`.
   Generated inputs include `reset_atoms image all` and a communication cutoff
   sized for network bond extents (≥ 115 Å for rim135-class systems).

## Validation

`validate_bond_geometry` SHALL reject systems when min-image max bond exceeds
`max(20 Å, 2 × interaction_cutoff)`. When image flags are assigned, Euclidean
folded-coordinate spans are not gated (network chord bonds).

Rim135 integration test SHALL assert Euclidean max bond < 20 Å after build.

## Rationale vs image flags

LAMMPS image flags require consistent unwrapped positions across all bonded paths;
network chord bonds (~69 in rim135) make global flag assignment fragile. Molecule
translation preserves intra-fragment geometry and matches the ES++ unfolded export
pattern end users already understand from GROMACS coordinates.
