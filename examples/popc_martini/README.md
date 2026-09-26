# MARTINI 3 POPC bilayer -> Slipids POPC + TIP3P

A MARTINI 3 lipid bilayer in MARTINI water backmapped to an atomistic
AMBER-family lipid force field. It shows that the CG side can be a standard,
transferable CG force field (not only a structure-based one) and that solvent
beads spanning several molecules (one W bead = four waters) are handled.

- CG: MARTINI 3.0.0, 128 POPC (12 beads) + 1752 W; nrexcl 1.
- AT: Slipids POPC (Urey-Bradley angles, multi-term periodic dihedrals,
  harmonic impropers) + flexible TIP3P; nrexcl 3; 38 176 atoms.
- Mapping: CG2AT `martini_3-0_slipids` fragments (every atom in one bead).

```bash
bash fetch_sources.sh                 # third-party force-field files, SHA-256 checked
uv run python make_cg_topology.py     # flattened CG topology
uv run backmap-prep build settings.yaml
```

The CG pair interactions are generated tables reproducing the MARTINI
settings in GROMACS (LJ potential-shift at 1.1 nm, reaction field,
epsilon_r 15); CG G96 angles are generated tables too. CG pairs exclude only
bonded neighbours while AT pairs exclude up to 1-4 (`exclusion_nrexcl_cg`).
During the lambda ramp the AT Coulomb interaction is a plain cut-off (1.4 nm);
continue the backmapped system with PME.

See `SOURCES.md` for the origin and licence of the third-party files.
