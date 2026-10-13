# Add fragment orientation and stereocenter restraint (`prep.fragment_orientation`, `simulation.stereo_restraint`)

## Why

`backmap-prep` places each atomistic fragment on its CG bead by translation only. A stereocenter
or a double bond whose atoms lie in two beads then gets its geometry from the CG configuration,
which carries no chirality. In the MARTINI 3 POPC bilayer, 71 of 128 glycerol carbons (C2) had the
handedness opposite to the independent CG2AT reference after backmapping, and 11 of 128 oleoyl
double bonds were trans after the ramp
(`research/experiments/20261009_popc-chirality-bonded.md`). The bilayer observables used so far
(area per lipid, thickness, order parameters, chain dihedrals) do not detect either.

## What changes

- New option `prep.fragment_orientation: translate | align` (default `translate`, output
  byte-identical to earlier releases).
- `align`: per residue, (1) rotate each fragment so that the template vectors between fragment
  centers of mass match the vectors between the bonded CG beads (Kabsch, proper rotation),
  (2) refine with the covalent bonds that join fragments, (3) for every stereocenter whose
  neighbors lie in two fragments, require the handedness of the source molecule and give a
  violating fragment the best-fitting orientation among a fixed number of seeded random rotations.
  Beads and the internal geometry of fragments do not change.
- New option `simulation.stereo_restraint` (kJ/mol, robust protocol, default off): `fix restrain
  dihedral` on the improper dihedral of every stereocenter of the atomistic part, from the first
  minimization until `<prefix>_hybrid.data` is written. A stereocenter is an atom with four bonded
  atomistic neighbors that are inequivalent (Morgan classes over atom type and environment); the target
  is the ideal tetrahedral improper (35.26 deg) with the sign it has in the hybrid structure.
  `stereo.lmp` and `stereo_release.lmp` are written next to the input script. Reason: with `align`
  the handedness is right at placement but 17 of 128 POPC lipids invert at the start of the ramp,
  when the bonds between fragments switch on from a nearly flat arrangement. A soft restraint
  (83.68 kJ/mol) leaves 6; stiffer restraints strain the center and a decaying restraint releases
  the same lipids (`research/experiments/20261010_popc-stereo-restraint.md`,
  `20261012_popc-aligned-restrained-tier-c.md`).
- Documentation page with a figure (`docs/fragment-orientation.md`), settings reference row,
  changelog entry.

## Out of scope

- A restraint that keeps double-bond isomers: with `align` all double bonds of the POPC example stay cis.
- A force-field improper term at the stereocenter: it is the same energy as the restraint kept on,
  and would hold the sign at the price of a strained center.
- Chirality that is not present in the source molecule.

## Tests

`python/tests/test_fragment_orientation.py`: rotation helper (known rotation, proper rotation for
a mirrored target, single and collinear vectors), schema default and validation, and an end-to-end
build of the dodecane example in both modes (beads fixed, fragments rigid, orientations differ).
The default path is covered by the existing byte-parity tests.
`python/tests/test_stereo_restraint.py`: stereocenter detection (handedness sign, hydrogen left out,
equivalent neighbors, same-type neighbors told apart by environment), the written files, default off.
