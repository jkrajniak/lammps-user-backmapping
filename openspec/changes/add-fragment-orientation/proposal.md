# Add fragment orientation (`prep.fragment_orientation`)

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
- Documentation page with a figure (`docs/fragment-orientation.md`), settings reference row,
  changelog entry.

## Out of scope

- A restraint that keeps double-bond isomers during the lambda ramp (separate question, see the
  experiment note).
- Chirality that is not present in the source molecule.

## Tests

`python/tests/test_fragment_orientation.py`: rotation helper (known rotation, proper rotation for
a mirrored target, single and collinear vectors), schema default and validation, and an end-to-end
build of the dodecane example in both modes (beads fixed, fragments rigid, orientations differ).
The default path is covered by the existing byte-parity tests.
