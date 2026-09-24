# Unify the hybrid builder into one engine

## Why

`backmap-prep` builds hybrid systems with two independent engines, selected by
`prep.engine: linear | network`. The method is one method; there is no physical
distinction between backmapping a linear melt and a network. Two builders let the
same input produce different systems, and a defect in one is invisible from the
other. The COMPHY-D-26-00942 revision analysis (research notes
`experiments/20260923_review-energetics-dodecane-pe.md` and
`experiments/20260924_pe75-melt-rerun.md`) found that the linear engine, which
built every dodecane and polyethylene example in the paper, is wrong in ways the
network engine is not:

1. **Fragment placement.** Each AT atom is placed at `bead + (template_atom -
   template_atoms[0])`, i.e. relative to the first atom of the whole template
   molecule, not with the fragment COM on its bead. For PE75 the bead-to-COM
   distance is 34 A (median); for dodecane 6 A. At the first MD step
   `fix backmap` moves every bead onto its fragment COM and the CG configuration
   collapses (CG bond energy 1.2e3 -> 4.9e5 kcal/mol in 0.2 fs). Present since
   the first commit (`f313f2b`).
2. **CG angles and dihedrals are dropped.** The linear PE75 data file contains no
   CG angles and no CG dihedrals, although the settings define them.
3. **Cross-bead AT dihedral parameters are shifted.** Settings `params` strings
   start with the GROMACS function code (`3 C0..C5`); the linear engine reads it
   as C0 (`3/4.184 = 0.717018`) and drops C5.

The network engine (vendored bakery) places fragments with their mass-weighted
COM on the bead, emits every bonded term defined in the settings, and parses the
function code. It is also the engine with byte-level parity against bakery's
published hybrid files (rim135, PET). It is the right basis for a single engine.

Separately, the committed dodecane and PE production inputs carry bond and angle
force constants 2x the GROMACS values (generated before the k/2 conversion fix of
2026-08-06 and never regenerated). This change regenerates every example from
the single engine, which removes those stale inputs.

## What Changes

- **BREAKING (settings):** `prep.engine` is removed. A value of `linear` or
  `network` is accepted with a deprecation warning and ignored. `hybrid` becomes
  optional with defaults (`hyb_conf.gro`, `hyb_topol.top`, molecule type `HYB`,
  exclusion 3). Beads defined with a plain `atoms:` list map at any degree.
- The network engine becomes the only builder for `build`, `rebuild`, and
  `cg-only`. `builder.build_system` and the linear-only code paths are deleted.
- Features that only the linear engine had are ported: LAMMPS-format
  `cg_system` and `molecules[].source` (converted to GROMACS-equivalent inputs
  before the bakery build, with an exact unit round trip), `CGPositionOverride`
  rebuilds.
- Network-engine defects exposed by linear systems are fixed: united-atom LJ
  parameters from `[ atomtypes ]` sigma/epsilon (currently written as 0), the
  communication cutoff (85 A for PE75), and already-qualified bead atom names
  (`1:PE:C1` was prefixed again and the fragment silently skipped).
- Every example is regenerated. Hand-maintained production inputs that restate
  the force field (`in.pe_robust`, `in.pe_at*`, `in.dodecane*`, ...) take their
  force-field block from the generated input instead of restating it.
- New validation: for each example, per-term energies of the generated LAMMPS
  system at lambda = 1 must match GROMACS (`gmx energy` on the hybrid topology's
  AT part, same coordinates) within tolerance; every bead must sit on its
  fragment COM.

## Impact

- Affected specs: `backmap-input-generator`.
- Affected code: `python/src/backmap_prep/{builder,cli,schema,writers}.py`,
  `network/*`, tests (`test_builder.py` and friends), all `examples/*`,
  docs (settings reference, tutorial, getting started), CHANGELOG, README.
- Scientific impact: the dodecane and PE Tier B/C results reported in the
  submitted manuscript were produced with the defects above. They must be rerun
  and the manuscript and response letter must say so.
