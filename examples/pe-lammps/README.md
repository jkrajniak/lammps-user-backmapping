# PE — fully LAMMPS-native input (parity example)

This example is a parity check for using `format: lammps` on **both**
`cg_system` and `molecules[].source` at once — a fully GROMACS-free
`backmap-prep` run. See [`../pe/README.md`](../pe/README.md) for the PE
(polyethylene) system description (100-carbon chain, 50 CG beads,
tabulated CG interactions, OPLS united-atom AT topology).

## What's different from `../pe/`

Both source formats, in `settings.yaml`:

```yaml
# ../pe/settings.yaml
cg_system: {coordinates: cg_conf.gro, topology: topol_cg.top, format: gromacs}
molecules[0].source: {coordinates: pe_single.gro, topology: topol_aa.top}
```

```yaml
# this example
cg_system: {format: lammps, data: pe_cg.data}
molecules[0].source: {format: lammps, data: pe_at.data, input_script: in.pe_at}
```

`cross_interactions`, `simulation`, and the tabulated CG potentials are
otherwise unchanged — those are configured the same way regardless of
either format (see `openspec/changes/lammps-native-cg-input/design.md` and
`.../lammps-native-at-fragment-input/design.md`).

## Fixture provenance

Both `pe_cg.data` and `pe_at.data`/`in.pe_at` are **converted directly
from `../pe/`'s own GROMACS source files** (`cg_conf.gro`/`topol_cg.top`
and `pe_single.gro`/`topol_aa.top` respectively) using `backmap_prep`'s
own GROMACS parser and unit-conversion functions — *not* extracted from a
built hybrid. Two reasons this matters, both found the hard way while
verifying this example on real LAMMPS (not just at the Python level):

- **AT fragments are a rigid single-molecule template** (relative atom
  offsets from atom 1), unlike the CG side. Extracting one from an
  already-backmapped hybrid instance captures that instance's actual
  (non-ideal) placement, not the idealized template — silently wrong
  geometry.
- **CG positions matter too, for a subtler reason**: `cg_conf.gro`
  (a real GROMACS trajectory) leaves some CG beads sitting just outside
  `[0, box)`, which is normal. A `pe_cg.data` derived by extracting from
  an *already-wrapped* hybrid has everything pre-wrapped inside
  `[0, box)` instead — physically equivalent, but it changes which
  raw/pre-wrap representative feeds into `network/pbc.py`'s per-molecule
  PBC unwrap, which can pick a different (though still self-consistent)
  overall image-flag offset per molecule. Converting directly from
  `cg_conf.gro` avoids the question entirely.

## The numeric-ID convention

Bead `type` and `atoms` entries reference LAMMPS numeric type/atom IDs
(as strings) rather than symbolic GROMACS names — `type: "1"` instead of
`type: A`, `atoms: ["1:PE:3"]` instead of `["1:PE:C3"]` — since neither a
`data` file's `Masses` section nor its `Atoms` section carries a symbolic
name. `cross_interactions` atom references (not bead references) are
rewritten the same way, e.g. `PE:C2` → `PE:2`; bead references like
`PE:A1`/`PE:B12` are untouched.

## Verifying parity

Building both `../pe/settings.yaml` and this example's `settings.yaml`
produces **physically identical** hybrid output: all 1500 atoms match
exactly on molecule ID, type, charge, position, *and* PBC image flags,
and `in.pe`'s coefficients match exactly (aside from the expected
`table_A_A.table` vs `table_1_1.table` filename difference from the
numeric-ID convention). This was verified both at the Python/text level
and by actually running both builds to completion on a real LAMMPS
binary with the validated 3-phase "robust" protocol (minimize →
`nve/limit` ramp → staged NVT — see `examples/pe/large/in.pe_robust`;
the default auto-generated `in.pe`'s Langevin protocol is known-fragile
even for correct topology, per `research/notebook/2026-07-09_pe-tier-b-pass.md`).

### A trail of two wrong "it's benign" conclusions, corrected in place

Getting to that clean result took three attempts, each corrected here
rather than silently rewritten, because each one taught something real:

1. **First claim: "pre-existing bug in `network/pbc.py`'s image-flag
   assignment, unrelated to this PR."** Wrong in its scope, though the
   underlying observation (LAMMPS logs `WARNING: Inconsistent image
   flags` for this system) was real — that warning is still present in
   the passing runs below and is confirmed benign. The conclusion "not
   my bug" was reached by testing with a stale, non-matching pair-table
   fixture and never actually finding what made the run crash.
2. **Second claim: the crash was caused by a genuine 34-vs-35
   inconsistent-bond difference, root-caused via unwrapped bond-length
   arithmetic.** This was a real, fixable difference (closed by the
   `pe_cg.data`/`pe_at.data` provenance fix above, which now gives 0/1500
   mismatches everywhere, including image flags) — but fixing it alone
   did **not** stop the crash. The image-flag story was a real, separate
   thing that turned out not to be what was killing the simulation.
3. **Actual root cause, found by testing the untouched pristine pair
   table directly and working backward from what changed**: `table_A_A.table`'s
   real IBI data starts at `r_min = 0.02 Å`. `builder.py`'s CG
   pair-table lookup (added for `cg_system.format: lammps`, but applied
   identically regardless of format) prefers a `.table` file over `.xvg`
   when both exist — and `table_converter.py` used to *re-extend* a
   `.table` source with a repulsive wall down to `r_floor = 1e-4 Å` even
   though it's already native LAMMPS format, exactly like every other
   `.table`-sourced table already skips. The extrapolated wall value
   for this large an `r_min`/`r_floor` gap reached ~1.4e31 kcal/mol —
   large enough to corrupt LAMMPS's `pair_style table` cubic-spline fit
   near that region, not just at the single extension point. **Capping
   the wall energy at a sane ceiling (1e6 kcal/mol) did not fix the
   crash either** — confirming the problem wasn't the wall's magnitude
   but the act of inserting one far-spaced extension point into an
   otherwise uniformly-spaced table at all. The actual fix: `.table`
   pair sources are now a pure passthrough, matching how bond/angle/
   dihedral `.table` sources already worked. Verified by testing the
   untouched pristine table directly (worked), then confirming the fix
   in code reproduces that exact byte-identical untouched file.

Both the image-flag warning and the wall-extrapolation bug were
pre-existing and equally reachable via the GROMACS path once it also
goes through the same `.table`-preferring lookup — neither is specific
to LAMMPS-native CG or AT-fragment input. This example was simply the
first time anyone ran the (small) `examples/pe/` all the way through a
real, validated LAMMPS protocol.

## AT-fragment dihedral coverage note

PE's beads are 2 atoms each, same as dodecane. A 4-atom dihedral spanning
consecutive backbone atoms is therefore *never* fully contained within one
bead for either system — the LAMMPS-native AT-fragment reader's
`dihedral_coeff`/`ryckaert` parsing path exists and is correct (verified
independently in `tests/test_lammps_data_parser.py` and
`tests/test_lammps_script_parser.py` with a synthetic fixture), but no
example in this repository's linear-engine set exercises it end-to-end,
because none currently maps a bead to 3+ atoms. This is a pre-existing
property of the example suite (the same is true of the GROMACS path's own
intra-bead dihedral code), not a gap introduced by this feature.
