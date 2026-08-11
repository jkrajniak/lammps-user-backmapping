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

- `pe_cg.data`: derived from a freshly-built `../pe/` hybrid via
  `../dodecane/prepare_cg.py` (generic, reused as-is).
- `pe_at.data` / `in.pe_at`: **converted directly from `../pe/pe_single.gro`
  and `../pe/topol_aa.top`** using `backmap_prep`'s own GROMACS parser and
  unit-conversion functions — *not* extracted from a built hybrid. AT
  fragments are a rigid single-molecule template (relative atom offsets
  from atom 1), unlike the CG side; extracting from an already-backmapped
  hybrid instance captures that instance's actual (non-ideal) placement
  rather than the idealized template, which silently produces the wrong
  geometry. See the design doc's correction note for the full explanation.

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
exactly on molecule ID, type, charge, and position, and `in.pe`'s
coefficients match exactly (aside from the expected `table_A_A.table` vs
`table_1_1.table` filename difference from the numeric-ID convention).

The one textual difference is per-atom PBC image flags (`ix iy iz`) for
atoms whose molecule crosses a periodic boundary: 79 of 1500 atoms get a
different assignment than the GROMACS path, because the post-build
image-flag canonicalization walks the bond graph in enumeration order,
and this example's AT-fragment bond enumeration order differs slightly
from the GROMACS path's. The *wrapped* position is identical either way.

> **IMPORTANT CORRECTION:** this section originally called the image-flag
> difference "self-consistent" and "not something LAMMPS's force
> computation depends on." That was wrong, and was not verified before
> being written — see the debugging-discipline lesson below. Running this
> example on a real LAMMPS build (the project's compute VM) crashes at step 2 with
> `WARNING: Inconsistent image flags` followed by `ERROR: Bond atoms ...
> missing`. Checking the unwrapped bond lengths directly (`x + ix*box`)
> confirms 34/1480 bonds have a genuinely inconsistent image-flag pair
> (one endpoint off by one box length from what its bonded partner
> implies) — not a benign canonicalization choice.
>
> Critically: **running the unmodified `../pe/` GROMACS-path build on the
> same LAMMPS binary reproduces the identical failure** (35/1480 bad
> bonds, crashes with the same warning a few steps into the second `run`
> block). This is a pre-existing bug in `network/pbc.py`'s image-flag
> assignment (`prepare_network_coordinates`/`_assign_image_flags_by_bond_tree`),
> not something this LAMMPS-native-input feature introduced or made
> materially worse (34 vs 35 bad bonds) — it was simply never exercised by
> a full, real LAMMPS run of `examples/pe/` before. Out of scope for this
> change to fix; tracked separately. The lesson: "positions match" is not
> the same claim as "image flags are self-consistent," and the original
> text asserted the second without checking it — exactly the kind of
> unverified, convenient-sounding conclusion the project's debugging
> discipline warns against.

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
