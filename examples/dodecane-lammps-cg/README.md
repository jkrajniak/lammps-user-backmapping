# Dodecane — LAMMPS-native CG system (parity example)

This example is a parity check for `cg_system.format: lammps`, not a new
physical system. See [`../dodecane/README.md`](../dodecane/README.md) for
the dodecane system description (12-carbon linear alkane, 6 CG beads,
tabulated CG interactions, GROMOS united-atom AT topology).

## What's different from `../dodecane/`

Only the `cg_system` block of `settings.yaml`:

```yaml
# ../dodecane/settings.yaml
cg_system:
  coordinates: cg_conf.gro
  topology: topol_cg.top
  format: gromacs
```

```yaml
# this example
cg_system:
  format: lammps
  data: dodecane_cg.data
```

`dodecane_cg.data` is a standalone LAMMPS `data` file (box, `Masses`,
`Atoms # full`) supplying the same CG system as `cg_conf.gro`/`topol_cg.top`,
generated with `../dodecane/prepare_cg.py` from a freshly-built
`../dodecane/dodecane.data` hybrid output. Everything else — the AT
fragment (`dodecane_single.gro`/`topol_aa.top`), `cross_interactions`,
`simulation`, and `table_groups` — is unchanged, since CG-CG bonded terms
and nonbonded tables are configured the same way regardless of CG format
(see `openspec/changes/lammps-native-cg-input/design.md`).

## The `type: "1"`/`"2"` convention

A LAMMPS `data` file has no symbolic atom-type name — only numeric type
IDs. `dodecane_cg.data`'s `Masses` section assigns type 1 = A (29.062
g/mol), type 2 = B (28.054 g/mol). Because of this, `molecules[].beads[].type`
and `simulation.table_groups` reference `"1"`/`"2"` here instead of the
symbolic `A`/`B` used in `../dodecane/settings.yaml`. The pair-table files
are named to match (`table_1_1.table`, `table_1_2.table`, `table_2_2.table`)
rather than `table_A_A.table` etc.

## Verifying parity

Running `backmap-prep build settings.yaml` here and in `../dodecane/`
(against a `cg_conf.gro` consistent with the `dodecane.data` that
`dodecane_cg.data` was derived from) produces numerically identical hybrid
output — the only differences are the cosmetic type-name label in
comments/table filenames. This was checked by hand during implementation;
see the research experiment note for the exact command trail and diff.
