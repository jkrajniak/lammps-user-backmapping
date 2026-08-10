# Design — crosslinked melamine (MF) network example

## Context

`research/experiments/20260802_melamine-bakery-rerun.md` (2026-08-06 update)
found bakery's actual melamine reference is a real, ~45%-crosslinked covalent
network (675 inter-molecular bonds across 1500 CG beads,
`bakery_full_unused/examples/network_backmapping/mf/backmapping/cg_topol.top`),
while `examples/melamine/` simulates 500 fully independent monomers. Comparing
the two is very likely the dominant reason L2(g) has never passed in the
melamine campaign, independent of force-field tuning.

`phase3-network-backmapping` already built and proved (on epoxy/RIM135 and PET)
a generator surface for exactly this class of system: degree-dependent bead
mapping, active-site crosslink placement, multi-molecule hybrid systems. This
change applies that proven surface to MF, per that change's own design.md
migration order (item 2, right after epoxy).

**Confirmed by source inspection** (not assumed): `network/api.py` imports
`from .bakery.structures import BackmapperSettings2` -- bakery's own
topology-building Python module (`structures.py`, pure Python, not the
ESPResSo++ C++ MD engine) is already vendored into this repo at
`backmap_prep/network/bakery/`. `network/v2_loader.py`'s `_emit_bead_extras`
already emits bakery-XML-format `<remove active_site="...">` elements matching
exactly the semantics MF's own `settings.xml` uses (`beads degree="2"` /
`degree="3"`, `active_site`, `remove`). `HybridBuildResult.missing_definitions_path`
already tracks a `missing_definitions.txt`-style output, the same mechanism
bakery's own MF example uses to flag unparameterized angle/dihedral terms. This
means the hard part (active-site atom removal, degree-dependent AT fragment
selection) is bakery's own proven code, already integrated -- not new code we
need to write.

## Goals

- Import bakery's exact 675-bond MF crosslink network as a fixed, static
  topology (no live reaction dynamics).
- Produce a stable eq/ramp/production LAMMPS run on the resulting crosslinked
  AT structure.
- Source real OPLS-AA parameters for the crosslink-site C-O-C angle / C-C-O-C
  dihedral terms bakery's own reference left unparameterized.
- Compute RDFs against the same reference `.xvg` files, for the first time as
  a fair, like-for-like comparison.

## Non-goals

(See proposal.md's Non-goals -- repeated here for the parts that shape the
architecture specifically.)

- No live/reactive bond formation -- the network is static, imported once.
- No changes to `examples/melamine/` or its capability spec.
- No ESPResSo++ build -- only the already-vendored pure-Python `bakery.structures`
  module is used, not the C++ engine.
- No general "arbitrary degree of cure" generator -- this imports one specific,
  fixed network.

## Architecture / data flow

```
bakery cg_topol.top (675 "chem" bonds, 1500 "static" bonds, 500 residues)
        |
        v  [new: one-time conversion step]
crosslink pair list (molecule-instance + bead pairs), + active-site atom
identities per bakery's settings.xml (which arm's O stays as the bridging
ether O, which arm's C1/O1/H01 group is fully removed)
        |
        v  [new: examples/melamine_network/settings.v2.yaml]
molecules: MF with beads A1/A2/A3, mapping: degree-1 (unreacted, matches
examples/melamine/'s existing 9-atom bead) + degree-2 (reacted: base=1,
remove=[...], active_sites=[...]) -- same shape as epoxy's A1/A2 beads in
examples/epoxy/settings.v2.yaml, NOT a new schema feature
        |
        v  [existing, proven: network/v2_loader.py -> bakery.structures -> network/lammps_builder.py]
hybrid GRO/TOP (Tier A parity check against bakery's own hyb_conf.gro/hyb_topol.top)
        |
        v  [existing: backmap-prep build-hybrid -> LAMMPS .data + in.*]
        |
        v  [new: OPLS-AA ether angle/dihedral coeffs added to shared forcefield_dir,
        |   resolving whatever `missing_definitions.txt` reports for this system]
        |
        v  [existing protocol pattern: bakery-faithful eq -> ramp -> production,
        |   VM-only execution, tmux -- same pattern as this session's melamine work]
LAMMPS production run -> dump.at_prod
        |
        v  [existing: compare_melamine_structure.py, same ref/*.xvg files]
RDF comparison report
```

The LAMMPS-level mechanics for a crosslink bond itself need no new C++ code:
a crosslink is a real, permanent covalent bond between two AT atoms, so it is
written as a plain always-on `harmonic` bond/angle/dihedral term (matching how
same-CG-bead intra-fragment terms are already treated -- "real intra-molecular
chemistry that exists independent of the CG/AT resolution ramp",
`backmap_lambda_weights.h`), not the lambda-weighted `backmap/harmonic`
inter-bead style. RIM135/epoxy already has real EPO-IPD-HDD crosslinks and
passes Tier B/C, so this path is already exercised and proven, not new.

## Testing (tiers, matching this project's existing convention)

- **Tier A (topology parity):** generated hybrid GRO/TOP atom/bond/angle/dihedral
  counts match bakery's own `hyb_conf.gro`/`hyb_topol.top` for this system;
  total charge conserved within 1e-6 e (mirrors the `network-backmap-prep`
  charge-conservation requirement already in `phase3-network-backmapping`'s spec
  delta, applied here without needing that change's still-open charge-transfer
  task -- see Risks).
- **Tier B (stability):** eq/ramp/production completes with zero/near-zero
  dangerous neighbor rebuilds, matching the bar already met by the uncrosslinked
  example and by RIM135/PET.
- **Tier C (structure):** RDF comparison against `ref_C_C.xvg`/`ref_C_N.xvg`/
  `ref_O_H.xvg`, reported as-is (see proposal.md Success criterion -- no target
  pass count committed to).

## Risks

| Risk | Mitigation |
|------|------------|
| Translating bakery's specific 675-bond crosslink pattern into `settings.v2.yaml` at scale (not hand-authorable) | Write a small one-time conversion script parsing `cg_topol.top`'s `; chem` bonds + `settings.xml`'s active-site/remove rules directly into the v2 YAML shape (or whatever intermediate format `network/v2_loader.py` expects) -- analogous in spirit to today's `symmetrize_ring_charges.py`/`ring_angle_variance_check.py` one-off scripts, but feeding the proper pipeline instead of hand-patching a data file |
| Missing OPLS-AA parameters for crosslink-site C-O-C angle / C-C-O-C dihedral (confirmed absent even in bakery's own reference) | Source real OPLS-AA ether parameters (approved decision, not placeholders); add to the shared `forcefield_dir` alongside melamine's existing amide/imide dihedral fill-ins; verify explicit params always take priority over table lookup (already true for the existing dihedral fill-in, per `research/experiments/20260802_melamine-bakery-rerun.md`) |
| `phase3-network-backmapping` task 2.4 (charge transfer rules) is still unchecked/incomplete | Not required here: inspected bakery's own `settings.xml` for MF -- the `degree="2"` (unreacted) and `degree="3"` (reacted) bead blocks use **identical** `charge_map` values; reacted sites just remove atoms (and their charges) rather than redistributing charge among survivors. Reusing bakery's own approach (no charge transfer, just removal) sidesteps this open task entirely. Charge-sum conservation is still checked (Tier A), just not via a "transfer" mechanism. |
| MF's bakery source layout (single `at_topol.top` + XML active-site/remove rules) may differ from epoxy's (separate pre-built per-degree `.itp`/`.gro` files) in a way the v2 loader doesn't yet handle | Not yet confirmed either way -- first implementation task is reading `network/v2_loader.py`'s mapping-resolution code end-to-end against both epoxy's and (reconstructed) MF's input shapes before writing `settings.v2.yaml`. If a small loader extension is needed, scope it as a minimal, additive change (new source-entry variant), not a rewrite -- the loader already vendors bakery's own structures.py, which natively understands MF's exact XML shape, so the gap (if any) is likely only in how our YAML expresses it, not in the underlying engine. |
| Density/box equilibration gap (no barostat, flagged separately this session for the uncrosslinked example) | Out of scope for this change specifically, but worth checking whether bakery's crosslinked `cg_conf.gro` box/density looks more physically reasonable than the uncrosslinked case's ~0.644 g/cm³ before assuming it needs the same fix |

## Execution

All LAMMPS compilation/execution happens on the VM (`sc@<vm_host>`), in tmux,
per established project practice -- same as every run this session. Local work
is limited to Python/YAML/force-field authoring and result analysis.

## Research cross-reference

When a run is validated: update `research/checkpoints.md` (new melamine-network
entry, per that file's registry practice -- do not add until a run is
confirmed), and `phase3-network-backmapping/tasks.md` item 6.1 (check it off,
link to this change).
