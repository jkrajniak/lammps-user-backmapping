# Design — Phase 3 network backmapping

## Context

Linear polymer migration (Phase 1/2) is complete for dodecane, PE variants, and
non-reactive melamine. Bakery network examples use XML features not yet in
`settings.yaml`:

| Feature | Bakery location | Current backmap-prep |
|---------|-----------------|----------------------|
| `atoms_by_degree` | rim135, epoxy XML | Not implemented |
| Multi-molecule types | EPO + IPD + HDD in one hybrid | Single molecule type |
| Active sites | Cross-link placement | Not implemented |
| Charge transfer | Reaction rules in structures.py | Not implemented |
| Exclusion lists | `exclusion_hyb_topol.list` | Partial / manual |

Research spec `jcc-2017-validation` defines acceptance before manuscript claims.

## Goals

- Implement minimum generator surface for **epoxy (RIM135)** Tier A parity
- Enable LAMMPS backmapping smoke (Tier B) on cured epoxy topology
- Path to Tier C RDF vs ESPResSo++ reference (±5% peak positions)

## Non-goals

- Live chemical reaction during λ ramp
- PET / hyperbranched in the first merge (follow-on changes)
- ESPResSo++ engine maintenance

## Schema extensions (proposed)

```yaml
molecules:
  - name: EPO
    beads:
      - name: bead_A
        atoms_by_degree:
          0: [atoms for degree 0]
          1: [atoms for degree 1]
          2: [atoms for degree 2]
    active_sites: [...]
    charge_rules: [...]
```

Exact YAML shape to be refined against `bakery/tests/rim135/settings.xml` during
task 1.x gap analysis.

**Update (2026-06-22):** Settings v2 **implemented** for rim135. See
[settings-v2-rim135-plan.md](./settings-v2-rim135-plan.md),
[settings-format-v2.md](./settings-format-v2.md), `network/v2_loader.py`.
`backmap-prep build-hybrid examples/epoxy/settings.v2.yaml` passes Tier A parity.
Next: wire unified `build` → LAMMPS; Tier B smoke on VM.

## Migration order

1. **Epoxy (rim135)** — flagship JCC paper system; bakery test exists
2. **MF polymerized** — `network_backmapping/mf/` (depends on melamine base)
3. **PET** — `network_backmapping/pete/` (multi-species tables)
4. **AB2/ABx** — hyperbranched stretch

## Validation linkage

| Tier | Epoxy criterion | Source spec |
|------|-----------------|-------------|
| A | Hybrid topology vs `ref_hyb_*` | migration-parity (research) |
| B | λ ramp completes | jcc-2017-validation |
| C | RDF peaks within 5% of ESPResSo++ | jcc-2017-validation |
| — | Charge conserved ±1e-6 e | jcc-2017-validation |

## Risks

| Risk | Mitigation |
|------|------------|
| XML complexity | Start from rim135 test case only; manual YAML |
| Charge bugs | Unit tests on charge sum before LAMMPS |
| Large topology | Reduced system size for CI; full size for paper |

## Research cross-reference

When this change merges, update:

- `user-backmapping/research/openspec/.../tasks.md` Phase 3 checkboxes
- `user-backmapping/research/STATUS.md` network row
