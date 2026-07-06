# Phase 3 — network backmapping migration

## Why

JCC 2017/2018 (Krajniak, Zhang et al., *J. Comput. Chem.* 2018, online 2017)
extended adaptive-resolution reverse mapping to **polymer networks**: epoxy
(RIM135), trimethylol melamine polymerization, PET, and hyperbranched systems.
These require features beyond Phase 1/2 linear-polymer support:

- Degree-dependent bead definitions (`atoms_by_degree`)
- Multi-molecule-type hybrid systems (EPO, IPD, HDD, etc.)
- Active site placement and charge transfer on bond formation

The research OpenSpec `bakery-lammps-migration-validation` defines acceptance
criteria in `research/openspec/.../specs/jcc-2017-validation/spec.md`. This
code-repo change implements the generator and LAMMPS example migration blocked
on those criteria.

## What changes

- Extend `backmap-prep` YAML schema and builder for network features
- Migrate bakery network examples starting with **rim135 / epoxy**
- Add `examples/epoxy/` (and later PET, hyperbranched) with Tier A parity vs
  `bakery/tests/rim135/`
- Document degree-dependent mapping in MkDocs

## Capabilities

### New capabilities

- `network-backmap-prep` — degree-dependent beads, multi-molecule YAML, charges
- `epoxy-example` — RIM135 cured network example and validation harness

### Modified capabilities

- `backmap-input-generator` — schema extensions for `atoms_by_degree`, active sites
- `integration-testing` — epoxy integration test (currently PRECONDITION_NOT_MET)

## Non-goals

- Porting ESPResSo++ runtime (LAMMPS only)
- Full migration of all `network_backmapping/*` subprojects in one change
- Chemical reaction dynamics during backmapping (static cured topology only)

## Impact

- **Python:** `backmap_prep/` builder, parsers, writers, schema validation
- **Examples:** new `examples/epoxy/` from `bakery/examples/rim135/`
- **Docs:** settings reference, network tutorial
- **Research:** unblocks JCC validation tasks Phase 3 in research OpenSpec

## Dependencies

- Research acceptance criteria: `user-backmapping/research/openspec/changes/bakery-lammps-migration-validation/specs/jcc-2017-validation/spec.md`
- Bakery reference: `bakery/tests/rim135/`, `bakery/examples/network_backmapping/epoxy/`
- Deferred from: archived `migrate-bakery-examples` (Decision 1)

## References

- DOI: 10.1002/jcc.25129
- Bakery XML patterns: `network_backmapping/epoxy/`, `rim135/settings.xml`
