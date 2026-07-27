# Tasks — Phase 3 network backmapping

Linked research plan: `user-backmapping/research/openspec/changes/bakery-lammps-migration-validation/`
Settings v2 plan: [settings-v2-rim135-plan.md](./settings-v2-rim135-plan.md)

---

## 1. Gap analysis and schema

- [x] 1.1 Read `bakery/tests/rim135/settings.xml` and document required YAML fields
- [x] 1.1b **Settings format v2** design — `settings-format-v2.md` + `examples/epoxy/settings.v2.yaml`
- [ ] 1.2 Compare with `bakery/src/structures.py` degree-dependent logic
- [ ] 1.3 Draft YAML schema extension in `docs/settings-reference.md` (draft section)
- [ ] 1.4 Add research decision record for schema choices (`research/decisions/`)

## 2. Generator implementation

- [x] 2.0 Implement Settings v2 loader (`mapping.base` / `add` / `remove` resolution) — `network/v2_loader.py`
- [x] 2.0b Schema extensions — delta mapping, unified `source`, `angles_file`, `data_dir`
- [x] 2.0c Wire CLI + `build_hybrid_gromacs(Settings)` + `BackmapperSettings2(Element)`
- [ ] 2.1 Implement `atoms_by_degree` parsing in builder (v2 uses `mapping` list; legacy alias optional)
- [ ] 2.2 Support multiple molecule entries in one `settings.yaml` run (rim135 already multi-species in one hybrid)
- [x] 2.3 Active site placement — via v2 loader + `predefined_active_sites.txt` (Tier A parity)
- [ ] 2.4 Implement charge transfer rules with conservation check
- [x] 2.5 Unit tests — `test_v2_loader.py`; integration `test_rim135_build_hybrid_v2_parity`

## 3. Epoxy example migration

- [x] 3.1 `examples/epoxy/` layout (assets stay in `bakery/tests/rim135`)
- [x] 3.2 `settings.v2.yaml` from bakery XML reference
- [x] 3.3 Tier A parity vs `ref_hyb_conf.gro` / `ref_hyb_topol.top` (`chain_rng_seed=42`)
- [x] 3.4 Generate LAMMPS input; Tier B smoke backmapping (CG angle tables — `archive/2026-06-22-cg-angle-tables`; CG dihedrals — `archive/2026-06-22-cg-dihedrals`; local smoke: 10+10 steps on rim135)
- [x] 3.5 `run_test.sh` (uses v2); `compare-topology` in test suite

## 4. Validation and docs

- [x] 4.1 Tier C lite: C–O/C–N RDF peak metrics vs paper AA ref (4/4 PASS, Jul 2026)
- [ ] 4.2 Update `integration-testing` spec scenarios for epoxy
- [ ] 4.3 MkDocs network backmapping tutorial
- [x] 4.4 Notify research repo to pin evidence in EVIDENCE_BASE
- [x] 4.5 `examples/epoxy/README.md` — v2 as primary entry point

## 5. Settings v2 — next steps

- [ ] 5.0 Wire `backmap-prep build` for network + v2 (LAMMPS `.data` / `in.*`)
- [ ] 5.1 Optional `migrate-settings` XML → YAML tool
- [ ] 5.2 Deprecate v1 `prep.bakery_xml` bridge in docs (keep compatibility)

## 6. Follow-on (separate changes)

- [ ] 6.1 MF polymerized network (`network_backmapping/mf/`)
- [ ] 6.2 PET (`network_backmapping/pete/`)
- [ ] 6.3 Hyperbranched AB2/ABx

---

## Blockers from research Phase 1

PE linear Tier B still **FAIL** (2026-06-22) — unrelated to Phase 3 scope but
may share force-weighting fixes. See `research/notebook/2026-06-22_pe-tier-b-smoke.md`.

Rim135 Tier B blocked on LAMMPS repo deploy to the validation VM.
