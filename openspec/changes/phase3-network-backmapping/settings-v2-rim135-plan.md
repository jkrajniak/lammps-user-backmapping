# Plan — Settings v2 for rim135

**Status:** Tier A complete (2026-06-22)
**Linked:** [settings-format-v2.md](./settings-format-v2.md), [tasks.md](./tasks.md)

---

## Goal

Epoxy (rim135) uses the same unified path as melamine/dodecane:

```text
settings.v2.yaml → backmap-prep build-hybrid → hyb_conf.gro + hyb_topol.top → (future) LAMMPS
```

No two-stage hybrid→LAMMPS bridge. No runtime dependency on bakery XML for authoring.

---

## Done

### Schema & loader

- [x] v2 types in `python/src/backmap_prep/schema.py` (`ActiveSiteSpec`, delta mapping, unified `source`, `angles_file`, `data_dir`)
- [x] `network/v2_loader.py` — YAML → bakery-compatible XML ElementTree
- [x] Delta mapping resolution (`base` / `add` / `remove`)
- [x] Inline comma-group flattening in atom lists
- [x] Active-site carry-forward on higher mapping degrees
- [x] Cross-bond `comment` → `params` default (`; EPO_C1/C2_IPD_N1/N2`)

### Engine wiring

- [x] `BackmapperSettings2` accepts `etree.Element` root (`isinstance` check, not `hasattr(..., "find")`)
- [x] `build_hybrid_gromacs(Path | Settings)` in `network/api.py`
- [x] CLI `build-hybrid`: native path when `prep.bakery_xml` absent + `resolve_data_dir()`

### Rim135 assets

- [x] `examples/epoxy/settings.v2.yaml` — full rim135 (EPO, HDD, IPD, JEF)
- [x] `examples/epoxy/cross_angles_rim135.yaml` — 51 CROSS_BCK triples
- [x] `prep.data_dir` → `bakery/tests/rim135` (assets not duplicated)

### Parity (Tier A)

- [x] `test_rim135_build_hybrid_v2_parity` — byte-identical `ref_hyb_conf.gro`, topology sections match
- [x] Legacy bridge still works: `examples/epoxy/settings.yaml` + `settings.xml`
- [x] `examples/epoxy/run_test.sh` switched to v2 entry point
- [x] `examples/epoxy/README.md` updated

### Bug fixes during v2 bring-up

| Issue | Fix |
|-------|-----|
| IPD K1 atom order (`H9` after `H13` vs ref `N1 H8 H9 N2 H12 H13`) | Reorder degree-0 list in `settings.v2.yaml` |
| YAML block lists as comma strings | `_flatten_atoms()` in v2_loader |
| String paths mistaken for XML Element | `isinstance(..., etree.Element)` in BackmapperSettings2 |
| Missing active_site on IPD degrees 3–4 | Carry-forward in `_emit_cg_bead` |

---

## Remaining

### P1 — Unified `build` (LAMMPS inputs from v2)

- [x] Wire `backmap-prep build` for `prep.engine: network` + v2 YAML
- [x] Emit LAMMPS `.data` + `in.*` from hybrid outputs (`network/lammps_builder.py`)
- [x] Integration test: `test_rim135_build_v2_lammps_smoke`

### P2 — Tier B (LAMMPS MD smoke)

- [ ] Deploy LAMMPS repo on the validation VM (`/home/<vm_user>/sc`)
- [ ] Short backmapping run on rim135 hybrid; confirm no immediate blowup
- [x] Document in research notebook (LAMMPS-only; no ESPResSo++) — `research/notebook/2026-06-22_rim135-tier-b-lammps-smoke.md`

### P3 — Polish & migration tooling

- [x] `backmap-prep migrate-settings settings.xml -o settings.yaml` (CLI stub; exit 2 — full converter pending)
- [ ] Mark v1 XML bridge deprecated in docs (keep read path)
- [x] Pydantic `ConfigDict` migration (schema deprecation warning)
- [ ] Charge transfer rules in v2 (currently commented in ref XML; not required for Tier A)

### P4 — Follow-on systems

- [ ] MF polymerized network
- [ ] PET / hyperbranched AB2/Abx

---

## Verification commands

```bash
cd lammps-user-backmapping

# v2 Tier A parity (integration)
uv run pytest python/tests/test_network.py::test_rim135_build_hybrid_v2_parity -m integration -q

# All rim135 network tests
uv run pytest python/tests/test_network.py -m integration -q

# Shell Tier A (from rim135 data dir)
./examples/epoxy/run_test.sh

# Manual build
uv run backmap-prep build-hybrid examples/epoxy/settings.v2.yaml
```

Requires `bakery/tests/rim135` at `../../../../bakery/tests/rim135` or `BACKMAP_RIM135`.

---

## Key paths

| Path | Role |
|------|------|
| `examples/epoxy/settings.v2.yaml` | Native rim135 config |
| `examples/epoxy/cross_angles_rim135.yaml` | External cross angles |
| `python/src/backmap_prep/network/v2_loader.py` | YAML → XML adapter |
| `python/tests/test_network.py` | Tier A integration tests |
| `bakery/tests/rim135/ref_hyb_*` | Reference outputs |

---

## Architecture (current)

```text
settings.v2.yaml
    → load_settings() + resolve_data_dir()
    → settings_to_xml_root()
    → BackmapperSettings2(Element)
    → prepare_hybrid (bakery structures.py)
    → hyb_conf.gro, hyb_topol.top
```

Legacy path: `settings.yaml` → `prep.bakery_xml` → same engine from file.
