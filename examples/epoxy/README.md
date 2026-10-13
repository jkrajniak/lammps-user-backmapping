# Epoxy (RIM135) network backmapping

Reactive epoxy-amine network: EPO + IPD + HDD (+ JEF) CG melt backmapped to a
GROMACS hybrid coordinate/topology pair for subsequent LAMMPS backmapping.

This example uses the **Phase 3 network engine** (ported bakery `structures.py`).

Two equivalent entry points:

- **Settings v2 (recommended):** `settings.v2.yaml` — native YAML with delta mapping
- **Legacy bridge:** `settings.yaml` — points at bakery `settings.xml` via `prep.bakery_xml`

## Prerequisites

- `uv` and this package installed (`uv sync` from repo root)
- Bakery rim135 test data at `../../../../bakery/tests/rim135` (sibling repo), **or**
  set `BACKMAP_RIM135` to that directory
- OPLS-AA forcefield for func-3 dihedraltypes lookup: bundled at
  `examples/forcefield/oplsaa.ff/` (also copied beside this example as `forcefield/`).
  Set `prep.forcefield_dir: forcefield` in `settings.v2.yaml`, or point at a GROMACS
  install via `GMXDATA` / `GROMACS_DATA`.

## Where the data comes from

Neither entry point is self-contained; both read assets of the earlier publications.

- **`settings.v2.yaml`** reads the CG coordinates, topology and tables from the published JCC 2018
  reproducibility archive ([`paper-reverse-mapping-polymer-networks`](https://github.com/jkrajniak/paper-reverse-mapping-polymer-networks),
  `preparation/rim135_small/backmapping/aa`). `data_dir` and `tables_dir` point at it as a sibling of
  this repository; change them if it is elsewhere.
- **`settings.yaml`** (legacy bridge) reads a bakery `settings.xml`, and the CG angle tables
  `table_a<N>.xvg` it refers to must be in the directory of the settings file or in the tables
  directory. Without them the build stops with `Missing angle table table_a1.xvg for func-8 angle`.
  The archive above holds a `settings.xml` and the tables in the same directory:
  set `prep.bakery_xml` to that `settings.xml` and copy its `*.xvg` files next to `settings.yaml`.

## Quick start

```bash
./run_test.sh
```

Or manually:

```bash
export BACKMAP_RIM135=/path/to/bakery/tests/rim135   # optional
cd "$BACKMAP_RIM135"
uv run --directory ../../.. backmap-prep build-hybrid /path/to/examples/epoxy/settings.v2.yaml
```

## Tier A parity

`run_test.sh` builds the hybrid system with `chain_rng_seed: 42` and checks:

- `hyb_conf.gro` byte-identical to `ref_hyb_conf.gro`
- `hyb_topol.top` section keys match via `backmap-prep compare-topology`
  (cross angles compared in canonical orientation)

Rim135 `settings.xml` includes `predefined_active_sites.txt` so cross-bond
active-site choices match the reference topology.

## Files

| File | Role |
|------|------|
| `settings.v2.yaml` | Native v2 rim135 config (delta mapping, external cross angles, forcefield dir) |
| `forcefield/oplsaa.ff/` | Bundled OPLS-AA for dihedraltypes lookup (symlink to shared `examples/forcefield/`) |
| `cross_angles_rim135.yaml` | Cross-bead angle triples for hybrid topology |
| `settings.yaml` | Legacy bridge: `prep.bakery_xml` → bakery `settings.xml` |
| `run_test.sh` | Runs `build-hybrid` in rim135 data dir and compares to refs |
| `large/in.rim135` | LAMMPS Tier B/C 3-phase production script (see `large/README.md`) |

All GROMACS/CG assets remain in `bakery/tests/rim135/` (not duplicated here).

## Supported build paths

| Path | Command | MD status |
|------|---------|-----------|
| **Bakery reference (recommended)** | `backmap-prep build` from static `cg_cl_conf.gro` | **PASS** — VM 10k+10k |
| Rebuild from equilibrated CG | `backmap-prep rebuild --cg-frame …` | **Experimental** — PBC/image flags differ; may explode in MD |

Use the bakery `build` path for production runs and paper validation until the
rebuild workflow is fixed (see `research/notebook/2026-06-22_rim135-tier-b-lammps-smoke.md`).

## Structural validation

After a successful VM backmap, compare C–O / C–N RDF peaks vs the 2017/2018 JCC
paper reference:

```bash
uv run examples/epoxy/compare_rim135_structure.py \
  --dump examples/epoxy/validation/dump.backmap \
  --gro "$BACKMAP_RIM135/hyb_conf.gro" \
  --final-data examples/epoxy/validation/rim135_final.data \
  --paper-dir ../paper-reverse-mapping-polymer-networks/paper/rim135_small \
  --plot examples/epoxy/validation/rdf_comparison.png
```

Pinned report: `validation/structural_validation_report.txt` (4/4 peak metrics PASS vs AA ref).

## Settings format v2

Rim135 is fully defined in `settings.v2.yaml` (no `bakery_xml` required). The
engine converts v2 YAML to bakery-compatible XML internally via `network/v2_loader.py`.

Design spec: `openspec/changes/phase3-network-backmapping/settings-format-v2.md`

Tier A parity is checked by `test_rim135_build_hybrid_v2_parity` (byte-identical
`hyb_conf.gro` and topology section parity vs `ref_hyb_*`).
