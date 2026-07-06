# Settings format v2 — design

Status: **implemented for rim135** (2026-06-22). Native loader: `network/v2_loader.py`.
Legacy v1 XML still supported via `prep.bakery_xml`.

## Motivation

Bakery `settings.xml` (v1) dates from ~2015. It works but is hard to review, diff, and validate:

- Duplicate `source_coordinate` / `source_topology` blocks
- Space-separated atom lists in text nodes
- Verbose `1:MOL:ATOM` references when molecule context is known
- Full atom lists repeated for every bonding degree (e.g. IPD `K1`: five near-identical blocks)
- Split `hybrid_configuration` + `hybrid_topology` with free-text cross-angle blobs
- Cryptic charge-transfer strings (`EPO:C1#H25,EPO:C21#H26,…`)

Roughly 250 of 390 lines in `bakery/tests/rim135/settings.xml` repeat information.

**v2** is YAML-native, aligned with linear `backmap-prep` (`schema.py`), and maps to the same internal structures as v1.

## Goals

- One configuration shape for linear and network prep
- Human-readable, diff-friendly, Pydantic-validated
- Delta-based degree mapping (inherit + add/remove atoms)
- Structured charge rules and cross interactions
- Legacy v1 XML supported via adapter until examples are migrated

## Non-goals (v2.0)

- Changing hybrid generation semantics
- Live reaction during λ ramp
- Replacing GROMACS `.itp` / `.gro` assets

## Top-level layout

```yaml
version: 2
prep:
  engine: network          # linear | network

molecules: [...]           # CG→AT mapping per species
cg_system: [...]           # cured CG melt (was cg_configuration)
hybrid: [...]              # output paths + includes (was hybrid_configuration + hybrid_topology shell)
cross_interactions: [...]  # cross-CG bonds/angles (shared with linear prep)
simulation: [...]          # LAMMPS backmapping parameters
output: [...]
```

## Section reference

### `molecules[]`

| Field | Required | Notes |
|-------|----------|-------|
| `name` | yes | CG molecule name in cured topology |
| `ident` | no | AT topology name; defaults to `name` |
| `source` | yes | List of degree-dependent AT files (see below) |
| `beads` | yes | CG bead definitions |
| `charge_management` | no | Equilibration + transfer rules |

#### Unified `source` (replaces `source_coordinate` + `source_topology`)

```yaml
source:
  - molecule_degree: 0
    coordinates: epon-828.gro
    topology: epon-828.itp
  - molecule_degree: 1
    when: A1
    coordinates: epon-828_deg1_A1.gro
    topology: epon-828_deg1_A1.itp
```

Each row is one fragment variant. Coordinate and topology paths cannot drift apart.

For simple molecules (IPD monomer):

```yaml
source:
  coordinates: ipd.gro
  topology: ipd.itp
```

#### Atom references

Inside a molecule, atoms are short names (`C1`, `H25`). Full bakery form `1:EPO:C1` remains valid for edge cases; `chain` defaults to `1`.

#### Bead `mapping` with deltas

```yaml
beads:
  - name: K1
    type: K
    mapping:
      - degree: 0
        atoms: [C1, C2, C3, ..., N1, H8, N2, H12, H13]

      - degree: 1
        base: 0
        remove: [H8]
        active_sites:
          - { atom: N1, max_degree: 3 }

      - degree: 2
        base: 0
        remove: [H8, H9, H12]
        active_sites:
          - { atom: N1, max_degree: 3 }
          - { atom: N2, max_degree: 3 }
```

Resolution order: start from `base` degree atom set → apply `remove` → apply `add` → attach `active_sites`.

Equivalent v1: multiple `<beads degree="…">` blocks with full atom lists.

#### Charge management

```yaml
charge_management:
  equilibrate: true
  transfers:
    - when: { atom: N1, degree: 2 }
      from: H8
      to:
        - { partner: EPO, atom: C1,  add_atom: H25 }
        - { partner: EPO, atom: C21, add_atom: H26 }
        - { partner: HDD, atom: C1,  add_atom: H24 }
```

v1 equivalent: `<charge_transfer when="IPD:N1:2" from_atom="IPD:H8" when_to_atom="…" />`

### `cg_system`

```yaml
cg_system:
  coordinates: cg_cl_conf.gro
  topology: cg_cl_topol.top
  predefined_active_sites: active_sites.dat   # optional
```

### `hybrid`

```yaml
hybrid:
  coordinates: hyb_conf.gro
  topology: hyb_topol.top
  includes:
    - oplsaa.ff/forcefield.itp
  molecule_type:
    name: HYB
    exclusion: 3
```

Cross bonded terms live in `cross_interactions`, not as unstructured text under `hybrid_topology`.

### `cross_interactions`

Same schema as linear prep (`schema.py`). Example:

```yaml
cross_interactions:
  bonds:
    - comment: EPO epoxide — IPD amine
      pairs:
        - [EPO:C1, IPD:N1]
        - [EPO:C21, IPD:N2]
  angles:
  # Long predefined lists may move to external YAML:
  # angles_file: cross_angles_rim135.yaml
```

## Optional XML serialization (v2)

If XML is required, v2 XML is **isomorphic to YAML** — not a third dialect:

```xml
<settings version="2">
  <molecule name="EPO" ident="EPO">
    <source molecule_degree="0" coordinates="epon-828.gro" topology="epon-828.itp"/>
    <bead name="A1" type="A">
      <mapping degree="1" molecule_degree="0">
        <atom>C1</atom><atom>O1</atom>
      </mapping>
      <mapping degree="2" base="1" molecule_degree="1,2">
        <add>H25</add>
        <active_site atom="C1" max_degree="4"/>
      </mapping>
    </bead>
  </molecule>
</settings>
```

**Authoring default: YAML.** XML is an export format only.

## Migration path

```
bakery settings.xml v1 ──► legacy adapter ──┐
                                              ├──► internal model ──► prepare_hybrid()
settings.yaml v2 ────────► native loader ────┘
```

1. **Now:** v1 XML → `BackmapperSettings2._parse()` (ported bakery)
2. **Next:** v2 YAML → `Settings` → populate same internal structures
3. **Tooling:** `backmap-prep migrate-settings settings.xml -o settings.yaml`
4. **Deprecation:** v1 read-only; new examples in v2 only

## Rim135 size estimate

| Section | v1 lines | v2 est. |
|---------|----------|---------|
| EPO + HDD + JEF mapping | ~290 | ~60 |
| IPD K1 degree variants | ~75 | ~15 |
| CG + hybrid I/O | ~100 | ~25 |
| Cross terms | ~50 | ~50 or external file |
| **Total** | **~390** | **~150** |

See `examples/epoxy/settings.v2.yaml` — wired to `backmap-prep build-hybrid`; Tier A parity green.

## Relation to existing code

| Component | v1 | v2 target |
|-----------|----|-----------|
| `python/src/backmap_prep/schema.py` | v2 models + delta resolution + `angles_file` merge | Docs + `ConfigDict` cleanup |
| `network/v2_loader.py` | YAML → XML ElementTree | Optional direct internal model (skip XML) |
| `examples/epoxy/settings.v2.yaml` | Primary rim135 config | Wire into `build` (LAMMPS) |
| `examples/epoxy/settings.yaml` | Legacy bridge to `settings.xml` | Deprecate when `build` uses v2 |

## Open questions

- External file for long cross-angle lists (rim135 has ~50 triples)?
- `base` + `add`/`remove` vs explicit `atoms` only — allow both?
- Version field: `version: 2` in YAML vs `prep.engine: network` alone?
