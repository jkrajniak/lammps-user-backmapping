# Polyethylene (2:1 AA) — Backmapping

Backmapping of polyethylene from CG to all-atom resolution using OPLS/AA
force field with explicit hydrogens. Same CG model as the PE (2:1 UA)
example but with an all-atom target.

## System

- **10 chains** (reduced from original 75), 3520 atoms total (500 CG + 3020 AT)
- CG model: 50 beads per chain (A1, B2–B49, A50), tabulated IBI interactions
- AT model: OPLS/AA all-atom (C: opls_135/136, H: opls_140)
- Temperature: 423 K

Per-chain mapping:

| CG bead | Type | AT atoms per bead  |
|---------|------|--------------------|
| A1      | A    | 7 (2C + 5H)        |
| B2–B49  | B    | 6 (2C + 4H)        |
| A50     | A    | 7 (2C + 5H)        |

## Quick Start

```bash
# 1. Generate LAMMPS files from settings
uv run backmap-prep examples/pe_aa/settings.yaml

# 2. Run LAMMPS
cd examples/pe_aa && lmp -in in.pe_aa
```

## Input Files

| File | Description |
|------|-------------|
| `settings.yaml` | backmap-prep configuration |
| `pe_single.gro` | Single-chain AT reference (302 atoms) |
| `topol_aa.top` | AT topology (OPLS/AA with explicit H) |
| `cg_conf.gro` | CG coordinates (10 chains, 500 beads) |
| `topol_cg.top` | CG topology (tabulated bonds, angles, dihedrals) |
| `table_b1.xvg` | CG B–B bond table (shared with PE UA) |
| `table_b2.xvg` | CG A–B bond table |
| `table_a1.xvg` | CG A–B–B angle table |
| `table_a2.xvg` | CG B–B–B angle table |
| `table_d1.xvg` | CG dihedral table |
| `table_A_A.xvg` | CG A–A pair table |
| `table_A_B.xvg` | CG A–B pair table |
| `table_B_B.xvg` | CG B–B pair table |

## Origin

Migrated from [bakery](https://github.com/bakery-cg2at/bakery) `examples/pe_aa/prepare/`.
Uses the same CG tables as the PE (2:1 UA) example — the CG representation is identical.
Original system: 75 chains at 423 K; reduced to 10 chains for quick testing.

## Large-scale variant

A larger-chain variant is in `large/`. Copy `conf_cg.gro` from `bakery/examples/pe_aa/prepare/` as `large/cg_conf.gro` if needed, then run `backmap-prep` with `large/settings.yaml`. See `large/README.md` and [Large-scale examples](../../docs/large-scale-examples.md).
