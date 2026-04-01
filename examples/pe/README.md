# Polyethylene (2:1 UA) — Backmapping

Backmapping of polyethylene from CG to atomistic resolution using OPLS
united-atom force field with 2 UA atoms per CG bead.

## System

- **10 chains** (reduced from original 75), 1500 atoms total (500 CG + 1000 AT)
- CG model: 50 beads per chain (A1, B1–B48, A2), tabulated IBI interactions
- AT model: OPLS united-atom (CH₃ / CH₂)
- Temperature: 423 K

Per-chain mapping:

| CG bead | Type | AT atoms   | AT types  |
|---------|------|------------|-----------|
| A1      | A    | C1, C2     | CH₃, CH₂ |
| B1      | B    | C3, C4     | CH₂, CH₂ |
| B2      | B    | C5, C6     | CH₂, CH₂ |
| ...     | B    | ...        | CH₂, CH₂ |
| B48     | B    | C97, C98   | CH₂, CH₂ |
| A2      | A    | C99, C100  | CH₂, CH₃ |

## Quick Start

```bash
# 1. Generate LAMMPS files from settings
uv run backmap-prep examples/pe/settings.yaml

# 2. Run LAMMPS
cd examples/pe && lmp -in in.pe
```

## Input Files

| File | Description |
|------|-------------|
| `settings.yaml` | backmap-prep configuration |
| `pe_single.gro` | Single-chain AT reference (100 atoms) |
| `topol_aa.top` | AT topology (OPLS/UA bonds, angles, dihedrals) |
| `cg_conf.gro` | CG coordinates (10 chains, 500 beads) |
| `topol_cg.top` | CG topology (tabulated bonds, angles, dihedrals) |
| `table_b1.xvg` | CG B–B bond table |
| `table_b2.xvg` | CG A–B bond table |
| `table_a1.xvg` | CG A–B–B / B–B–A angle table |
| `table_a2.xvg` | CG B–B–B angle table |
| `table_d1.xvg` | CG dihedral table |
| `table_A_A.xvg` | CG A–A pair table |
| `table_A_B.xvg` | CG A–B pair table |
| `table_B_B.xvg` | CG B–B pair table |

## Generated Files

| File | Description |
|------|-------------|
| `pe.data` | LAMMPS data file |
| `in.pe` | LAMMPS input script |
| `table_b1.table` | Converted CG B–B bond table |
| `table_b2.table` | Converted CG A–B bond table |

## Origin

Migrated from [bakery](https://github.com/bakery-cg2at/bakery) `examples/pe/`.
Original system: 75 chains at 423 K; reduced to 10 chains for quick testing.

## Large-scale variant

A 75-chain variant is in `large/`. Copy `cg_conf.gro` from `bakery/examples/pe/` into `large/` if needed, then run `backmap-prep` with `large/settings.yaml`. See `large/README.md` and [Large-scale examples](../../docs/large-scale-examples.md).
