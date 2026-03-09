# Polyethylene (10:1 AA) — Backmapping

Backmapping of polyethylene from CG to all-atom resolution using OPLS/AA
force field with ~30 atoms per CG bead.

## System

- **10 chains** (reduced from original 75), 3120 atoms total (100 CG + 3020 AT)
- CG model: 10 beads per chain (A1, B2–B9, A10), tabulated IBI interactions
- AT model: OPLS/AA all-atom (C: opls_135/136, H: opls_140)
- Temperature: 423 K

Per-chain mapping (30 atoms per bead):

| CG bead | Type | Heavy atoms | AT atoms per bead |
|---------|------|-------------|-------------------|
| A1      | A    | C1–C9+C99   | 31 (10C + 21H)    |
| B2      | B    | C10–C19     | 30 (10C + 20H)    |
| B3      | B    | C20–C29     | 30 (10C + 20H)    |
| ...     | B    | ...         | 30                |
| B9      | B    | C80–C89     | 30 (10C + 20H)    |
| A10     | A    | C90–C98+C100| 31 (10C + 21H)    |

## Quick Start

```bash
# 1. Generate LAMMPS files from settings
cd python && uv run backmap-prep ../examples/pe_10/settings.yaml

# 2. Run LAMMPS
cd ../examples/pe_10 && lmp -in in.pe_10
```

## Input Files

| File | Description |
|------|-------------|
| `settings.yaml` | backmap-prep configuration |
| `pe_single.gro` | Single-chain AT reference (302 atoms) |
| `topol_aa.top` | AT topology (OPLS/AA with explicit H) |
| `cg_conf.gro` | CG coordinates (10 chains, 100 beads) |
| `topol_cg.top` | CG topology (tabulated bonds) |
| `table_b1.xvg` | CG B–B bond table |
| `table_b2.xvg` | CG A–B bond table |
| `table_A_A.xvg` | CG A–A pair table |
| `table_A_B.xvg` | CG A–B pair table |
| `table_B_B.xvg` | CG B–B pair table |

## Origin

Migrated from [bakery](https://github.com/bakery-cg2at/bakery) `examples/pe_10/`.
Original system: 75 chains at 423 K; reduced to 10 chains for quick testing.

## Large-scale variant

A 75-chain variant is in `large/`. Copy `conf_cg.gro` from `bakery/examples/pe_10/` as `large/cg_conf.gro` if needed, then run `backmap-prep` with `large/settings.yaml`. See `large/README.md` and [Large-scale examples](../../docs/large-scale-examples.md).
