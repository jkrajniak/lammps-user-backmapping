# Melamine-Formaldehyde — Backmapping

Backmapping of melamine-formaldehyde (MF) molecules from CG to atomistic
resolution. Triangular CG topology with 3 beads per molecule.

## System

- **50 molecules** (reduced from original 500), 1500 atoms total (150 CG + 1350 AT)
- CG model: 3 beads per molecule (A1, A2, A3), triangular bonding, tabulated bonds
- AT model: OPLS/AA (triazine ring + methylol groups)
- Temperature: 300 K

Per-molecule mapping:

| CG bead | Type | AT atoms (9 per bead) |
|---------|------|-----------------------|
| A1      | A    | N11, C11, N21, C1, O1, H01, H011, H012, H21 |
| A2      | A    | N12, C12, N22, C2, O2, H02, H021, H022, H31 |
| A3      | A    | N13, C13, N23, C3, O3, H03, H031, H032, H11 |

CG topology: A1–A2–A3–A1 (triangle, 3 bonds per molecule).

## Quick Start

```bash
# 1. Generate LAMMPS files from settings
uv run backmap-prep examples/melamine/settings.yaml

# 2. Run LAMMPS
cd examples/melamine && lmp -in in.melamine
```

## Input Files

| File | Description |
|------|-------------|
| `settings.yaml` | backmap-prep configuration |
| `single_mf.gro` | Single-molecule AT reference (27 atoms) |
| `topol_aa.top` | AT topology (OPLS/AA, 8 atom types) |
| `cg_conf.gro` | CG coordinates (50 molecules, 150 beads) |
| `topol_cg.top` | CG topology (tabulated bonds, triangle) |
| `table_b1.xvg` | CG bond table |
| `table_A_A.xvg` | CG A–A pair table |

## Generated Files

| File | Description |
|------|-------------|
| `melamine.data` | LAMMPS data file |
| `in.melamine` | LAMMPS input script |
| `table_b1.table` | Converted CG bond table |

## Origin

Migrated from [bakery](https://github.com/bakery-cg2at/bakery) `examples/melamine/`.
Original system: 500 molecules; reduced to 50 for quick testing.

## Large-scale variant

A 500-molecule variant is in `large/`. Copy `cg_conf_500.gro` from `bakery/examples/melamine/` as `large/cg_conf.gro` if needed, then run `backmap-prep` with `large/settings.yaml`. See `large/README.md` and [Large-scale examples](../../docs/large-scale-examples.md).
