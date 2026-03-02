# Polyethylene (4:1 UA) — Backmapping

Backmapping of polyethylene from CG to atomistic resolution using OPLS
united-atom force field with 4 UA atoms per CG bead.

## System

- **10 chains** (reduced from original 75), 1250 atoms total (250 CG + 1000 AT)
- CG model: 25 beads per chain (A1, B1–B23, A2), tabulated IBI interactions
- AT model: OPLS united-atom (CH₃ / CH₂)
- Temperature: 423 K

Per-chain mapping (4 atoms per bead):

| CG bead | Type | AT atoms       | AT types          |
|---------|------|----------------|-------------------|
| A1      | A    | C1–C4          | CH₃, CH₂, CH₂, CH₂ |
| B1      | B    | C5–C8          | CH₂, CH₂, CH₂, CH₂ |
| B2      | B    | C9–C12         | CH₂, CH₂, CH₂, CH₂ |
| ...     | B    | ...            | CH₂, CH₂, CH₂, CH₂ |
| B23     | B    | C93–C96        | CH₂, CH₂, CH₂, CH₂ |
| A2      | A    | C97–C100       | CH₂, CH₂, CH₂, CH₃ |

## Quick Start

```bash
# 1. Generate LAMMPS files from settings
cd python && uv run backmap-prep ../examples/pe4/settings.yaml

# 2. Run LAMMPS
cd ../examples/pe4 && lmp -in in.pe4
```

## Input Files

| File | Description |
|------|-------------|
| `settings.yaml` | backmap-prep configuration |
| `pe_single.gro` | Single-chain AT reference (100 atoms) |
| `topol_aa.top` | AT topology (OPLS/UA bonds, angles, dihedrals) |
| `cg_conf.gro` | CG coordinates (10 chains, 250 beads) |
| `topol_cg.top` | CG topology (tabulated bonds, angles, dihedrals) |
| `table_b1.xvg` | CG bond table |
| `table_a1.xvg` | CG angle table |
| `table_d1.xvg` | CG dihedral table |
| `table_A_A.xvg` | CG A–A pair table |
| `table_A_B.xvg` | CG A–B pair table |
| `table_B_B.xvg` | CG B–B pair table |

## Origin

Migrated from [bakery](https://github.com/bakery-cg2at/bakery) `examples/pe4/`.
Original system: 75 chains at 423 K; reduced to 10 chains for quick testing.
