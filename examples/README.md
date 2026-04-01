# Backmapping Examples

Each example demonstrates CG-to-atomistic backmapping for a different
molecular system using `backmap-prep` and LAMMPS.

## Available Examples

| Directory | System | CG beads | AT atoms/bead | Force field | Mapping |
|-----------|--------|----------|---------------|-------------|---------|
| `dodecane/` | Dodecane (C12) | 6 | 2 | GROMOS UA | Linear chain |
| `pe/` | Polyethylene (C100) | 50 | 2 | OPLS UA | Linear chain, 2:1 |
| `pe4/` | Polyethylene (C100) | 25 | 4 | OPLS UA | Linear chain, 4:1 |
| `pe_10/` | Polyethylene (C100) | 10 | ~30 | OPLS AA | Linear chain, 10:1 |
| `pe_aa/` | Polyethylene (C100) | 50 | 6-7 | OPLS AA | Linear chain, 2:1 with H |
| `melamine/` | Melamine-formaldehyde | 3 | 9 | OPLS AA | Triangular CG topology |

## Quick Start

All examples follow the same workflow:

```bash
# 1. Generate LAMMPS input files
uv run backmap-prep examples/<name>/settings.yaml

# 2. Run the backmapping simulation
cd examples/<name> && lmp -in in.<name>
```

## Choosing an Example

- **Getting started**: Start with `dodecane/` — smallest and simplest system
- **Linear polymers (UA)**: `pe/` and `pe4/` show different CG-to-AT mapping ratios
- **All-atom with H**: `pe_10/` and `pe_aa/` demonstrate OPLS/AA backmapping
- **Non-linear topology**: `melamine/` shows a triangular CG bonding pattern

## Large-scale variants

Each example has a `large/` subdirectory for production-scale systems (e.g. 75 chains for PE, 500 molecules for melamine). Inputs are sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project; see `large/README.md` in each directory and the [Large-scale examples](https://jkrajniak.github.io/lammps-user-backmapping/large-scale-examples/) doc for how to obtain files and run `backmap-prep`. To validate file generation without running LAMMPS:

```bash
./scripts/validate-large-scale-prep.sh dodecane
```

## Origin

The PE and melamine examples were migrated from the
[bakery](https://github.com/bakery-cg2at/bakery) project. System sizes
have been reduced for quick testing (10 chains / 50 molecules vs. originals
with 75 chains / 500 molecules).
