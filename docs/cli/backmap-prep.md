# backmap-prep CLI

`backmap-prep` reads a YAML settings file describing a CG-to-AT backmapping
setup and generates all files needed to run the simulation in LAMMPS.

## Installation

```bash
cd python
uv sync
```

The `backmap-prep` command is then available in the virtual environment:

```bash
uv run backmap-prep --help
```

## Syntax

```
backmap-prep SETTINGS [--settings SETTINGS] [--output-prefix PREFIX]
```

## Arguments

### `SETTINGS` (positional, required)

Path to the YAML settings file. See the
[Settings Reference](../settings-reference.md) for the full file format.

```bash
uv run backmap-prep settings.yaml
```

### `--settings SETTINGS`

Alternative way to specify the settings file (for backward compatibility).
If both the positional argument and `--settings` are provided, `--settings`
takes precedence.

```bash
uv run backmap-prep --settings settings.yaml
```

### `--output-prefix PREFIX`

Override the `output.prefix` value from the YAML settings file. This changes
the names of the generated files.

```bash
uv run backmap-prep settings.yaml --output-prefix myrun
# Produces: myrun.data, in.myrun
```

## Generated Files

`backmap-prep` produces three types of files in the same directory as the
settings file:

| File | Description |
|------|-------------|
| `<prefix>.data` | LAMMPS data file with CG + AT atoms, bonds, angles |
| `in.<prefix>` | LAMMPS input script with complete simulation setup |
| `*.table` | Tabulated potential files (converted from GROMACS XVG) |

### Data File

Contains:

- All CG beads and AT atoms with correct types and masses
- Complete bond, angle, and dihedral topology
- Simulation box dimensions from the CG system

### Input Script

Contains:

- `pair_style backmap` with AT and CG sub-styles
- All `pair_coeff` entries (atomistic, cg, none)
- `bond_style hybrid` with standard and lambda-weighted bond styles
- `angle_style` definitions
- `fix backmap` with parameters from the settings file
- Thermostat and integrator setup
- Three-phase simulation protocol (equilibration, backmapping, production)
- Dump and thermo output commands

### Table Files

CG tabulated potentials referenced in `cross_interactions.bonds[].table`
are converted from GROMACS XVG format to LAMMPS table format.

## Examples

Basic usage:

```bash
uv run backmap-prep settings.yaml
```

With custom output prefix:

```bash
uv run backmap-prep settings.yaml --output-prefix dodecane_run1
```

Running from a different directory:

```bash
uv run backmap-prep examples/dodecane/settings.yaml
```

## Exit Codes

| Code | Meaning |
|------|---------|
| 0 | Success |
| 1 | Settings file not found |
| (exception) | Invalid settings file or processing error |

## Related

- [Settings Reference](../settings-reference.md) -- complete YAML file
  documentation
- [Getting Started](../getting-started.md) -- installation and first run
- [Tutorial](../tutorial-new-system.md) -- end-to-end setup guide
