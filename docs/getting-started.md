# Getting Started

This guide walks you through installing the LAMMPS Backmapping Package and
running the included dodecane example.

## Prerequisites

- **LAMMPS** source tree (C++17 compatible)
- **CMake** &ge; 3.16
- **Python** &ge; 3.10
- **uv** (Python dependency manager)

## Installation

### 1. Install the C++ package into LAMMPS

Clone the backmapping repository and use the install script to copy the source
files into your LAMMPS tree:

```bash
git clone https://github.com/jkrajniak/lammps-user-backmapping.git
cd lammps-user-backmapping

./install.sh /path/to/lammps
```

Then rebuild LAMMPS with the package enabled:

=== "CMake (recommended)"

    ```bash
    cd /path/to/lammps/build
    cmake -D PKG_BACKMAP=yes ../cmake
    cmake --build . -j$(nproc)
    ```

=== "Traditional make"

    ```bash
    cd /path/to/lammps/src
    make yes-backmap
    make mpi -j$(nproc)
    ```

To uninstall later:

```bash
./install.sh --uninstall /path/to/lammps
```

### 2. Install the Python CLI

```bash
uv sync
```

This installs the `backmap-prep` command-line tool into the virtual environment.
You can also install with development dependencies:

```bash
make install-dev
```

Verify the installation:

```bash
uv run backmap-prep --help
```

## Running the Dodecane Example

The repository includes a complete working example with a dodecane system
(6 CG beads mapped to 12 united-atom carbons).

### 1. Generate LAMMPS input files

From the repository root:

```bash
uv run backmap-prep examples/dodecane/settings.yaml
```

Output is written to `examples/dodecane/`. This produces:

- `dodecane.data` -- LAMMPS data file with both CG and AT atoms
- `in.dodecane` -- LAMMPS input script with all interaction definitions
- `table_*.table` -- tabulated CG bond and pair potentials

### 2. Run the simulation

```bash
cd examples/dodecane && lmp -in in.dodecane
```

The simulation runs three phases:

1. **CG equilibration** (10,000 steps) -- lambda ramp is inactive, system
   equilibrates at CG level
2. **Backmapping** (10,000 steps) -- lambda ramps from 0 to 1, gradually
   introducing AT interactions
3. **AT production** (10,000 steps) -- fully atomistic simulation

### 3. Check the output

- `dump.backmap` -- trajectory with per-atom lambda values (column `f_bm`)
- `log.lammps` -- thermodynamic output showing energy convergence

!!! tip
    The per-atom lambda value is accessible via `f_bm` in dump commands and
    can be used for visualization or post-processing.

## More Examples

Beyond dodecane, the repository includes several other systems of increasing
complexity:

| Example | CG beads | AT atoms/bead | Force field | Description |
|---------|----------|---------------|-------------|-------------|
| `examples/dodecane/` | 6 | 2 | GROMOS UA | Linear alkane, simplest case |
| `examples/pe/` | 50 | 2 | OPLS UA | Polyethylene (100 C), 2:1 mapping |
| `examples/pe4/` | 25 | 4 | OPLS UA | Polyethylene (100 C), 4:1 mapping |
| `examples/pe_10/` | 10 | ~30 | OPLS AA | Polyethylene (100 C), 10:1 mapping with H |
| `examples/pe_aa/` | 50 | 6-7 | OPLS AA | Polyethylene (100 C), 2:1 mapping with H |
| `examples/melamine/` | 3 | 9 | OPLS AA | Melamine-formaldehyde, triangular CG |

Each example directory contains a `README.md` with system-specific
instructions.

### Large-scale variants

Each example has a `large/` subdirectory with (or instructions for) production-scale
inputs (e.g. 75 chains, 500 molecules) sourced from the [bakery](https://github.com/bakery-cg2at/bakery) project. See [Large-scale examples](large-scale-examples.md) for how to obtain inputs and run them.

## Next Steps

- Read the [Theory](theory.md) section to understand how the method works
- See the [Settings Reference](settings-reference.md) for all configuration
  options
- Follow the [Tutorial](tutorial-new-system.md) to set up backmapping for
  your own molecular system
