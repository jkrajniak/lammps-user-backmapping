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
cd python
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

```bash
cd examples/dodecane
uv run backmap-prep settings.yaml
```

This produces:

- `dodecane.data` -- LAMMPS data file with both CG and AT atoms
- `in.dodecane` -- LAMMPS input script with all interaction definitions
- `table_b1.table` -- tabulated CG-CG bond potential

### 2. Run the simulation

```bash
lmp -in in.dodecane
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

## Next Steps

- Read the [Theory](theory.md) section to understand how the method works
- See the [Settings Reference](settings-reference.md) for all configuration
  options
- Follow the [Tutorial](tutorial-new-system.md) to set up backmapping for
  your own molecular system
