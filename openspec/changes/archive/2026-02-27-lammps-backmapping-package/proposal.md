## Why

The current AdResS backmapping workflow relies on ESPResSo++ (espressopp-bakery) as the MD engine, driven by the AdResSLab Python framework. ESPResSo++ is no longer actively maintained, has a Python 2 dependency, and requires a custom fork with `DynamicResolution`, `VelocityVerletHybrid`, and `FixedVSList` components. LAMMPS is a better long-term platform: actively maintained, widely available, and supports user packages that can be distributed independently.

The existing `lammps-user-adress` package targets full spatial AdResS (regions, thermodynamic force, bidirectional constraints) which is far broader than what backmapping requires. A focused LAMMPS user package dedicated to backmapping — where the entire simulation box transitions from CG to atomistic resolution over time — is both achievable and sufficient for the intended use case.

## What Changes

- **New LAMMPS user package** (`USER-BACKMAP`): A minimal C++ package implementing time-dependent resolution switching for CG→AT backmapping. Contains:
  - `fix backmap/resolution`: manages time-dependent lambda ramp, CG-AT molecule mapping, COM position tracking, and CG force distribution to atomistic particles
  - `pair_style hybrid/backmap`: hybrid pair style that interpolates between AT and CG pair forces weighted by per-atom lambda (or adapt existing `pair_adress` with minimal changes)

- **Python input generator** (`backmap-prep`): A standalone Python tool that takes GROMACS topology (`.top`) and coordinates (`.gro`) and produces LAMMPS input files (`.data` + `.in`) ready for backmapping simulation. Replaces the ESPResSo++ `start_simulation` script and relevant parts of AdResSLab.

- **BREAKING**: This is a new standalone package, not modifications to `lammps-user-adress`. The existing `lammps-user-adress` spatial AdResS code is not affected.

## Capabilities

### New Capabilities

- `fix-backmap-resolution`: LAMMPS fix that implements time-dependent resolution switching — ramping per-atom lambda from 0 (CG) to 1 (AT) over the course of a simulation. Handles CG virtual site position updates (COM of AT atoms), CG force distribution to AT atoms (mass-weighted), and prevention of CG particle integration. Replaces ESPResSo++'s `DynamicResolution` + `VelocityVerletHybrid` + `FixedVSList`.

- `pair-backmap-hybrid`: LAMMPS pair style that computes `F = λ·F_AT + (1−λ)·F_CG` using per-atom lambda values. CG forces are computed only between CG particles. AT forces are computed between AT particles. Supports standard LAMMPS sub-styles (lj/cut, table, etc.) for both AT and CG interactions.

- `backmap-input-generator`: Python tool to convert GROMACS topology + coordinates into LAMMPS data and input files for backmapping. Handles AT+CG particle generation, molecule ID assignment, unit conversion (GROMACS nm/kJ·mol⁻¹ → LAMMPS Å/kcal·mol⁻¹), and tabulated potential format conversion (`.xvg` → LAMMPS `.table`).

### Modified Capabilities

_(none — this is a new standalone package)_

## Impact

- **New code**: ~250 lines C++ (fix), ~200 lines C++ (pair style or adapter), ~300 lines Python (input generator)
- **Dependencies**: LAMMPS stable release (2024+), Python 3.10+ with numpy for the generator
- **Existing code**: `lammps-user-adress` is unaffected. AdResSLab topology parsing (`gromacs_topology.py`, `files_io.py`) can be reused in the Python generator
- **Testing**: Validated against ESPResSo++ backmapping results for the water2 system (3456 SPC water molecules, WCG coarse-grained sites)
