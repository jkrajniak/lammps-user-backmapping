# AGENTS.md

## Project Overview

Time-dependent backmapping from coarse-grained to atomistic resolution,
implemented as a LAMMPS user package. Two main components:

- **C++ LAMMPS styles** (`src/`) — `fix backmap`, `pair_style backmap`,
`bond_style backmap/{harmonic,table}`, `angle_style backmap/{harmonic,table}`.
- **Python CLI** (`python/`) — `backmap-prep` generates LAMMPS data files,
input scripts, and interaction tables from GROMACS topologies.
- **Docs** (`docs/`) — MkDocs Material site; **examples** in `examples/`.

Key papers: Krajniak et al., JCTC 2016 (10.1021/acs.jctc.6b00595) — linear
polymers; Krajniak, Zhang et al., J. Comput. Chem. 2018, *39*, 648--664
(10.1002/jcc.25129; **online Dec. 2017**) — **networks** (epoxy, melamine),
hyperbranched, PET. Project goal: LAMMPS + `backmap-prep` instead of legacy
ESPResSo++ bonded code for broader adoption.

All rules in `openspec/project.md` are mandatory — read it before changing code.

## Repository Layout

```
src/                    # C++ LAMMPS styles
python/src/backmap_prep/# Python package source
python/tests/           # pytest unit tests
pyproject.toml          # Python project metadata, ruff, mypy, pytest config
examples/               # Example simulations (dodecane, polyethylene, melamine)
docs/                   # MkDocs Material documentation site
openspec/               # Specifications and change tracking
```

---

## Build, Lint, and Test Commands

### Python (backmap-prep)

```bash
# Setup
make install              # uv sync
make install-dev          # uv sync --extra dev
make install-hooks        # Install pre-commit hooks

# Linting & formatting
make lint                 # ruff check python/src/ python/tests/
make format               # ruff format (in-place)
make format-check         # ruff format --check (dry run)
make typecheck            # mypy python/src/

# Testing
make test                 # Run all tests
make test-cov             # Tests + coverage report
uv run pytest python/tests/test_foo.py                # Single file
uv run pytest python/tests/test_foo.py::test_bar      # Single test
uv run pytest python/tests/test_foo.py -k "pattern"   # By name pattern

# Pre-commit
make pre-commit           # Staged files only
make pre-commit-all       # All files

# Docs
make docs                 # Build MkDocs site (strict)
make docs-serve           # Serve locally

# Cleanup
make clean                # Remove caches and build artifacts
```

### C++ (LAMMPS styles)

The LAMMPS source tree is at `/Users/jakubkrajniak/Work/Science/lammps` with a
CMake build in `build/`. To compile the backmap package:

```bash
# Copy source files into the LAMMPS tree
cp src/*.cpp src/*.h /Users/jakubkrajniak/Work/Science/lammps/src/

# Rebuild LAMMPS (from the build directory)
cmake --build /Users/jakubkrajniak/Work/Science/lammps/build -j$(nproc)

# The binary is at build/lmp; optionally install:
cmake --install /Users/jakubkrajniak/Work/Science/lammps/build
```

Lint with:

```bash
clang-format --style=file --fallback-style=Google src/*.cpp src/*.h
```

---

## Code Style — Python

- **Python ≥ 3.10**, full type annotations, Pydantic v2 for schemas.
- **Formatter/linter**: ruff (line length 100, config in `pyproject.toml`).
- **Type checker**: mypy strict mode, Pydantic plugin enabled.
- **Dependencies**: `uv` exclusively — never `pip`. Use `uv add`, `uv sync`.
- **Imports**: `from __future__ import annotations`; group stdlib → third-party → local.
- **Naming**: `snake_case` modules/functions, `PascalCase` classes, `UPPER_SNAKE_CASE` constants.
- **Structure**: functions under 50 lines; `__all__` in public modules; Pydantic models or
dataclasses for data containers.
- **Errors**: guard clauses at top, early returns, specific exceptions (never bare `except:`).
- **Style**: f-strings, trailing commas in multi-line calls, no narrating comments.
- **Tests**: every module needs a `python/tests/test_<module>.py`; use parametrised tests
and fixtures; cover edge cases.

### Import Example

```python
from __future__ import annotations

import sys
from pathlib import Path

import yaml
from pydantic import BaseModel

from backmap_prep.parsers import parse_gro
```

## Code Style — C++

- **C++17**, LAMMPS coding style: 2-space indent, `ClassName`, `lower_snake_case`.
- **Formatter**: clang-format with Google style (see `.pre-commit-config.yaml`).
- **Includes**: own header first, then `<cstdlib>` etc., then LAMMPS headers.
- **Members**: trailing underscore (`cut_global_`); init to `nullptr` in constructor.
- **Memory**: `memory->create()`/`destroy()`/`grow()` — no raw `new`/`delete`.
- **Errors**: `error->all()`, `error->one()`, `error->warning()`; `utils::sfmt()` for messages.
- **Input validation**: all user input checked in `coeff()`/`settings()`.
- **Constants**: `static constexpr` over macros.
- **Naming**: `PascalCase` classes, `snake_case` functions, `UPPER_SNAKE_CASE` constants.

---

## Documentation & Changelog Policy

Stale documentation is treated as a defect. When a change affects behaviour,
CLI options, settings, examples, citations, or author list, update **all three**
in the same commit or PR:

1. **Docs site** (`docs/`, `mkdocs.yml`) — relevant MkDocs pages.
2. **README.md** — top-level README.
3. **CHANGELOG.md** — under `[Unreleased]`, using Keep a Changelog format.

## Commit Conventions

- [Conventional Commits](https://www.conventionalcommits.org/):
`feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`.
- One logical change per commit; under 300 lines of diff when possible.
- Never add `Co-Authored-By` lines to commit messages.

## Working Principles

1. **Scientific rigour first.** Clarify the physical basis before implementing.
  Validate against known limits and reference data.
2. **Reproducibility.** Every numerical claim backed by a re-runnable script or test.
3. **Test-driven for numerics.** Check energy conservation, force symmetry, or
  agreement with analytical expressions before merging force-computation changes.
4. **Performance awareness.** Profile before optimising; document O(N) trade-offs.

## Literature Research

Use available MCP tools to ground decisions in published literature:

- **Google Scholar** (`user-google-scholar`) — keyword search, citation metrics.
- **arXiv** (`user-arxiv-mcp-server`) — preprint search, abstracts, PDFs.

Always cite: authors, journal/venue, year, DOI or arXiv ID.

## Ruff Lint Rules (reference)

`E W F I N UP B SIM TCH RUF PT C4 PIE RET` — see `pyproject.toml` for details.
`E501` (line-too-long) is ignored; ruff formatter handles wrapping.
