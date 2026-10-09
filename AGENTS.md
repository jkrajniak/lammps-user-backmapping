# AGENTS.md

## Project Overview

Adaptive-resolution backmapping from coarse-grained to atomistic resolution,
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

---

## Build, Lint, and Test Commands

### Python (backmap-prep)

Targets are in the `Makefile` (`make install-dev`, `make lint`, `make test`, `make docs`).
Single test: `uv run pytest python/tests/test_foo.py::test_bar`.

### C++ (LAMMPS styles)

Build and run LAMMPS only on the remote VM, never locally. Lint with
`clang-format --style=file --fallback-style=Google src/*.cpp src/*.h`.

---

## Code Style — Python

- **Python ≥ 3.10**, full type annotations, Pydantic v2 for schemas.
- **Dependencies**: `uv` exclusively — never `pip`. Use `uv add`, `uv sync`.
- **Structure**: functions under 50 lines; `__all__` in public modules; Pydantic models or
dataclasses for data containers.
- **Errors**: guard clauses at top, early returns, specific exceptions (never bare `except:`).
- **Tests**: every module needs a `python/tests/test_<module>.py`; use parametrised tests
and fixtures; cover edge cases.

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
