# LAMMPS Backmapping Package — Project Constitution

## Overview

This project implements time-dependent backmapping (reverse mapping) from
coarse-grained to atomistic resolution as a LAMMPS user package. It consists of:

- **C++ LAMMPS styles** (`src/`): fix, pair, bond, and angle styles that drive
  the lambda-ramp resolution change.
- **Python tooling** (`python/`): the `backmap-prep` CLI that generates LAMMPS
  input files from GROMACS topologies and coordinate files.

Reference:
> Krajniak et al., *J. Chem. Theory Comput.* 2016, 12, 5549–5562.
> DOI: 10.1021/acs.jctc.6b00595

---

## Tech Stack

| Layer | Technology | Notes |
|-------|-----------|-------|
| Simulation engine | LAMMPS (C++17) | User package compiled into LAMMPS |
| Input generator | Python ≥ 3.10 | `backmap-prep` CLI |
| Build (C++) | CMake | Integrated with LAMMPS build system |
| Build (Python) | Hatchling + uv | `uv` for all dependency management |
| Linting / Formatting | ruff | Python linting and formatting |
| Type checking | mypy | Strict mode |
| Testing | pytest | All Python code must have unit tests |
| Pre-commit | pre-commit | Enforces style before every commit |

---

## C++ Guidelines (LAMMPS Styles)

### Style & Conventions

- Follow the [LAMMPS coding style](https://docs.lammps.org/Modify_style.html):
  two-space indentation, `ClassName` for classes, `lower_snake_case` for
  variables and functions.
- Every source file must include the LAMMPS copyright header and the package
  author attribution.
- Use C++17 features where LAMMPS supports them; avoid C++20-only constructs.

### Memory Safety

- **No raw `new`/`delete` for arrays.** Use LAMMPS `memory->create()` /
  `memory->destroy()` or `memory->grow()` for per-atom and per-type arrays.
- Every allocation in the constructor or `allocate()` must have a matching
  deallocation in the destructor.
- Guard every pointer with a `nullptr` initialisation in the constructor
  initialiser list.
- When resizing arrays on `atom->nmax` changes, always use `memory->grow()`
  with the new size; never assume the old pointer is still valid.
- Validate all user inputs in `coeff()` / `settings()` and call
  `error->all()` with a clear message on failure.

### Robustness

- Check return values from LAMMPS utility functions (`force->numeric()`,
  `force->inumeric()`, etc.) and handle errors immediately.
- Use `comm->forward_comm(this)` / `comm->reverse_comm(this)` correctly;
  document which per-atom quantities are communicated and their size.
- In `compute()`, always guard against zero-length neighbor lists and
  atoms that have migrated out of the local domain.

---

## Python Guidelines

### Modern Python

- Target Python ≥ 3.10. Use modern syntax: `match` statements,
  `X | Y` union types, `list[T]` / `dict[K, V]` built-in generics.
- Use `def` for pure functions, `async def` only when I/O-bound.
- All functions must have complete type annotations.
- Use Pydantic v2 models for structured data (input schemas, configs).
- Prefer functional, declarative code; use classes only when state is
  genuinely needed.
- Use descriptive variable names with auxiliary verbs where appropriate
  (e.g., `is_valid`, `has_bonds`).

### Dependency Management

- **Use `uv` exclusively.** Never use `pip` or `pip-tools` directly.
  ```bash
  uv add <package>      # add dependency
  uv remove <package>   # remove dependency
  uv sync               # reinstall from lock
  uv run script.py      # run with correct env
  ```

### Testing

- **Every Python module must have corresponding pytest unit tests.**
- Tests live alongside the source in a `tests/` directory mirroring the
  package structure: `python/tests/test_<module>.py`.
- Test files must be importable and self-contained (no reliance on
  working-directory tricks).
- Use `pytest` fixtures for shared setup; prefer parametrised tests for
  covering multiple input variants.
- Aim for meaningful coverage of edge cases, not just happy paths.

### Linting, Formatting & Type Checking

- **ruff** handles both linting and formatting. Configuration lives in
  `pyproject.toml` under `[tool.ruff]`.
- **mypy** in strict mode checks all type annotations. Configuration
  lives in `pyproject.toml` under `[tool.mypy]`.
- Pre-commit hooks enforce both on every commit; CI must also run them.

### Error Handling

- Guard clauses at the top of functions; early returns for invalid input.
- Raise domain-specific exceptions with clear messages.
- Never silently swallow exceptions.

---

## Commit Conventions

- Use [Conventional Commits](https://www.conventionalcommits.org/):
  `feat:`, `fix:`, `docs:`, `refactor:`, `test:`, `chore:`.
- Keep commits atomic: one logical change per commit.

---

## Documentation & Change Log

- **CHANGELOG.md** — Every user-visible change (features, fixes, breaking
  changes) must be recorded in `CHANGELOG.md` at the repository root. Follow
  the [Keep a Changelog](https://keepachangelog.com/) format with sections
  `Added`, `Changed`, `Deprecated`, `Removed`, `Fixed`, and `Security` under
  each release heading. Update the `[Unreleased]` section as part of the same
  commit or PR that introduces the change.
- **README.md** — The top-level `README.md` must always reflect the current
  state of the project. When a change affects installation steps, usage
  instructions, supported features, CLI options, or repository layout, the
  README must be updated in the same commit or PR. Stale documentation is
  treated as a defect.

---

## Repository Layout

```
src/                    # C++ LAMMPS styles
python/
  src/backmap_prep/     # Python package source
  tests/                # pytest unit tests
  pyproject.toml        # Python project metadata & tool config
examples/               # Example simulations
openspec/               # Specifications and change tracking
```
