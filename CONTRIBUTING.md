# Contributing

Thanks for your interest in contributing to the LAMMPS Backmapping Package!

## Getting Started

1. Fork the repository and clone your fork.
2. Install the C++ package into a local LAMMPS tree (see [README](README.md#installation)).
3. Install the Python development environment:

   ```bash
   make install-dev
   make install-hooks
   ```

## Development Workflow

1. Create a feature branch from `main`.
2. Make your changes.
3. Run the checks before submitting:

   ```bash
   # C++ formatting
   clang-format --style=file --fallback-style=Google -i src/*.cpp src/*.h

   # Python checks
   make lint
   make format
   make typecheck
   make test
   ```

4. Open a pull request against `main`.

## Code Style

- **C++**: Follow the existing LAMMPS coding style. Use `clang-format` with
  the project configuration (Google fallback).
- **Python**: Code is formatted and linted with [Ruff](https://docs.astral.sh/ruff/)
  and type-checked with [mypy](https://mypy-lang.org/) in strict mode. See
  `python/pyproject.toml` for the full configuration.

## Reporting Issues

Please open a GitHub issue with:

- A clear description of the problem or feature request.
- Steps to reproduce (for bugs).
- Relevant LAMMPS version and OS information.

## License

By contributing, you agree that your contributions will be licensed under
GPL-3.0-or-later, consistent with the rest of the project.
