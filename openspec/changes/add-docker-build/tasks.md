## 1. Docker Build Files

- [x] 1.1 Create `Dockerfile` with multi-stage build (builder + runtime stages), `LAMMPS_VERSION` build arg defaulting to `stable_22Jul2025_update2`, MPI-enabled cmake build with `PKG_BACKMAP=yes`, and `backmap-prep` Python CLI installation
- [x] 1.2 Create `.dockerignore` excluding `.git/`, `build/`, `lammps/`, `__pycache__/`, `.cursor/`, `openspec/`, `*.pyc`, `.ruff_cache/`, `.mypy_cache/`

## 2. Documentation

- [x] 2.1 Add `docs/docker.md` page with build instructions, run examples (volume mounts, MPI), LAMMPS version override, and Singularity/Apptainer conversion guide
- [x] 2.2 Add Docker entry to `mkdocs.yml` navigation
- [x] 2.3 Update `README.md` with a Docker installation/usage section
- [x] 2.4 Add entry to `CHANGELOG.md` under `[Unreleased]` → `Added`
