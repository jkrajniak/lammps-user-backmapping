## Context

The LAMMPS Backmapping Package currently documents itself through a README.md with installation, usage, and development instructions. The YAML settings file schema is defined in Python (Pydantic models in `schema.py`) but has no user-facing reference documentation. Users must read the example `settings.yaml` and source code to understand available options. There is no tutorial explaining how to set up backmapping for a new molecular system.

## Goals / Non-Goals

**Goals:**
- Provide a browsable, searchable documentation site published via GitHub Pages
- Document every YAML settings field with type, default, and description
- Walk users through setting up a new system end-to-end (tutorial)
- Explain the theoretical background of the backmapping method
- Automate deployment so docs stay in sync with the codebase

**Non-Goals:**
- API reference for internal Python functions (may come later)
- C++ developer documentation (LAMMPS style internals)
- Hosting on a custom domain
- Internationalization

## Decisions

### 1. Static site generator: MkDocs with Material theme

**Choice**: MkDocs + mkdocs-material

**Rationale**: MkDocs is Python-native (fits the project's tooling), has first-class GitHub Pages deployment support via `mkdocs gh-deploy`, and the Material theme provides search, dark mode, admonitions, and content tabs out of the box. Alternatives considered:
- Sphinx: more powerful for API docs but heavier configuration and RST-based
- Jekyll: Ruby-based, doesn't fit the Python ecosystem
- Docusaurus: React-based, overkill for a scientific package

### 2. Documentation structure

```
docs/
  index.md                  # Landing page / overview
  getting-started.md        # Quickstart guide
  theory.md                 # Theoretical background
  settings-reference.md     # Complete YAML settings reference
  tutorial-new-system.md    # Step-by-step tutorial
  components/
    fix-backmap.md          # fix backmap documentation
    pair-backmap.md         # pair_style backmap documentation
    bond-styles.md          # bond_style backmap/* documentation
    angle-styles.md         # angle_style backmap/* documentation
  cli/
    backmap-prep.md         # CLI usage reference
mkdocs.yml                  # Site configuration
```

**Rationale**: Separating reference, tutorial, and theory content follows the Diátaxis documentation framework (tutorials, how-to guides, reference, explanation).

### 3. Deployment: GitHub Actions

**Choice**: GitHub Actions workflow deploying to `gh-pages` branch on push to `main`.

**Rationale**: Standard approach, zero external dependencies, uses `mkdocs gh-deploy` which handles the branch creation and CNAME setup automatically.

### 4. Dependency management

MkDocs and plugins will be added as `docs` extras in `pyproject.toml` or managed as a standalone `docs/requirements.txt`. Since the project uses `uv`, a simple `uv pip install mkdocs mkdocs-material` in the CI workflow is sufficient.

**Choice**: Add a `docs/requirements.txt` to keep documentation dependencies separate from the Python package dependencies.

## Risks / Trade-offs

- **[Docs drift]** Documentation may become stale if the settings schema changes without updating the reference. → Mitigation: Add a CI check or note in the project constitution to update docs alongside schema changes.
- **[Build dependency]** MkDocs adds a build-time dependency. → Mitigation: It's only needed for docs builds, not runtime. Isolated in `docs/requirements.txt`.
- **[Settings reference maintenance]** The settings reference is manually written, not auto-generated from Pydantic models. → Mitigation: Keep it manageable by structuring it to mirror the schema. Auto-generation can be added later.

## Open Questions

- Should the settings reference eventually be auto-generated from Pydantic model docstrings? (Deferred — manual is fine for initial version.)
