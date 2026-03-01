## Why

The project currently has only a top-level README with basic installation and usage instructions. Users who want to set up backmapping for a new molecular system lack step-by-step guidance, detailed settings reference, and conceptual explanations of the method. Publishing documentation as GitHub Pages makes the project more discoverable and lowers the barrier to adoption.

## What Changes

- Add a `docs/` directory with MkDocs-based documentation site
- Create a comprehensive YAML settings reference documenting every field, its type, default, and purpose
- Write a tutorial walking through how to set up backmapping for a new molecular system from scratch
- Add a "Getting Started" quickstart guide
- Document the theoretical background (lambda ramp, cross interactions, force weighting)
- Configure GitHub Actions to build and deploy documentation to GitHub Pages
- Add `docs` Make target for local preview

## Capabilities

### New Capabilities
- `documentation-site`: MkDocs-based documentation site with GitHub Pages deployment, settings reference, tutorials, and theoretical background

### Modified Capabilities

(none — no existing spec requirements change)

## Impact

- **New files**: `docs/` directory tree, `mkdocs.yml`, `.github/workflows/docs.yml`
- **Modified files**: `Makefile` (new `docs` target), `README.md` (link to published docs)
- **Dependencies**: `mkdocs`, `mkdocs-material` (dev/docs dependencies only, not runtime)
- **No C++ or Python runtime code changes**
