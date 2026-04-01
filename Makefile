.PHONY: help install install-dev install-hooks lint format typecheck test test-cov clean docs docs-serve

help: ## Show this help
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
		awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

# ── Setup ───────────────────────────────────────────────────────────

install: ## Install Python package
	uv sync

install-dev: ## Install Python package with dev dependencies
	uv sync --extra dev

install-hooks: install-dev ## Install pre-commit hooks into the repo
	uv run pre-commit install

# ── Quality ─────────────────────────────────────────────────────────

lint: ## Run ruff linter
	uv run ruff check python/src/ python/tests/

format: ## Run ruff formatter (in-place)
	uv run ruff format python/src/ python/tests/

format-check: ## Check formatting without modifying files
	uv run ruff format --check python/src/ python/tests/

typecheck: ## Run mypy type checker
	uv run mypy python/src/

# ── Tests ───────────────────────────────────────────────────────────

test: ## Run pytest
	uv run pytest

test-cov: ## Run pytest with coverage report
	uv run pytest --cov=backmap_prep --cov-report=term-missing

# ── Pre-commit ──────────────────────────────────────────────────────

pre-commit: ## Run all pre-commit hooks on staged files
	uv run pre-commit run

pre-commit-all: ## Run all pre-commit hooks on every file
	uv run pre-commit run --all-files

# ── Documentation ───────────────────────────────────────────────────

docs: ## Build documentation site
	NO_MKDOCS_2_WARNING=1 uv run --with mkdocs --with mkdocs-material mkdocs build --strict

docs-serve: ## Serve documentation locally for preview
	NO_MKDOCS_2_WARNING=1 uv run --with mkdocs --with mkdocs-material mkdocs serve

# ── Cleanup ─────────────────────────────────────────────────────────

clean: ## Remove build artifacts and caches
	rm -rf python/.mypy_cache python/.ruff_cache python/.pytest_cache
	find python -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find python -type d -name '*.egg-info' -exec rm -rf {} + 2>/dev/null || true
