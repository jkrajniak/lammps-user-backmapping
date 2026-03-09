#!/usr/bin/env bash
# Optional validation: run backmap-prep for at least one large-scale example.
# Does NOT run LAMMPS (simulations). Use for CI or manual checks.
# Exit 0 if backmap-prep succeeds; non-zero otherwise.

set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
EXAMPLE="${1:-dodecane}"
SETTINGS="$REPO_ROOT/examples/$EXAMPLE/large/settings.yaml"

if [[ ! -f "$SETTINGS" ]]; then
  echo "Settings not found: $SETTINGS" >&2
  exit 1
fi

cd "$REPO_ROOT/python"
uv run backmap-prep "$SETTINGS"
echo "OK: backmap-prep succeeded for $EXAMPLE large-scale"
