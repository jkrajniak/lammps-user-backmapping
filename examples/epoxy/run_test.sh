#!/bin/sh
# Tier A generator test for rim135 epoxy network hybrid prep.
# Copyright (C) 2016 Jakub Krajniak — GPLv3

set -e

ROOT="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
PKG="$(CDPATH= cd -- "$ROOT/../.." && pwd)"
RIM135="${BACKMAP_RIM135:-$ROOT/../../../../bakery/tests/rim135}"

if [ ! -f "$RIM135/settings.xml" ]; then
    echo "rim135 data not found at: $RIM135" >&2
    echo "Clone bakery or set BACKMAP_RIM135" >&2
    exit 1
fi

cd "$RIM135"

uv run --directory "$PKG" backmap-prep build-hybrid "$ROOT/settings.v2.yaml"

COMPARE_OUT=$(uv run --directory "$PKG" backmap-prep compare-topology ref_hyb_topol.top hyb_topol.top || true)
echo "$COMPARE_OUT"
TOPOL=0
echo "$COMPARE_OUT" | grep -q "False" && TOPOL=1

diff ref_hyb_conf.gro hyb_conf.gro
CONF=$?

if [ "$TOPOL" = "0" ] && [ "$CONF" = "0" ]; then
    echo "$0 OK"
    rm -f hyb_conf.gro hyb_topol.top missing_definitions.txt
    exit 0
fi

echo "$0 FAIL (topol=$TOPOL conf=$CONF)"
exit 1
