#!/usr/bin/env bash
# Download the third-party force-field files of this example at pinned commits
# and check their SHA-256 (see SOURCES.md for origin and licence).
set -euo pipefail
cd "$(dirname "$0")"
M=784591ebdc91d762ed4df986c4650546c938f776   # marrink-lab/martini-forcefields (Apache-2.0)
C=5fd704204067782e2ed01bec9c36491f088c12c4   # owenvickery/cg2at (GPL-3.0)
get() {  # url dest sha256
  [ -f "$2" ] || curl -fsSL -o "$2" "$1"
  echo "$3  $2" | shasum -a 256 -c -
}
MF=https://raw.githubusercontent.com/marrink-lab/martini-forcefields/$M/martini_forcefields/regular/v3.0.0/gmx_files
FF=https://raw.githubusercontent.com/owenvickery/cg2at/$C/database/forcefields/amber99sb-ildn_slipids.ff
mkdir -p amber99sb-ildn_slipids.ff
get $MF/martini_v3.0.0.itp martini_v3.0.0.itp 0ecbf0ff334b2ec21dc5e3151aca42167b977bc2aed41eb2a44bd8c997754894
get $MF/martini_v3.0.0_phospholipids_v1.itp martini_v3.0.0_phospholipids_v1.itp 21fb2eb6d8a1d48854d179dc71981c1210f608486efb3e6644347391ba9cf392
get $MF/martini_v3.0.0_solvents_v1.itp martini_v3.0.0_solvents_v1.itp e6a78414b317c38a1a816eded4aee6082a5da1ed229a126359f0cc92314dcb4e
get $FF/forcefield.itp amber99sb-ildn_slipids.ff/forcefield.itp a73df8ee82c58f106364276622f7bfab8130cdc2cb698342030ef84545569d0b
get $FF/ffnonbonded.itp amber99sb-ildn_slipids.ff/ffnonbonded.itp 18c9513cde1c6b5baa8dc79716cdaed68997f4042dd4c40fb570810b51382803
get $FF/ffbonded.itp amber99sb-ildn_slipids.ff/ffbonded.itp e872e75bc8a3a106e4156eeb0566775460b0609ef8c81a113dc813b2d675f58b
get $FF/gbsa.itp amber99sb-ildn_slipids.ff/gbsa.itp 3325ccf477d38d5c9bbf70810e6386c74951791bb0307e740733168627f63def
