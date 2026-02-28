#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG_SRC_DIR="${SCRIPT_DIR}/src"
PKG_NAME="BACKMAP"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

info()  { echo -e "${BLUE}[INFO]${NC}  $*"; }
ok()    { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}  $*"; }
err()   { echo -e "${RED}[ERROR]${NC} $*" >&2; }

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS] <lammps-dir>

Install the LAMMPS backmapping package into a LAMMPS source tree.

Arguments:
  lammps-dir          Path to the LAMMPS root directory (containing src/ and cmake/)

Options:
  -l, --link          Symlink sources instead of copying (dev mode).
                      Edits in either tree are immediately visible to the other.
  -u, --uninstall     Remove the package from LAMMPS
  -h, --help          Show this help message

The script copies (or symlinks with --link) package sources into
lammps/src/${PKG_NAME}/ and configures whichever build system is detected:

  CMake            Patches cmake/CMakeLists.txt to register the package.
                   Rebuild with:  cmake -D PKG_${PKG_NAME}=yes <build-dir>

  Traditional make The bundled Install.sh is used by "make yes-backmap".
                   Rebuild with:  cd src && make yes-backmap && make mpi
EOF
    exit "${1:-0}"
}

detect_build_system() {
    local lammps_dir="$1"
    local has_cmake=false
    local has_make=false

    if [[ -f "${lammps_dir}/cmake/CMakeLists.txt" ]]; then
        has_cmake=true
    fi

    if [[ -d "${lammps_dir}/src/MAKE" ]]; then
        has_make=true
    fi

    if $has_cmake && $has_make; then
        echo "both"
    elif $has_cmake; then
        echo "cmake"
    elif $has_make; then
        echo "make"
    else
        echo "unknown"
    fi
}

patch_cmake() {
    local cmake_file="$1"

    if grep -q "^  ${PKG_NAME}\$" "$cmake_file" 2>/dev/null; then
        info "CMake: ${PKG_NAME} already in STANDARD_PACKAGES — skipping patch"
        return 0
    fi

    # Find the STANDARD_PACKAGES block boundaries
    local list_start
    list_start=$(grep -n "^set(STANDARD_PACKAGES" "$cmake_file" | head -1 | cut -d: -f1)
    if [[ -z "$list_start" ]]; then
        warn "CMake: could not find STANDARD_PACKAGES list — manual patching required"
        return 1
    fi

    local list_end
    list_end=$(awk -v s="$list_start" 'NR > s && /\)/ {print NR; exit}' "$cmake_file")
    if [[ -z "$list_end" ]]; then
        warn "CMake: could not locate end of STANDARD_PACKAGES list — manual patching required"
        return 1
    fi

    # Insert alphabetically: find first entry that sorts after PKG_NAME
    local insert_at=""
    while IFS= read -r line; do
        local lineno entry
        lineno=$(echo "$line" | cut -d: -f1)
        entry=$(echo "$line" | sed 's/^[0-9]*://' | tr -d ' ')
        if [[ -n "$entry" && "$entry" > "$PKG_NAME" ]]; then
            insert_at="$lineno"
            break
        fi
    done < <(sed -n "$((list_start+1)),$((list_end-1))p" "$cmake_file" | grep -n '^\s\+[A-Z]' | awk -v off="$list_start" -F: '{printf "%d:%s\n", $1+off, $2}')

    # Fallback: insert before closing paren
    if [[ -z "$insert_at" ]]; then
        insert_at="$list_end"
    fi

    sed -i.bak "${insert_at}i\\
  ${PKG_NAME}
" "$cmake_file"
    rm -f "${cmake_file}.bak"

    ok "CMake: added ${PKG_NAME} to STANDARD_PACKAGES in cmake/CMakeLists.txt"
}

unpatch_cmake() {
    local cmake_file="$1"

    if ! grep -q "^  ${PKG_NAME}\$" "$cmake_file" 2>/dev/null; then
        info "CMake: ${PKG_NAME} not in STANDARD_PACKAGES — nothing to remove"
        return 0
    fi

    sed -i.bak "/^  ${PKG_NAME}\$/d" "$cmake_file"
    rm -f "${cmake_file}.bak"

    ok "CMake: removed ${PKG_NAME} from STANDARD_PACKAGES"
}

install_package() {
    local lammps_dir="$1"
    local use_links="$2"
    local dest_dir="${lammps_dir}/src/${PKG_NAME}"

    if [[ -d "$dest_dir" ]]; then
        warn "Package directory already exists: ${dest_dir}"
        info "Updating files..."
    else
        mkdir -p "$dest_dir"
    fi

    if [[ "$use_links" == "true" ]]; then
        for f in "${PKG_SRC_DIR}"/*.cpp "${PKG_SRC_DIR}"/*.h; do
            ln -sf "$f" "$dest_dir/"
        done
        ln -sf "${PKG_SRC_DIR}/Install.sh" "$dest_dir/"
        ln -sf "${PKG_SRC_DIR}/CMakeLists.txt" "$dest_dir/"
        ln -sf "${PKG_SRC_DIR}/README" "$dest_dir/"
        ok "Symlinked package sources to ${dest_dir}"
    else
        cp "${PKG_SRC_DIR}"/*.cpp "${PKG_SRC_DIR}"/*.h "$dest_dir/"
        cp "${PKG_SRC_DIR}/Install.sh" "$dest_dir/"
        cp "${PKG_SRC_DIR}/CMakeLists.txt" "$dest_dir/"
        cp "${PKG_SRC_DIR}/README" "$dest_dir/"
        ok "Copied package sources to ${dest_dir}"
    fi

    local build_sys
    build_sys=$(detect_build_system "$lammps_dir")

    local pkg_lower
    pkg_lower=$(echo "$PKG_NAME" | tr '[:upper:]' '[:lower:]')

    if [[ "$build_sys" == "cmake" || "$build_sys" == "both" ]]; then
        patch_cmake "${lammps_dir}/cmake/CMakeLists.txt"
    fi

    echo ""
    case "$build_sys" in
        both)
            ok "Both build systems detected. Next steps:"
            echo ""
            echo "  CMake (recommended):"
            echo "    cd <build-dir>"
            echo "    cmake -D PKG_${PKG_NAME}=yes <lammps-dir>/cmake"
            echo "    cmake --build ."
            echo ""
            echo "  Traditional make:"
            echo "    cd ${lammps_dir}/src"
            echo "    make yes-${pkg_lower}"
            echo "    make mpi"
            ;;
        cmake)
            ok "CMake build system detected. Next steps:"
            echo ""
            echo "    cd <build-dir>"
            echo "    cmake -D PKG_${PKG_NAME}=yes <lammps-dir>/cmake"
            echo "    cmake --build ."
            ;;
        make)
            ok "Traditional make build system detected. Next steps:"
            echo ""
            echo "    cd ${lammps_dir}/src"
            echo "    make yes-${pkg_lower}"
            echo "    make mpi"
            ;;
        *)
            warn "Could not detect LAMMPS build system"
            echo ""
            echo "  You can still build manually:"
            echo "    CMake:  cmake -D PKG_${PKG_NAME}=yes <lammps-dir>/cmake"
            echo "    Make:   cd src && make yes-${pkg_lower} && make mpi"
            ;;
    esac
}

uninstall_package() {
    local lammps_dir="$1"
    local dest_dir="${lammps_dir}/src/${PKG_NAME}"

    if [[ -d "$dest_dir" ]]; then
        if [[ -f "${dest_dir}/Install.sh" ]]; then
            info "Running Install.sh uninstall..."
            (cd "$dest_dir" && bash Install.sh 0)
        fi
        rm -rf "$dest_dir"
        ok "Removed package directory: ${dest_dir}"
    else
        info "Package directory not found: ${dest_dir} — nothing to remove"
    fi

    if [[ -f "${lammps_dir}/cmake/CMakeLists.txt" ]]; then
        unpatch_cmake "${lammps_dir}/cmake/CMakeLists.txt"
    fi

    ok "Uninstall complete"
}

# ── Argument parsing ─────────────────────────────────────────────────

MODE="install"
USE_LINKS="false"
LAMMPS_DIR=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -l|--link)
            USE_LINKS="true"
            shift
            ;;
        -u|--uninstall)
            MODE="uninstall"
            shift
            ;;
        -h|--help)
            usage 0
            ;;
        -*)
            err "Unknown option: $1"
            usage 1
            ;;
        *)
            if [[ -n "$LAMMPS_DIR" ]]; then
                err "Unexpected argument: $1"
                usage 1
            fi
            LAMMPS_DIR="$1"
            shift
            ;;
    esac
done

if [[ -z "$LAMMPS_DIR" ]]; then
    err "Missing required argument: lammps-dir"
    usage 1
fi

LAMMPS_DIR="$(cd "$LAMMPS_DIR" && pwd)"

if [[ ! -d "${LAMMPS_DIR}/src" ]]; then
    err "Not a valid LAMMPS directory (missing src/): ${LAMMPS_DIR}"
    exit 1
fi

# ── Run ──────────────────────────────────────────────────────────────

info "LAMMPS directory: ${LAMMPS_DIR}"
info "Build system:     $(detect_build_system "$LAMMPS_DIR")"
if [[ "$USE_LINKS" == "true" ]]; then
    info "Install mode:     symlink (dev)"
else
    info "Install mode:     copy"
fi
echo ""

case "$MODE" in
    install)   install_package "$LAMMPS_DIR" "$USE_LINKS" ;;
    uninstall) uninstall_package "$LAMMPS_DIR" ;;
esac
