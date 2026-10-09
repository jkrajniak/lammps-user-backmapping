# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Structural validation for rim135 LAMMPS backmap vs JCC 2018 paper reference data.

Computes C-O and C-N radial distribution functions from a LAMMPS custom dump
(AA production frames) and compares to GROMACS `.xvg` references in
`paper-reverse-mapping-polymer-networks/paper/rim135_small/`.

Also reports total charge, mass density, and bonded distance sanity checks on
`rim135_final.data`.

Usage:
  uv run compare_rim135_structure.py \\
    --dump /path/to/dump.backmap \\
    --gro /path/to/hyb_conf.gro \\
    --final-data /path/to/rim135_final.data \\
    --paper-dir /path/to/paper-reverse-mapping-polymer-networks/paper/rim135_small \\
    [--plot rdf_comparison.png]
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PEAK_POS_TOL_NM = 0.03  # ~0.3 Å
PEAK_HT_TOL = 0.35
L2_TOL = 0.35
R_MAX_NM = 1.0
N_BINS = 500
R_MIN_PEAK_NM = 0.20  # skip bonded peaks when locating first solvation shell


@dataclass
class RDFCurve:
    r_nm: np.ndarray
    g: np.ndarray
    label: str


@dataclass
class DumpFrame:
    step: int
    box_ang: np.ndarray  # Lx, Ly, Lz
    ids: np.ndarray
    xyz_ang: np.ndarray  # (N, 3)


def load_xvg(path: Path) -> RDFCurve:
    rows: list[tuple[float, float]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "@")):
            continue
        parts = line.split()
        rows.append((float(parts[0]), float(parts[1])))
    data = np.array(rows, dtype=float)
    return RDFCurve(r_nm=data[:, 0], g=data[:, 1], label=path.name)


def element_from_atom_name(name: str) -> str | None:
    """Match GROMACS selections name \"C*\", \"O*\", \"N*\" on hybrid atom names."""
    if not name:
        return None
    letter = name[0].upper()
    if letter in {"C", "O", "N", "H"}:
        return letter
    return None


def load_element_map(gro_path: Path) -> dict[int, str]:
    lines = gro_path.read_text().splitlines()
    natoms = int(lines[1].strip())
    mapping: dict[int, str] = {}
    for line in lines[2 : 2 + natoms]:
        name = line[10:15].strip()
        index = int(line[15:20].strip())
        element = element_from_atom_name(name)
        if element is not None:
            mapping[index] = element
    return mapping


def parse_lammps_dump(path: Path, min_step: int = 0) -> list[DumpFrame]:
    frames: list[DumpFrame] = []
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip() != "ITEM: TIMESTEP":
            i += 1
            continue
        step = int(lines[i + 1].strip())
        i += 2
        while i < len(lines) and lines[i].strip() != "ITEM: NUMBER OF ATOMS":
            i += 1
        n_atoms = int(lines[i + 1].strip())
        i += 2
        while i < len(lines) and lines[i].strip() != "ITEM: BOX BOUNDS pp pp pp":
            i += 1
        xlo, xhi = map(float, lines[i + 1].split()[:2])
        ylo, yhi = map(float, lines[i + 2].split()[:2])
        zlo, zhi = map(float, lines[i + 3].split()[:2])
        box = np.array([xhi - xlo, yhi - ylo, zhi - zlo], dtype=float)
        i += 4
        while i < len(lines) and not lines[i].startswith("ITEM: ATOMS"):
            i += 1
        i += 1
        ids = np.empty(n_atoms, dtype=int)
        xyz = np.empty((n_atoms, 3), dtype=float)
        for j in range(n_atoms):
            parts = lines[i + j].split()
            ids[j] = int(parts[0])
            xyz[j, 0] = float(parts[3])
            xyz[j, 1] = float(parts[4])
            xyz[j, 2] = float(parts[5])
        i += n_atoms
        if step >= min_step:
            frames.append(DumpFrame(step=step, box_ang=box, ids=ids, xyz_ang=xyz))
    return frames


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.round(delta / box)


def compute_rdf(
    ref_xyz: np.ndarray,
    sel_xyz: np.ndarray,
    box_ang: np.ndarray,
    *,
    r_max_nm: float = R_MAX_NM,
    n_bins: int = N_BINS,
    chunk: int = 64,
) -> RDFCurve:
    r_max_ang = r_max_nm * 10.0
    dr = r_max_ang / n_bins
    hist = np.zeros(n_bins, dtype=float)
    volume = float(np.prod(box_ang))
    n_sel = len(sel_xyz)
    rho = n_sel / volume

    for start in range(0, len(ref_xyz), chunk):
        ref_chunk = ref_xyz[start : start + chunk]
        delta = ref_chunk[:, None, :] - sel_xyz[None, :, :]
        delta = minimum_image(delta, box_ang)
        dist = np.linalg.norm(delta, axis=2)
        dist = dist.reshape(-1)
        dist = dist[dist > 1e-6]
        bins = (dist / dr).astype(int)
        bins = bins[(bins >= 0) & (bins < n_bins)]
        np.add.at(hist, bins, 1.0)

    n_ref = len(ref_xyz)
    hist /= max(n_ref, 1)
    r_ang = (np.arange(n_bins) + 0.5) * dr
    shell = 4.0 * np.pi * r_ang**2 * dr
    g = hist / (shell * rho + 1e-30)
    return RDFCurve(r_nm=r_ang / 10.0, g=g, label="computed")


def average_rdf(curves: list[RDFCurve]) -> RDFCurve:
    if not curves:
        raise ValueError("no RDF curves to average")
    r = curves[0].r_nm
    g = np.mean(np.stack([c.g for c in curves], axis=0), axis=0)
    return RDFCurve(r_nm=r, g=g, label="computed-mean")


def first_peak(curve: RDFCurve, r_min_nm: float = R_MIN_PEAK_NM) -> tuple[float, float]:
    mask = curve.r_nm > r_min_nm
    if not np.any(mask):
        return 0.0, 0.0
    g = curve.g[mask]
    r = curve.r_nm[mask]
    idx = int(np.argmax(g))
    return float(r[idx]), float(g[idx])


def l2(a: RDFCurve, b: RDFCurve) -> float:
    g_a = np.interp(b.r_nm, a.r_nm, a.g, left=0.0, right=0.0)
    return float(np.sqrt(np.mean((g_a - b.g) ** 2)))


def compare_pair(name: str, computed: RDFCurve, reference: RDFCurve) -> tuple[bool, list[str]]:
    lines: list[str] = []
    ok_all = True
    pos_c, ht_c = first_peak(computed)
    pos_r, ht_r = first_peak(reference)
    pos_err = abs(pos_c - pos_r)
    ok = pos_err <= PEAK_POS_TOL_NM
    ok_all &= ok
    lines.append(
        f"  {'PASS' if ok else 'FAIL'}: {name} peak position "
        f"computed={pos_c:.3f} nm ref={pos_r:.3f} nm Δ={pos_err:.3f} nm (tol={PEAK_POS_TOL_NM})"
    )
    ht_err = abs(ht_c - ht_r) / ht_r if ht_r > 0 else 0.0
    ok = ht_err <= PEAK_HT_TOL
    ok_all &= ok
    lines.append(
        f"  {'PASS' if ok else 'FAIL'}: {name} peak height "
        f"computed={ht_c:.3f} ref={ht_r:.3f} rel err={ht_err:.2f} (tol={PEAK_HT_TOL})"
    )
    l2_val = l2(computed, reference)
    ok = l2_val <= L2_TOL
    ok_all &= ok
    lines.append(f"  {'PASS' if ok else 'FAIL'}: {name} L2(g)={l2_val:.3f} (tol={L2_TOL})")
    return ok_all, lines


def _parse_lammps_data_simple(
    data_path: Path,
) -> tuple[
    dict[int, tuple[float, float, float, int, int, int]],
    list[tuple[int, int]],
    np.ndarray,
    dict[int, float],
    dict[int, float],
    dict[int, int],
]:
    """Return atoms, bonds, box Lxyz, masses, charges, type_by_id."""
    atoms: dict[int, tuple[float, float, float, int, int, int]] = {}
    bonds: list[tuple[int, int]] = []
    masses: dict[int, float] = {}
    charges: dict[int, float] = {}
    type_by_id: dict[int, int] = {}
    box = np.zeros(3)
    section: str | None = None
    for line in data_path.read_text().splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("Masses"):
            section = "masses"
            continue
        if stripped.startswith("Atoms"):
            section = "atoms"
            continue
        if stripped.startswith("Bonds"):
            section = "bonds"
            continue
        if stripped.startswith(("Velocities", "Angles", "Dihedrals", "Impropers", "Pair Coeffs")):
            section = None
            continue
        if "xlo" in stripped:
            parts = stripped.split()
            box[0] = float(parts[1]) - float(parts[0])
        elif "ylo" in stripped:
            parts = stripped.split()
            box[1] = float(parts[1]) - float(parts[0])
        elif "zlo" in stripped:
            parts = stripped.split()
            box[2] = float(parts[1]) - float(parts[0])
        elif section == "masses":
            if stripped[0].isalpha():
                section = None
                continue
            parts = stripped.split()
            masses[int(parts[0])] = float(parts[1])
        elif section == "atoms":
            if stripped[0].isalpha():
                section = None
                continue
            parts = stripped.split()
            aid = int(parts[0])
            type_by_id[aid] = int(parts[2])
            charge = float(parts[3])
            charges[aid] = charge
            if len(parts) >= 10:
                x, y, z = map(float, parts[4:7])
                ix, iy, iz = map(int, parts[7:10])
            else:
                x, y, z = map(float, parts[4:7])
                ix = iy = iz = 0
            atoms[aid] = (x, y, z, ix, iy, iz)
        elif section == "bonds":
            if stripped[0].isalpha():
                section = None
                continue
            parts = stripped.split()
            bonds.append((int(parts[2]), int(parts[3])))
    return atoms, bonds, box, masses, charges, type_by_id


def parse_final_data_checks(data_path: Path) -> dict[str, float]:
    atoms, bonds, box, masses, charges, type_by_id = _parse_lammps_data_simple(data_path)
    total_charge = sum(charges.values())
    volume_nm3 = float(np.prod(box)) / 1000.0

    def atom_xyz(aid: int) -> np.ndarray:
        x, y, z, ix, iy, iz = atoms[aid]
        return np.array([x + ix * box[0], y + iy * box[1], z + iz * box[2]], dtype=float)

    max_bond = 0.0
    max_eucl = 0.0
    for i, j in bonds:
        pi = atom_xyz(i)
        pj = atom_xyz(j)
        delta = pj - pi
        delta_mi = minimum_image(delta, box)
        max_bond = max(max_bond, float(np.linalg.norm(delta_mi)))
        max_eucl = max(max_eucl, float(np.linalg.norm(delta)))

    total_mass = sum(masses.get(type_by_id.get(aid, 0), 0.0) for aid in atoms)
    density_g_cm3 = total_mass / (volume_nm3 * 602.214076)

    return {
        "total_charge": total_charge,
        "max_bond_ang": max_bond,
        "max_eucl_bond_ang": max_eucl,
        "volume_nm3": volume_nm3,
        "density_g_cm3": density_g_cm3,
        "n_atoms": float(len(atoms)),
        "n_bonds": float(len(bonds)),
    }


def plot_rdgs(curves: list[tuple[RDFCurve, RDFCurve, str]], output: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(curves)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), squeeze=False)
    for ax, (calc, ref, title) in zip(axes[0], curves, strict=True):
        ax.plot(calc.r_nm, calc.g, color="C0", lw=2.4, label="backmapped")
        ax.plot(ref.r_nm, ref.g, "--", color="C1", lw=2.4, label="reference")
        ax.set_xlim(0, R_MAX_NM)
        ax.set_xlabel("r (nm)", fontsize=14)
        ax.set_ylabel("g(r)", fontsize=14)
        ax.set_title(title, fontsize=15)
        ax.tick_params(axis="both", labelsize=12)
        ax.legend(fontsize=12)
    fig.tight_layout()
    # Line widths, font sizes and vector output match the conventions applied
    # to the dodecane (PR #7) and melamine_network (PR #11) figures. Writing a
    # .pdf path yields vector output; dpi only affects raster formats.
    fig.savefig(output)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description="Compare rim135 LAMMPS output vs paper RDF")
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--gro", type=Path, required=True)
    parser.add_argument("--final-data", type=Path, required=True)
    parser.add_argument("--paper-dir", type=Path, required=True)
    parser.add_argument(
        "--min-step",
        type=int,
        default=10000,
        help="First dump timestep for AT production averaging (default: 10000)",
    )
    parser.add_argument("--plot", type=Path, default=None)
    args = parser.parse_args()

    for path in (args.dump, args.gro, args.final_data, args.paper_dir):
        if not path.exists():
            print(f"ERROR: missing {path}", file=sys.stderr)
            return 1

    print("=== rim135 structural validation vs paper reference ===\n")

    elem_map = load_element_map(args.gro)
    frames = parse_lammps_dump(args.dump, min_step=args.min_step)
    if not frames:
        print("ERROR: no dump frames at or after min-step", file=sys.stderr)
        return 1
    print(f"Dump: {len(frames)} frames from step {frames[0].step} to {frames[-1].step}")

    pair_specs = [
        ("C", "O", "C-O"),
        ("C", "N", "C-N"),
    ]
    computed_by_pair: dict[str, list[RDFCurve]] = {label: [] for _, _, label in pair_specs}

    for frame in frames:
        elements = np.array([elem_map.get(int(aid), "") for aid in frame.ids], dtype=object)
        for ref_elem, sel_elem, label in pair_specs:
            ref_mask = elements == ref_elem
            sel_mask = elements == sel_elem
            if not ref_mask.any() or not sel_mask.any():
                continue
            computed_by_pair[label].append(
                compute_rdf(frame.xyz_ang[ref_mask], frame.xyz_ang[sel_mask], frame.box_ang)
            )

    ref_aa = {
        "C-O": args.paper_dir / "aa/rdf/rdf_s6_0.000001_C_O.xvg",
        "C-N": args.paper_dir / "aa/rdf/rdf_s6_0.000001_C_N.xvg",
    }
    ref_ua = {
        "C-O": args.paper_dir / "ua/rdf/rdf_15_0.000001_C_O.xvg",
        "C-N": args.paper_dir / "ua/rdf/rdf_15_0.000001_C_N.xvg",
    }

    n_pass = 0
    n_fail = 0
    plot_pairs: list[tuple[RDFCurve, RDFCurve, str]] = []

    for _, _, label in pair_specs:
        if not computed_by_pair[label]:
            print(f"\n--- {label}: no atoms matched ---")
            n_fail += 3
            continue
        computed = average_rdf(computed_by_pair[label])
        aa_ref = load_xvg(ref_aa[label])
        print(f"\n--- {label} vs paper AA (rdf_s6_0.000001) ---")
        ok, lines = compare_pair(label, computed, aa_ref)
        for line in lines:
            print(line)
        n_pass += sum(1 for line in lines if line.startswith("  PASS"))
        n_fail += sum(1 for line in lines if line.startswith("  FAIL"))
        plot_pairs.append((computed, aa_ref, f"{label} vs AA"))

        if ref_ua[label].exists():
            ua_ref = load_xvg(ref_ua[label])
            print(f"\n--- {label} vs paper UA (rdf_15_0.000001) ---")
            ok_u, lines_u = compare_pair(label, computed, ua_ref)
            for line in lines_u:
                print(line)
            n_pass += sum(1 for line in lines_u if line.startswith("  PASS"))
            n_fail += sum(1 for line in lines_u if line.startswith("  FAIL"))

    checks = parse_final_data_checks(args.final_data)
    print("\n--- final structure checks (rim135_final.data) ---")
    charge_ok = abs(checks["total_charge"]) < 1e-3
    bond_ok = checks["max_bond_ang"] < 20.0
    print(
        f"  {'PASS' if charge_ok else 'FAIL'}: total charge = {checks['total_charge']:.6f} e (tol 1e-3)"
    )
    print(
        f"  {'PASS' if bond_ok else 'FAIL'}: max min-image bond = {checks['max_bond_ang']:.2f} Å (tol 20 Å)"
    )
    print(f"  INFO: max euclidean bond = {checks['max_eucl_bond_ang']:.2f} Å")
    print(
        f"  INFO: density ≈ {checks['density_g_cm3']:.3f} g/cm³ (box {checks['volume_nm3']:.1f} nm³)"
    )
    n_pass += int(charge_ok) + int(bond_ok)
    n_fail += int(not charge_ok) + int(not bond_ok)

    if args.plot and plot_pairs:
        plot_rdgs(plot_pairs, args.plot)
        print(f"\nPlot saved to {args.plot}")

    print(f"\n{'=' * 50}")
    print(f"Results: {n_pass} passed, {n_fail} failed")
    print(
        "\nNote: only 10 ps AT production (10 dump frames) — RDF statistics are noisy vs "
        "long GROMACS runs in the 2018 paper (-b 2500 ps)."
    )
    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
