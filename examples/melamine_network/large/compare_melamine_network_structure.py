# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Structural validation for the crosslinked MF network LAMMPS backmap vs
bakery's real crosslinked-network reference RDFs.

Adapted from `examples/melamine/large/compare_melamine_structure.py`. The
reference `.xvg` files (`ref/ref_C_C.xvg`, `ref_C_N.xvg`, `ref_O_H.xvg`,
vendored unchanged from that example's `ref/`) were confirmed in
`research/experiments/20260802_melamine-bakery-rerun.md` to already be
bakery's own crosslinked `network_backmapping/mf` reference data -- they
were a chemical-system mismatch when compared against the uncrosslinked
`examples/melamine/` example, not against this one. This is the first
genuinely like-for-like comparison against them.

Only the AT-type -> element mapping differs from the uncrosslinked
example's script: `melamine_network`'s own `backmap-prep` build assigned a
different internal type order (types 6/7 are swapped -- verified directly
against this example's own `Masses` section, not assumed), so a separate,
independently-checked mapping is used rather than the uncrosslinked
example's hardcoded one. Everything else (RDF computation, peak/L2
tolerances, C-C bonded/melt-shell windows) is unchanged.

Usage:
  uv run compare_melamine_network_structure.py \\
    --dump dump.at_network_prod \\
    --final-data melamine_network_hybrid_final.data \\
    --ref-dir ../../melamine/large/ref \\
    --min-step 12000 \\
    --plot rdf_comparison.png
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PEAK_POS_TOL_NM = 0.03  # 0.3 Å
PEAK_HT_TOL = 0.20
L2_TOL = 0.20
R_MAX_NM = 1.2
N_BINS = 600
R_MIN_PEAK_NM = 0.28  # default melt-shell lower bound (cross pairs)

# C-C windows: paper ref has bonded ~0.23 nm and melt double-peak ~0.42-0.48 nm
CC_BONDED_R_NM = (0.20, 0.28)
CC_MELT_R_NM = (0.35, 0.55)
LONG_RANGE_R_NM = 1.0  # diagnostic; GROMACS C*-C* ref plateaus ~1.5

# melamine_network's own AT-type -> element map, verified against this
# example's melamine_network.data Masses section (NOT copied from the
# uncrosslinked example -- its internal type order differs: there, type 6
# is H and type 7 is O; here it is the reverse).
#   2  14.0067  opls_902  N
#   3  12.0110  opls_145  C
#   4  14.0067  opls_901  N
#   5  12.0110  opls_157  C
#   6  15.9994  opls_154  O
#   7   1.0080  opls_156  H
#   8   1.0080  opls_910  H
#   9   1.0080  opls_155  H
ELEMENT_BY_TYPE_NETWORK: dict[int, str] = {
    2: "N",
    3: "C",
    4: "N",
    5: "C",
    6: "O",
    7: "H",
    8: "H",
    9: "H",
}

AT_TYPES_NETWORK = frozenset(ELEMENT_BY_TYPE_NETWORK)

# ref_elem, sel_elem, label, ref_file, peak_r_min_nm (ignored for C-C; uses windows)
PAIR_REFS: list[tuple[str, str, str, str, float]] = [
    ("C", "C", "C-C", "ref_C_C.xvg", R_MIN_PEAK_NM),
    ("C", "N", "C-N", "ref_C_N.xvg", R_MIN_PEAK_NM),
    ("O", "H", "O-H", "ref_O_H.xvg", 0.45),  # paper H-bond shell ~0.84 nm
]


@dataclass
class PeakWindow:
    label: str
    r_min_nm: float
    r_max_nm: float


CC_BONDED = PeakWindow("bonded", *CC_BONDED_R_NM)
CC_MELT = PeakWindow("melt shell", *CC_MELT_R_NM)


@dataclass
class RDFCurve:
    r_nm: np.ndarray
    g: np.ndarray
    label: str


@dataclass
class DumpFrame:
    step: int
    box_ang: np.ndarray
    ids: np.ndarray
    mol: np.ndarray
    types: np.ndarray
    xyz_ang: np.ndarray


def elements_for_frame(frame: DumpFrame, *, type_map: dict[int, str]) -> np.ndarray:
    """Element groups C*/N*/O*/H* via LAMMPS AT type (matches topol_aa atom types)."""
    return np.array([type_map.get(int(t), "") for t in frame.types], dtype=object)


def load_xvg(path: Path) -> RDFCurve:
    rows: list[tuple[float, float]] = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith(("#", "@")):
            continue
        parts = line.split()
        rows.append((float(parts[0]), float(parts[1])))
    data = np.array(rows, dtype=float)
    return RDFCurve(r_nm=data[:, 0], g=data[:, 1], label=path.name)


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
        header = lines[i].strip().split()[2:]
        i += 1
        try:
            id_idx = header.index("id")
            mol_idx = header.index("mol")
            type_idx = header.index("type")
            xyz_idx = [header.index("x"), header.index("y"), header.index("z")]
        except ValueError as exc:
            raise ValueError(f"dump missing id/mol/type/x/y/z columns: {header}") from exc
        ids = np.empty(n_atoms, dtype=int)
        mol = np.empty(n_atoms, dtype=int)
        types = np.empty(n_atoms, dtype=int)
        xyz = np.empty((n_atoms, 3), dtype=float)
        for j in range(n_atoms):
            parts = lines[i + j].split()
            ids[j] = int(parts[id_idx])
            mol[j] = int(parts[mol_idx])
            types[j] = int(parts[type_idx])
            xyz[j, 0] = float(parts[xyz_idx[0]])
            xyz[j, 1] = float(parts[xyz_idx[1]])
            xyz[j, 2] = float(parts[xyz_idx[2]])
        i += n_atoms
        if step >= min_step:
            frames.append(
                DumpFrame(step=step, box_ang=box, ids=ids, mol=mol, types=types, xyz_ang=xyz)
            )
    return frames


def minimum_image(delta: np.ndarray, box: np.ndarray) -> np.ndarray:
    return delta - box * np.round(delta / box)


def _accumulate_distances(
    hist: np.ndarray,
    dist: np.ndarray,
    *,
    dr: float,
    n_bins: int,
    r_max_ang: float,
) -> None:
    dist = dist.reshape(-1)
    dist = dist[np.isfinite(dist) & (dist > 1e-6) & (dist <= r_max_ang)]
    if dist.size == 0:
        return
    bins = (dist / dr).astype(int)
    bins = bins[(bins >= 0) & (bins < n_bins)]
    np.add.at(hist, bins, 1.0)


def compute_rdf(
    ref_xyz: np.ndarray,
    sel_xyz: np.ndarray,
    box_ang: np.ndarray,
    *,
    ref_mol: np.ndarray | None = None,
    sel_mol: np.ndarray | None = None,
    same_group: bool = False,
    r_max_nm: float = R_MAX_NM,
    n_bins: int = N_BINS,
    chunk: int = 64,
) -> RDFCurve:
    """Compute g(r). When mol arrays are set, count intermolecular pairs only.

    For same_group (C*-C*), uses GROMACS-style **ordered** pairs (each atom as
    reference, all j≠i), not unique i<j counting.
    """
    r_max_ang = r_max_nm * 10.0
    dr = r_max_ang / n_bins
    hist = np.zeros(n_bins, dtype=float)
    volume = float(np.prod(box_ang))
    n_sel = len(sel_xyz)
    rho = n_sel / volume
    inter_only = ref_mol is not None and sel_mol is not None

    n_ref = len(ref_xyz)
    if same_group:
        xyz = ref_xyz
        mol = ref_mol
        for i_start in range(0, n_ref, chunk):
            i_end = min(i_start + chunk, n_ref)
            for j_start in range(0, n_ref, chunk):
                j_end = min(j_start + chunk, n_ref)
                delta = xyz[i_start:i_end, None, :] - xyz[j_start:j_end][None, :, :]
                delta = minimum_image(delta, box_ang)
                dist = np.linalg.norm(delta, axis=2)
                ii = np.arange(i_start, i_end)[:, None]
                jj = np.arange(j_start, j_end)[None, :]
                mask = ii != jj
                if inter_only and mol is not None:
                    mask &= mol[i_start:i_end, None] != mol[j_start:j_end][None, :]
                dist = np.where(mask, dist, np.nan)
                _accumulate_distances(hist, dist, dr=dr, n_bins=n_bins, r_max_ang=r_max_ang)
    else:
        for start in range(0, n_ref, chunk):
            ref_chunk = ref_xyz[start : start + chunk]
            delta = ref_chunk[:, None, :] - sel_xyz[None, :, :]
            delta = minimum_image(delta, box_ang)
            dist = np.linalg.norm(delta, axis=2)
            if inter_only:
                dist = np.where(
                    ref_mol[start : start + chunk, None] != sel_mol[None, :],
                    dist,
                    np.nan,
                )
            _accumulate_distances(hist, dist, dr=dr, n_bins=n_bins, r_max_ang=r_max_ang)

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


def peak_in_window(curve: RDFCurve, window: PeakWindow) -> tuple[float, float]:
    mask = (curve.r_nm >= window.r_min_nm) & (curve.r_nm <= window.r_max_nm)
    if not np.any(mask):
        return 0.0, 0.0
    g = curve.g[mask]
    r = curve.r_nm[mask]
    idx = int(np.argmax(g))
    return float(r[idx]), float(g[idx])


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


def _compare_peak_metrics(
    name: str,
    computed: RDFCurve,
    reference: RDFCurve,
    *,
    pos_c: float,
    ht_c: float,
    pos_r: float,
    ht_r: float,
) -> tuple[bool, list[str]]:
    lines: list[str] = []
    ok_all = True
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
    return ok_all, lines


def compare_cc_pair(
    computed: RDFCurve,
    reference: RDFCurve,
) -> tuple[bool, list[str], int, int]:
    """C-C: separate bonded and melt-shell peaks plus long-range diagnostic."""
    lines: list[str] = []
    ok_all = True
    n_pass = 0
    n_fail = 0

    for window in (CC_BONDED, CC_MELT):
        pos_c, ht_c = peak_in_window(computed, window)
        pos_r, ht_r = peak_in_window(reference, window)
        lines.append(f"  [{window.label} r={window.r_min_nm:.2f}-{window.r_max_nm:.2f} nm]")
        ok, wlines = _compare_peak_metrics(
            f"C-C {window.label}",
            computed,
            reference,
            pos_c=pos_c,
            ht_c=ht_c,
            pos_r=pos_r,
            ht_r=ht_r,
        )
        ok_all &= ok
        lines.extend(wlines)
        n_pass += sum(1 for line in wlines if line.startswith("  PASS"))
        n_fail += sum(1 for line in wlines if line.startswith("  FAIL"))

    g_c = float(np.interp(LONG_RANGE_R_NM, computed.r_nm, computed.g))
    g_r = float(np.interp(LONG_RANGE_R_NM, reference.r_nm, reference.g))
    lines.append(
        f"  INFO: C-C g(r={LONG_RANGE_R_NM:.1f} nm) computed={g_c:.3f} ref={g_r:.3f} "
        "(GROMACS C*-C* self-RDF baseline ≠ 1)"
    )

    l2_val = l2(computed, reference)
    ok = l2_val <= L2_TOL
    ok_all &= ok
    line = f"  {'PASS' if ok else 'FAIL'}: C-C L2(g)={l2_val:.3f} (tol={L2_TOL})"
    lines.append(line)
    n_pass += int(ok)
    n_fail += int(not ok)
    return ok_all, lines, n_pass, n_fail


def compare_pair(
    name: str,
    computed: RDFCurve,
    reference: RDFCurve,
    *,
    peak_r_min_nm: float = R_MIN_PEAK_NM,
) -> tuple[bool, list[str]]:
    lines: list[str] = []
    pos_c, ht_c = first_peak(computed, peak_r_min_nm)
    pos_r, ht_r = first_peak(reference, peak_r_min_nm)
    _, lines = _compare_peak_metrics(
        name, computed, reference, pos_c=pos_c, ht_c=ht_c, pos_r=pos_r, ht_r=ht_r
    )
    ok_all = all("PASS" in line for line in lines if "peak" in line)
    l2_val = l2(computed, reference)
    ok = l2_val <= L2_TOL
    ok_all &= ok
    lines.append(f"  {'PASS' if ok else 'FAIL'}: {name} L2(g)={l2_val:.3f} (tol={L2_TOL})")
    return ok_all, lines


def parse_final_data_checks(
    data_path: Path,
    *,
    at_types: frozenset[int] | None = None,
) -> dict[str, float]:
    charges: dict[int, float] = {}
    masses: dict[int, float] = {}
    type_by_id: dict[int, int] = {}
    atoms: dict[int, tuple[float, float, float, int, int, int]] = {}
    bonds: list[tuple[int, int]] = []
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
        if stripped.startswith(("Velocities", "Angles", "Dihedrals", "Impropers")):
            section = None
            continue
        if "xlo" in stripped:
            box[0] = float(stripped.split()[1]) - float(stripped.split()[0])
        elif "ylo" in stripped:
            box[1] = float(stripped.split()[1]) - float(stripped.split()[0])
        elif "zlo" in stripped:
            box[2] = float(stripped.split()[1]) - float(stripped.split()[0])
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
            charges[aid] = float(parts[3])
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

    def atom_xyz(aid: int) -> np.ndarray:
        x, y, z, ix, iy, iz = atoms[aid]
        return np.array([x + ix * box[0], y + iy * box[1], z + iz * box[2]], dtype=float)

    max_bond = 0.0
    for a1, a2 in bonds:
        if at_types is not None and (
            type_by_id.get(a1) not in at_types or type_by_id.get(a2) not in at_types
        ):
            continue
        d = np.linalg.norm(minimum_image(atom_xyz(a2) - atom_xyz(a1), box))
        max_bond = max(max_bond, float(d))

    atom_ids = [aid for aid in atoms if at_types is None or type_by_id.get(aid) in at_types]
    total_mass = sum(masses[type_by_id[aid]] for aid in atom_ids)
    total_charge = sum(charges[aid] for aid in atom_ids)
    volume_cm3 = float(np.prod(box)) * 1e-24
    density = (total_mass / 6.022e23) / volume_cm3 if volume_cm3 > 0 else 0.0
    return {
        "total_charge": total_charge,
        "max_bond_ang": max_bond,
        "density_g_cm3": density,
        "volume_nm3": float(np.prod(box)) / 1000.0,
    }


def plot_rdgs(
    curves: list[tuple[RDFCurve, RDFCurve, str]],
    output: Path,
    *,
    cc_windows: tuple[PeakWindow, PeakWindow] | None = None,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = len(curves)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4), squeeze=False)
    for ax, (calc, ref, title) in zip(axes[0], curves, strict=True):
        ax.plot(calc.r_nm, calc.g, label="backmapped", linewidth=2.4)
        ax.plot(ref.r_nm, ref.g, "--", label="reference", linewidth=2.4)
        if cc_windows is not None and title == "C-C":
            for window, color in zip(cc_windows, ("#cccccc", "#e8e8ff"), strict=False):
                ax.axvspan(
                    window.r_min_nm, window.r_max_nm, alpha=0.25, color=color, label=window.label
                )
        ax.set_xlim(0, R_MAX_NM)
        ax.set_xlabel("r (nm)", fontsize=14)
        ax.set_ylabel("g(r)", fontsize=14)
        ax.set_title(title, fontsize=15)
        ax.tick_params(axis="both", labelsize=12)
        ax.legend(fontsize=12)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare melamine_network LAMMPS RDF vs bakery's crosslinked reference"
    )
    parser.add_argument("--dump", type=Path, required=True)
    parser.add_argument("--final-data", type=Path, required=True)
    parser.add_argument("--ref-dir", type=Path, default=Path("ref"))
    parser.add_argument("--min-step", type=int, default=12000)
    parser.add_argument(
        "--stride",
        type=int,
        default=1,
        help="Use every Nth frame (after --min-step) instead of all of them. "
        "The per-frame pairwise-distance computation is O(n^2) and doesn't "
        "scale comfortably past ~1000 frames x 11475 atoms (8h25m for the "
        "full 996-frame, 1ns production run) -- use a stride for a faster, "
        "visually-representative rerun (e.g. figure restyling) when the "
        "full-precision numbers aren't being recomputed.",
    )
    parser.add_argument(
        "--inter",
        action="store_true",
        help="Intermolecular pairs only (default: all pairs, GROMACS g_rdf style)",
    )
    parser.add_argument("--plot", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=Path("structural_validation_report.txt"))
    args = parser.parse_args()

    type_map = ELEMENT_BY_TYPE_NETWORK
    at_types = AT_TYPES_NETWORK

    for path in (args.dump, args.final_data, args.ref_dir):
        if not path.exists():
            print(f"ERROR: missing {path}", file=sys.stderr)
            return 1

    pair_mode = "intermolecular only" if args.inter else "all pairs (GROMACS g_rdf)"
    print("=== melamine_network structural validation vs bakery crosslinked reference ===\n")
    print("Element groups: C*/N*/O*/H* via LAMMPS AT types (this example's own Masses map)")
    print(f"RDF pair selection: {pair_mode}\n")

    frames = parse_lammps_dump(args.dump, min_step=args.min_step)
    if not frames:
        print("ERROR: no dump frames at or after min-step", file=sys.stderr)
        return 1
    if args.stride > 1:
        frames = frames[:: args.stride]
    print(
        f"Dump: {len(frames)} frames from step {frames[0].step} to {frames[-1].step}"
        + (f" (stride={args.stride})" if args.stride > 1 else "")
    )

    computed_by_pair: dict[str, list[RDFCurve]] = {label: [] for _, _, label, _, _ in PAIR_REFS}
    for frame in frames:
        elements = elements_for_frame(frame, type_map=type_map)
        mol_arg = frame.mol if args.inter else None
        for ref_elem, sel_elem, label, _, _ in PAIR_REFS:
            ref_mask = elements == ref_elem
            sel_mask = elements == sel_elem
            if not ref_mask.any() or not sel_mask.any():
                continue
            ref_xyz = frame.xyz_ang[ref_mask]
            sel_xyz = frame.xyz_ang[sel_mask]
            same_group = ref_elem == sel_elem
            ref_mol = frame.mol[ref_mask] if mol_arg is not None else None
            sel_mol = frame.mol[sel_mask] if mol_arg is not None else None
            computed_by_pair[label].append(
                compute_rdf(
                    ref_xyz,
                    sel_xyz,
                    frame.box_ang,
                    ref_mol=ref_mol,
                    sel_mol=sel_mol if not same_group else ref_mol,
                    same_group=same_group,
                )
            )

    n_pass = 0
    n_fail = 0
    plot_pairs: list[tuple[RDFCurve, RDFCurve, str]] = []
    report_lines = [
        "melamine_network Tier C structural validation vs bakery crosslinked reference",
        f"RDF pair selection: {pair_mode}",
        "",
    ]

    for _, _, label, ref_name, peak_r_min in PAIR_REFS:
        ref_path = args.ref_dir / ref_name
        if not ref_path.exists():
            print(f"ERROR: missing reference {ref_path}", file=sys.stderr)
            return 1
        if not computed_by_pair[label]:
            print(f"\n--- {label}: no atoms matched ---")
            report_lines.append(f"{label}: FAIL (no atoms matched)")
            n_fail += 5 if label == "C-C" else 3
            continue
        computed = average_rdf(computed_by_pair[label])
        reference = load_xvg(ref_path)
        print(f"\n--- {label} vs {ref_name} ---")
        report_lines.append(f"--- {label} vs {ref_name} ---")
        if label == "C-C":
            _, lines, p_pass, p_fail = compare_cc_pair(computed, reference)
            n_pass += p_pass
            n_fail += p_fail
        else:
            _, lines = compare_pair(label, computed, reference, peak_r_min_nm=peak_r_min)
            n_pass += sum(1 for line in lines if line.startswith("  PASS"))
            n_fail += sum(1 for line in lines if line.startswith("  FAIL"))
        for line in lines:
            print(line)
            report_lines.append(line.strip())
        plot_pairs.append((computed, reference, label))

    checks = parse_final_data_checks(args.final_data, at_types=at_types)
    print("\n--- final structure checks (AT atoms only) ---")
    report_lines.extend(["", "--- final structure checks (AT atoms only) ---"])
    charge_ok = abs(checks["total_charge"]) < 1e-2
    bond_ok = checks["max_bond_ang"] < 20.0
    print(
        f"  {'PASS' if charge_ok else 'FAIL'}: total charge = {checks['total_charge']:.6f} e (tol 1e-2)"
    )
    print(
        f"  {'PASS' if bond_ok else 'FAIL'}: max min-image bond = {checks['max_bond_ang']:.2f} Å (tol 20 Å)"
    )
    print(f"  INFO: density ≈ {checks['density_g_cm3']:.3f} g/cm³")
    report_lines.append(
        f"total charge = {checks['total_charge']:.6f} e; max AT bond = {checks['max_bond_ang']:.2f} Å; "
        f"density ≈ {checks['density_g_cm3']:.3f} g/cm³"
    )

    summary = f"\nSummary: {n_pass} PASS, {n_fail} FAIL (RDF peak metrics)"
    print(summary)
    report_lines.extend(["", summary.strip()])
    args.report.write_text("\n".join(report_lines) + "\n")
    print(f"Report: {args.report}")

    if args.plot and plot_pairs:
        plot_rdgs(plot_pairs, args.plot, cc_windows=(CC_BONDED, CC_MELT))
        print(f"Plot: {args.plot}")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
