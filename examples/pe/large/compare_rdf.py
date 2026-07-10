# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "matplotlib"]
# ///
"""Compare RDF from backmapped system against a reference AT simulation.

Reads LAMMPS fix ave/time RDF output files and compares g(r) curves.
Reports peak position/height agreement and overall RMSD.

Usage:
  uv run compare_rdf.py --backmap rdf_backmap.dat --reference rdf_reference.dat [--plot rdf_comparison.png]

Note on type mapping:
  The backmapping RDF uses the hybrid type numbering (3=CH3, 4=CH2),
  while the reference uses renumbered types (1=CH3, 2=CH2).
  Both files contain the same 3 pair columns: CH3-CH3, CH3-CH2, CH2-CH2.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PAIR_LABELS = ["CH3-CH3", "CH3-CH2", "CH2-CH2"]

PEAK_POS_TOL = 0.3  # Å tolerance on first-peak position
PEAK_HT_TOL = 0.35  # relative tolerance on first-peak height
RMSD_TOL = 0.25  # absolute RMSD tolerance on g(r)


@dataclass
class RDFData:
    r: np.ndarray
    gr: list[np.ndarray]  # one g(r) array per pair


def parse_lammps_rdf(path: Path) -> RDFData:
    """Parse LAMMPS fix ave/time RDF output (mode vector).

    Format: comment lines starting with #, then repeating blocks of:
      timestep nrows     (2-column header)
      1 r g1 c1 g2 c2 g3 c3   (nrows data lines, 8 columns each)
    We take the last complete block (final time-average).
    """
    blocks: list[list[str]] = []
    current_block: list[str] = []
    nrows = 0
    rows_read = 0

    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()

        if len(parts) == 2 and rows_read == nrows:
            if current_block:
                blocks.append(current_block)
            current_block = []
            nrows = int(parts[1])
            rows_read = 0
            continue

        current_block.append(s)
        rows_read += 1

    if current_block:
        blocks.append(current_block)

    if not blocks:
        raise ValueError(f"No RDF data blocks found in {path}")

    last_block = blocks[-1]
    data = np.array([[float(x) for x in line.split()] for line in last_block])

    r = data[:, 1]
    n_pairs = (data.shape[1] - 2) // 2
    return RDFData(r=r, gr=[data[:, 2 + 2 * i] for i in range(n_pairs)])


def find_first_peak(r: np.ndarray, gr: np.ndarray, r_min: float = 2.0) -> tuple[float, float]:
    """Find position and height of the first peak above r_min."""
    mask = r > r_min
    if not np.any(mask):
        return 0.0, 0.0

    r_cut = r[mask]
    g_cut = gr[mask]

    peak_idx = np.argmax(g_cut)
    return float(r_cut[peak_idx]), float(g_cut[peak_idx])


def compute_rmsd(gr1: np.ndarray, gr2: np.ndarray) -> float:
    """RMSD between two g(r) curves (same r grid assumed)."""
    n = min(len(gr1), len(gr2))
    diff = gr1[:n] - gr2[:n]
    return float(np.sqrt(np.mean(diff**2)))


def interpolate_to_common_grid(
    rdf1: RDFData, rdf2: RDFData
) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
    """Interpolate both RDFs onto a common r grid."""
    r_min = max(rdf1.r[0], rdf2.r[0])
    r_max = min(rdf1.r[-1], rdf2.r[-1])
    n_points = min(len(rdf1.r), len(rdf2.r))
    r_common = np.linspace(r_min, r_max, n_points)

    gr1_interp = [np.interp(r_common, rdf1.r, g) for g in rdf1.gr]
    gr2_interp = [np.interp(r_common, rdf2.r, g) for g in rdf2.gr]

    return r_common, gr1_interp, gr2_interp


def compare_rdfs(
    backmap: RDFData,
    reference: RDFData,
    peak_pos_tol: float = PEAK_POS_TOL,
    peak_ht_tol: float = PEAK_HT_TOL,
    rmsd_tol: float = RMSD_TOL,
) -> tuple[int, int]:
    """Compare RDF curves and print results. Returns (n_pass, n_fail)."""
    n_pairs = min(len(backmap.gr), len(reference.gr))
    r_common, gr_bm, gr_ref = interpolate_to_common_grid(backmap, reference)

    n_pass = 0
    n_fail = 0

    for i in range(n_pairs):
        label = PAIR_LABELS[i] if i < len(PAIR_LABELS) else f"pair_{i + 1}"
        print(f"\n  --- {label} ---")

        pos_bm, ht_bm = find_first_peak(r_common, gr_bm[i])
        pos_ref, ht_ref = find_first_peak(r_common, gr_ref[i])

        pos_err = abs(pos_bm - pos_ref)
        ok = pos_err < peak_pos_tol
        status = "PASS" if ok else "FAIL"
        print(
            f"  {status}: 1st peak position: backmap={pos_bm:.2f} Å, ref={pos_ref:.2f} Å (Δ={pos_err:.2f}, tol={peak_pos_tol})"
        )
        n_pass += ok
        n_fail += not ok

        ht_rel = abs(ht_bm - ht_ref) / ht_ref if ht_ref > 0 else 0.0
        ok = ht_rel < peak_ht_tol
        status = "PASS" if ok else "FAIL"
        print(
            f"  {status}: 1st peak height: backmap={ht_bm:.3f}, ref={ht_ref:.3f} (rel err={ht_rel:.2f}, tol={peak_ht_tol})"
        )
        n_pass += ok
        n_fail += not ok

        rmsd = compute_rmsd(gr_bm[i], gr_ref[i])
        ok = rmsd < rmsd_tol
        status = "PASS" if ok else "FAIL"
        print(f"  {status}: RMSD(g(r)) = {rmsd:.4f} (tol={rmsd_tol})")
        n_pass += ok
        n_fail += not ok

    return n_pass, n_fail


def plot_comparison(
    backmap: RDFData,
    reference: RDFData,
    output: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n_pairs = min(len(backmap.gr), len(reference.gr))
    fig, axes = plt.subplots(1, n_pairs, figsize=(5 * n_pairs, 4), squeeze=False)

    for i in range(n_pairs):
        ax = axes[0, i]
        label = PAIR_LABELS[i] if i < len(PAIR_LABELS) else f"pair_{i + 1}"
        ax.plot(backmap.r, backmap.gr[i], label="backmapped", linewidth=1.5)
        ax.plot(reference.r, reference.gr[i], "--", label="reference AT", linewidth=1.5)
        ax.set_xlabel("r (Å)")
        ax.set_ylabel("g(r)")
        ax.set_title(label)
        ax.legend()
        ax.set_xlim(0, min(backmap.r[-1], reference.r[-1]))

    fig.tight_layout()
    fig.savefig(output, dpi=150)
    print(f"\n  Plot saved to {output}")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compare backmapped RDF against reference AT simulation"
    )
    parser.add_argument(
        "--backmap", type=Path, required=True, help="LAMMPS RDF output from backmapping run"
    )
    parser.add_argument(
        "--reference", type=Path, required=True, help="LAMMPS RDF output from standalone AT run"
    )
    parser.add_argument(
        "--plot",
        type=Path,
        default=None,
        help="Save comparison plot to file (e.g. rdf_comparison.png)",
    )
    parser.add_argument("--peak-pos-tol", type=float, default=PEAK_POS_TOL)
    parser.add_argument("--peak-ht-tol", type=float, default=PEAK_HT_TOL)
    parser.add_argument("--rmsd-tol", type=float, default=RMSD_TOL)
    args = parser.parse_args()

    for p in (args.backmap, args.reference):
        if not p.exists():
            print(f"ERROR: {p} not found", file=sys.stderr)
            return 1

    print("=== RDF Comparison: backmapped vs reference AT ===")

    backmap = parse_lammps_rdf(args.backmap)
    reference = parse_lammps_rdf(args.reference)

    print(f"\n  Backmapped: {len(backmap.r)} bins, {len(backmap.gr)} pairs")
    print(f"  Reference:  {len(reference.r)} bins, {len(reference.gr)} pairs")

    n_pass, n_fail = compare_rdfs(
        backmap,
        reference,
        peak_pos_tol=args.peak_pos_tol,
        peak_ht_tol=args.peak_ht_tol,
        rmsd_tol=args.rmsd_tol,
    )

    if args.plot:
        plot_comparison(backmap, reference, args.plot)

    print(f"\n{'=' * 50}")
    print(f"Results: {n_pass} passed, {n_fail} failed")

    if n_fail == 0:
        print("\nBackmapped RDF matches reference AT simulation.")
    else:
        print("\nRDF mismatch — possible causes:")
        print("  - Insufficient equilibration in Phase 3 (increase run length)")
        print("  - Too few molecules for good statistics (increase system size)")
        print("  - Force field mismatch between backmapped and reference")
        print("  - Lambda ramp too fast (decrease alpha)")

    return 1 if n_fail > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
