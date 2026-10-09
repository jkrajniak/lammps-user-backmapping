"""Tier C for the MARTINI POPC example: backmapped bilayer vs the Slipids reference.

For each trajectory (GROMACS tpr + xtc, Slipids POPC naming), over the frames
after ``--skip`` ps: area per lipid (Lx*Ly / lipids per leaflet), bilayer
thickness (mean z distance between the phosphorus planes of the two leaflets)
and the C-H order parameter profiles S_CH = <3 cos^2(theta) - 1> / 2 of the
sn-1 (palmitoyl, C32..C316) and sn-2 (oleoyl, C22..C218) chains, theta the angle
between C-H and the bilayer normal (z), averaged over the hydrogens of each
carbon and all lipids. Means with standard errors from ``--blocks`` blocks.

    uv run --with MDAnalysis --with matplotlib popc_tierc.py \\
        --traj backmapped prod.tpr prod.xtc --traj reference ref/prod.tpr ref/prod.xtc \\
        --skip 2000 --out tierc_popc [--plot]
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np

SN1 = [f"C3{i}" for i in range(2, 17)]
SN2 = [f"C2{i}" for i in range(2, 19)]


def hydrogens_of(carbon: str) -> str:
    """Selection of the hydrogens bonded to a chain carbon (Slipids/CHARMM names)."""
    m = re.fullmatch(r"C([23])(\d+)", carbon)
    chain, pos = m.group(1), m.group(2)
    if chain == "3":
        return " ".join(f"H{pos}{s}" for s in "XYZ")
    if pos in ("9", "10"):  # oleoyl double bond: one H each
        return f"H{pos}1"
    return " ".join(f"H{pos}{s}" for s in "RST")


def analyse(tpr: str, xtc: str, skip_ps: float, n_blocks: int) -> dict:
    from MDAnalysis import Universe

    u = Universe(tpr, xtc)
    lipids = u.select_atoms("resname POPC")
    n_lipid = lipids.n_residues
    p_atoms = u.select_atoms("resname POPC and name P")
    ends = u.select_atoms("resname POPC and name C218 C316")
    pairs = {}
    for chain, carbons in (("sn1", SN1), ("sn2", SN2)):
        for c in carbons:
            cs = u.select_atoms(f"resname POPC and name {c}")
            hs = u.select_atoms(f"resname POPC and name {hydrogens_of(c)}")
            # pair each hydrogen with the carbon of its residue
            c_by_res = {a.resindex: a.index for a in cs}
            pairs[(chain, c)] = (
                np.array([c_by_res[h.resindex] for h in hs]),
                hs.indices,
            )
    apl, thick, times = [], [], []
    s_ch = {k: [] for k in pairs}
    t0 = u.trajectory[0].time
    for ts in u.trajectory:
        if ts.time < t0 + skip_ps:
            continue
        lx, ly, lz = ts.dimensions[:3]
        times.append(ts.time)
        apl.append(lx * ly / (n_lipid / 2) / 100.0)  # nm^2
        # midplane: circular mean of the terminal methyl z (robust to PBC)
        ang = 2 * np.pi * ends.positions[:, 2] / lz
        mid = (np.arctan2(np.sin(ang).mean(), np.cos(ang).mean()) % (2 * np.pi)) * lz / (2 * np.pi)
        dz = p_atoms.positions[:, 2] - mid
        dz -= lz * np.round(dz / lz)
        thick.append((dz[dz > 0].mean() - dz[dz < 0].mean()) / 10.0)  # nm
        x = u.atoms.positions
        for key, (ci, hi) in pairs.items():
            v = x[hi] - x[ci]
            v -= ts.dimensions[:3] * np.round(v / ts.dimensions[:3])
            cos2 = (v[:, 2] / np.linalg.norm(v, axis=1)) ** 2
            s_ch[key].append(float(np.mean(1.5 * cos2 - 0.5)))

    def mean_sem(values: list[float]) -> tuple[float, float]:
        arr = np.asarray(values)
        blocks = np.array([b.mean() for b in np.array_split(arr, n_blocks)])
        return float(arr.mean()), float(blocks.std(ddof=1) / np.sqrt(n_blocks))

    return {
        "frames": len(times),
        "t_range": (times[0], times[-1]) if times else (0, 0),
        "apl": mean_sem(apl),
        "thickness": mean_sem(thick),
        "s_ch": {k: mean_sem(v) for k, v in s_ch.items()},
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument(
        "--traj", nargs=3, action="append", metavar=("LABEL", "TPR", "XTC"), required=True
    )
    ap.add_argument(
        "--skip", type=float, default=2000.0, help="ps discarded at the start of each trajectory"
    )
    ap.add_argument("--blocks", type=int, default=5)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    res = {label: analyse(tpr, xtc, args.skip, args.blocks) for label, tpr, xtc in args.traj}
    labels = list(res)
    lines = [f"{'quantity':24s}" + "".join(f"{lab:>22s}" for lab in labels)]
    lines.append(
        f"{'frames (t, ps)':24s}"
        + "".join(
            f"{res[lab]['frames']:>8d} ({res[lab]['t_range'][0]:.0f}-{res[lab]['t_range'][1]:.0f})".rjust(
                22
            )
            for lab in labels
        )
    )
    for q, unit in (("apl", "nm^2"), ("thickness", "nm")):
        lines.append(
            f"{q + ' / ' + unit:24s}"
            + "".join(f"{res[lab][q][0]:12.4f} +- {res[lab][q][1]:.4f}".rjust(22) for lab in labels)
        )
    for key in res[labels[0]]["s_ch"]:
        lines.append(
            f"{'S_CH ' + key[0] + ' ' + key[1]:24s}"
            + "".join(
                f"{res[lab]['s_ch'][key][0]:12.4f} +- {res[lab]['s_ch'][key][1]:.4f}".rjust(22)
                for lab in labels
            )
        )
    text = "\n".join(lines)
    Path(f"{args.out}.txt").write_text(text + "\n")
    print(text)

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(6.8, 2.8), sharey=True)
        for ax, chain, carbons in ((axes[0], "sn1", SN1), (axes[1], "sn2", SN2)):
            pos = [int(c[2:]) for c in carbons]
            for k, lab in enumerate(labels):
                m = [-res[lab]["s_ch"][(chain, c)][0] for c in carbons]
                e = [res[lab]["s_ch"][(chain, c)][1] for c in carbons]
                ax.errorbar(
                    pos,
                    m,
                    yerr=e,
                    marker="o",
                    ms=3,
                    lw=1,
                    capsize=2,
                    ls="-" if k == 0 else "--",
                    label=lab,
                )
            ax.set_xlabel(f"carbon ({chain})")
        axes[0].set_ylabel("$-S_{CH}$")
        axes[0].legend(frameon=False, fontsize=7)
        fig.tight_layout()
        fig.savefig(f"{args.out}_scd.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
