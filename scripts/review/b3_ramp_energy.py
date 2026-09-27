"""Review item B3 (R1-1): temperature and energy terms through the backmapping.

Walks a hybrid-run LAMMPS log (the generated protocol) in order, keeps every
thermo row with its stage (minimize, relax, ramp, nvt...), the timestep in force
and the elapsed simulated time, and writes a CSV and a stage summary (max AT
temperature and energies at the end of each stage). With AT production logs of
the backmapped system and of the independent reference, it also compares their
potential energy per atom (block mean +- SEM over the production thermo rows).

    uv run --with numpy --with matplotlib b3_ramp_energy.py --hybrid log.pe.lammps \\
        --backmap-at log.pe_at.lammps --reference-at log.pe_at_ref.lammps \\
        --prod-steps 5000000 --out b3_pe [--plot]
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np


def walk(log: Path) -> list[dict]:
    rows: list[dict] = []
    dt, time, stage, cols, last_step = 1.0, 0.0, "setup", None, None
    n_runs = 0
    for line in log.read_text().splitlines():
        s = line.strip()
        if s.startswith("timestep "):
            dt = float(s.split()[1])
        elif s.startswith("minimize "):
            stage, last_step = "minimize", None
        elif re.match(r"^run \d+", s):
            n_runs += 1
            stage, last_step = f"run{n_runs}", None
        elif s.startswith("Step "):
            cols = s.split()
        elif s.startswith("Loop time"):
            cols = None
        elif cols and s and s.split()[0].isdigit() and len(s.split()) == len(cols):
            vals = [float(v) for v in s.split()]
            row = dict(zip(cols, vals, strict=True))
            step = int(row["Step"])
            if stage != "minimize" and last_step is not None:
                time += (step - last_step) * dt
            last_step = step
            row.update(stage=stage, dt_fs=dt, time_fs=time)
            rows.append(row)
    return rows


def label_stages(rows: list[dict]) -> None:
    """Name run stages from the timestep and lambda (generated robust protocol)."""
    for r in rows:
        if r["stage"] == "minimize":
            continue
        lam = r.get("lambda", r.get("f_bm", float("nan")))
        dt = r["dt_fs"]
        if dt <= 0.011:
            r["stage"] = "relax"
        elif abs(dt - 0.1) < 1e-9:
            r["stage"] = "ramp" if lam < 1.0 else "ramp_hold"
        else:
            r["stage"] = f"nvt_dt{dt:g}"


def production_pe(log: Path, prod_steps: int, n_blocks: int = 5) -> tuple[float, float, int]:
    """PE per atom over the last prod_steps of an AT log: mean and block SEM."""
    rows = walk(log)
    atoms = int(re.search(r"^\s*(\d+) atoms$", log.read_text(), re.MULTILINE).group(1))
    last = rows[-1]["Step"]
    pe = np.array([r["PotEng"] for r in rows if r["Step"] > last - prod_steps]) / atoms
    blocks = np.array_split(pe, n_blocks)
    means = np.array([b.mean() for b in blocks])
    return float(pe.mean()), float(means.std(ddof=1) / np.sqrt(n_blocks)), atoms


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--hybrid", type=Path, required=True)
    ap.add_argument("--backmap-at", type=Path)
    ap.add_argument("--reference-at", type=Path)
    ap.add_argument("--prod-steps", type=int, default=5000000)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--plot", action="store_true")
    args = ap.parse_args()

    rows = walk(args.hybrid)
    label_stages(rows)
    keys = sorted({k for r in rows for k in r if k not in ("stage",)})
    with open(f"{args.out}.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["stage", *keys])
        w.writeheader()
        w.writerows(rows)

    report = [f"hybrid log {args.hybrid}: {len(rows)} thermo rows"]
    report.append(
        f"{'stage':12s} {'rows':>5s} {'t_end/ps':>9s} {'lambda_end':>10s} {'Tmax/K':>9s} {'T_end/K':>8s} {'PE_end':>12s}"
    )
    for stage in dict.fromkeys(r["stage"] for r in rows):
        sub = [r for r in rows if r["stage"] == stage]
        lam = sub[-1].get("lambda", sub[-1].get("f_bm", float("nan")))
        report.append(
            f"{stage:12s} {len(sub):5d} {sub[-1]['time_fs'] / 1000:9.3f} {lam:10.4f} "
            f"{max(r['Temp'] for r in sub):9.0f} {sub[-1]['Temp']:8.0f} {sub[-1]['PotEng']:12.1f}"
        )
    if args.backmap_at and args.reference_at:
        pb, sb, nb = production_pe(args.backmap_at, args.prod_steps)
        pr, sr, nr = production_pe(args.reference_at, args.prod_steps)
        report.append(
            f"AT production PE per atom: backmapped {pb:.5f} +- {sb:.5f} ({nb} atoms), "
            f"reference {pr:.5f} +- {sr:.5f} ({nr} atoms), difference {pb - pr:+.5f} kcal/mol"
        )
    text = "\n".join(report)
    Path(f"{args.out}.txt").write_text(text + "\n")
    print(text)

    if args.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        dyn = [r for r in rows if r["stage"] != "minimize"]
        t = np.array([r["time_fs"] for r in dyn]) / 1000
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(3.4, 4.0), sharex=True)
        a1.plot(t, [r["Temp"] for r in dyn], color="C3", lw=1)
        a1.set_yscale("log")
        a1.set_ylabel("AT temperature / K")
        for term, color in (
            ("PotEng", "k"),
            ("E_bond", "C0"),
            ("E_angle", "C1"),
            ("E_dihed", "C2"),
            ("E_vdwl", "C4"),
        ):
            if term in dyn[0]:
                a2.plot(t, [r[term] for r in dyn], color=color, lw=1, label=term)
        a2.set_ylabel("energy / kcal mol$^{-1}$")
        a2.set_xlabel("time / ps")
        a2.legend(fontsize=6, frameon=False)
        fig.tight_layout()
        fig.savefig(f"{args.out}.pdf")
    return 0


if __name__ == "__main__":
    sys.exit(main())
