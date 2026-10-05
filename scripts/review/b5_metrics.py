"""Review item B5 (R1-2): metrics of the ramp-time sweep, one row per job directory of b5_ramp_sweep.sh.

For every ``<job>/b5run`` (alpha, seed from b5.info) it reports

- rc: exit status of the Tier B/C script (a failed run is a result);
- the ramp: temperature maximum of the atomistic sites, temperature at the end of the ramp, maximum of the
  potential energy, and the time after the ramp until the temperature stays within 5 % of the target (b3 walker);
- the structure of the atomistic production against the independent reference: RDF L2 norms (mean of the pairs),
  CH2-CH2 first-peak height difference, and the Jensen-Shannon divergence of the bond-angle and dihedral
  distributions and the trans fraction of the dihedrals (b6 functions);
- ramp length in steps and ps (from the log).

    uv run --with numpy b5_metrics.py --reference <dir with rdf_reference.dat angle_hist_reference.dat dihedral_hist_reference.dat> \\
        --target 298 --csv b5_dodecane.csv <jobdir> [<jobdir> ...]

The averages over the seeds of one alpha, with their spread, are printed at the end.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import re
import statistics
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent


def load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


b3 = load("b3_ramp_energy", HERE / "b3_ramp_energy.py")
b6 = load("b6_bonded_distributions", HERE / "b6_bonded_distributions.py")


def info(job: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    f = job / "b5.info"
    if not f.exists():
        return out
    for ln in f.read_text().splitlines():
        m = re.match(r"alpha=(\S+) seed=(\S+)", ln)
        if m:
            out["alpha"], out["seed"] = m.group(1), m.group(2)
        m = re.match(r"rc=(\d+)", ln)
        if m:
            out["rc"] = m.group(1)
    return out


def ramp_metrics(log: Path, target: float) -> dict[str, float]:
    rows = b3.walk(log)
    b3.label_stages(rows)
    ramp = [r for r in rows if r.get("stage") == "ramp"]
    if not ramp:
        return {}
    t_ramp_end = max(r["time_fs"] for r in ramp) / 1000.0
    t_start = min(r["time_fs"] for r in ramp) / 1000.0
    after = [r for r in rows if r["time_fs"] / 1000.0 >= t_ramp_end and r.get("Temp", 0) > 0]
    settle = float("nan")
    for i, r in enumerate(after):
        if all(abs(x["Temp"] - target) <= 0.05 * target for x in after[i:]):
            settle = r["time_fs"] / 1000.0 - t_ramp_end
            break
    return {
        "ramp_ps": t_ramp_end - t_start,
        "Tmax_K": max(r["Temp"] for r in ramp),
        "T_end_ramp_K": ramp[-1]["Temp"],
        "PEmax": max(r["PotEng"] for r in ramp),
        "settle_ps": settle,
    }


def rdf_metrics(bm: Path, ref: Path, skip_first: int = 1) -> dict[str, float]:
    cmp = load(
        "compare_rdf_blocks",
        HERE.parent.parent / "examples" / "dodecane" / "large" / "compare_rdf_blocks.py",
    )
    a, b = cmp.parse_blocks(bm), cmp.parse_blocks(ref)
    l2s, dh = [], []
    for ga, gb in zip(a.gr_blocks, b.gr_blocks, strict=True):
        ma, _ = cmp.mean_sem(ga[skip_first:])
        mb, _ = cmp.mean_sem(gb[skip_first:])
        l2s.append(cmp.l2(ma, mb))
        pa = cmp.first_peak(a.r, ma)
        pb = cmp.first_peak(b.r, mb)
        dh.append(abs(pa[1] - pb[1]) / pb[1])
    return {
        "rdf_L2_mean": float(np.mean(l2s)),
        "rdf_L2_max": float(np.max(l2s)),
        "peak_height_relerr_max": float(np.max(dh)),
    }


def hist_metrics(bm: Path, ref: Path, dihedral: bool, skip_first: int = 1) -> dict[str, float]:
    ca, ba = b6.parse_blocks(bm)
    _, bb = b6.parse_blocks(ref)
    pa = b6.normalize(ba[skip_first:].sum(axis=0))
    pb = b6.normalize(bb[skip_first:].sum(axis=0))
    key = "dih" if dihedral else "ang"
    out = {f"{key}_JS": b6.js_divergence(pa, pb)}
    if dihedral:
        trans = np.abs(ca) > 120.0
        out["trans_bm"] = float(pa[trans].sum())
        out["trans_ref"] = float(pb[trans].sum())
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("jobs", nargs="+", type=Path)
    ap.add_argument("--reference", type=Path, required=True)
    ap.add_argument("--target", type=float, default=298.0)
    ap.add_argument("--hybrid-log", default="log.dodecane.lammps")
    ap.add_argument("--csv", type=Path)
    args = ap.parse_args()

    rows: list[dict] = []
    for job in args.jobs:
        run = job / "b5run" if (job / "b5run").exists() else job
        row: dict = {"job": job.name, **info(run)}
        try:
            row.update(ramp_metrics(run / args.hybrid_log, args.target))
            if row.get("rc") in (None, "0"):
                row.update(
                    rdf_metrics(run / "rdf_backmap.dat", args.reference / "rdf_reference.dat")
                )
                row.update(
                    hist_metrics(
                        run / "angle_hist_backmap.dat",
                        args.reference / "angle_hist_reference.dat",
                        False,
                    )
                )
                row.update(
                    hist_metrics(
                        run / "dihedral_hist_backmap.dat",
                        args.reference / "dihedral_hist_reference.dat",
                        True,
                    )
                )
        except (OSError, ValueError, KeyError) as exc:
            row["error"] = str(exc)[:80]
        rows.append(row)

    fields = sorted(
        {k for r in rows for k in r}, key=lambda k: (k not in ("job", "alpha", "seed", "rc"), k)
    )
    if args.csv:
        with args.csv.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            w.writerows(rows)
    by_alpha: dict[str, list[dict]] = {}
    for r in rows:
        by_alpha.setdefault(r.get("alpha", "?"), []).append(r)
    cols = [
        "ramp_ps",
        "Tmax_K",
        "settle_ps",
        "rdf_L2_mean",
        "peak_height_relerr_max",
        "ang_JS",
        "dih_JS",
        "trans_bm",
    ]
    print("alpha      n  ok  " + "  ".join(f"{c:>22s}" for c in cols))
    for alpha in sorted(by_alpha, key=float):
        g = by_alpha[alpha]
        ok = sum(1 for r in g if r.get("rc") == "0")
        cells = []
        for c in cols:
            v = [r[c] for r in g if isinstance(r.get(c), float) and np.isfinite(r[c])]
            cells.append(
                f"{statistics.mean(v):10.4g} +-{statistics.pstdev(v):8.2g}"
                if v
                else f"{'n/a':>22s}"
            )
        print(f"{alpha:8s} {len(g):3d} {ok:3d}  " + "  ".join(cells))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
