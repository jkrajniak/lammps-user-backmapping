#!/usr/bin/env python3
"""Direct ring C-N-C angle/distance variance check, bypassing
compare_melamine_structure.py's RDF binning and box-density normalization
entirely, to discriminate between two competing hypotheses for why the
melamine C-C(ring 1-3) and C-N RDF peaks come out ~30% too tall at the
correct position:

  (A) genuine dynamical under-broadening (Langevin damping, fix bm
      restraint, or charge-asymmetry-driven rigidity actually narrows the
      angle distribution), vs.
  (B) an RDF-normalization artifact (box was never barostatted/relaxed from
      the CG mapping's density, inflating density-independent intramolecular
      peaks via the 1/rho term).

Method: read the ring C-N-C angle atom triplets (angle types 37/38/41,
K=40.000044 kcal/mol/rad^2 for all three, theta0 differs slightly per
triplet -- confirms the force field's per-atom charges/angles were never
symmetrized across the molecule's 3-fold-symmetric arms) directly from the
final data file's topology, then measure the *raw* angle distribution
straight from dump.at_prod coordinates (minimum-image corrected; the dump
uses wrapped x/y/z, not unwrapped) for every ring instance across every
production frame (step >= 70000, matching the validated RDF comparison
window exactly). Compare the measured variance to the equipartition
prediction Var(theta) = kT/(2K) -- this measurement never touches the RDF
script's binning or box-density normalization, so it cleanly isolates
hypothesis (A) from (B):

  - measured variance well below kT/(2K)  -> real dynamical narrowing (A)
  - measured variance matches kT/(2K)     -> RDF peak inflation must be a
                                              normalization artifact (B)

Also reports the raw (non-density-normalized) C-C ring-flanking distance
distribution directly, as a secondary cross-check on the same question.

Usage (run on the VM, in the melamine_bakery dir where dump.at_prod and
melamine_hybrid_final.data live):
  .venv_rdf/bin/python ring_angle_variance_check.py
"""

from __future__ import annotations

import argparse
import math
import sys

import numpy as np

DATA_FILE = "melamine_hybrid_final.data"
KB_KCAL_PER_MOL_K = 0.0019872041
TEMP_K = 300.0

# angle_coeff values, from in.melamine_bakery_faithful.lammps (verified
# against the GROMACS source k via units.py's SPRING_ANGLE conversion):
#   type  K (kcal/mol/rad^2)  theta0 (deg)
RING_ANGLE_COEFFS = {
    37: (40.000044, 114.8430),
    38: (40.000044, 114.8070),
    41: (40.000044, 115.1170),
}


def parse_ring_angles(data_path: str) -> dict[int, list[tuple[int, int, int]]]:
    """type -> list of (atom_i, atom_j[vertex], atom_k) for ring C-N-C angles."""
    triples: dict[int, list[tuple[int, int, int]]] = {t: [] for t in RING_ANGLE_COEFFS}
    in_angles = False
    with open(data_path) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("Angles"):
                in_angles = True
                continue
            if in_angles:
                if not stripped:
                    continue
                parts = stripped.split()
                if len(parts) < 5 or not parts[0].isdigit():
                    if any(len(v) for v in triples.values()):
                        break
                    continue
                atype = int(parts[1])
                if atype in triples:
                    i, j, k = int(parts[2]), int(parts[3]), int(parts[4])
                    triples[atype].append((i, j, k))
    return triples


def load_needed_atom_positions(dump_path: str, needed_ids: set[int], min_step: int):
    """Yield (timestep, box_lengths(3,), {id: (x,y,z)}) for frames with
    timestep >= min_step, only keeping coordinates for needed_ids."""
    with open(dump_path) as fh:
        while True:
            line = fh.readline()
            if not line:
                return
            if line.strip() != "ITEM: TIMESTEP":
                continue
            timestep = int(fh.readline().strip())
            assert fh.readline().strip() == "ITEM: NUMBER OF ATOMS"
            natoms = int(fh.readline().strip())
            assert fh.readline().strip().startswith("ITEM: BOX BOUNDS")
            box = np.zeros(3)
            for d in range(3):
                lo, hi = map(float, fh.readline().split())
                box[d] = hi - lo
            header = fh.readline()
            assert header.strip().startswith("ITEM: ATOMS")
            cols = header.strip().split()[2:]
            id_idx = cols.index("id")
            x_idx, y_idx, z_idx = cols.index("x"), cols.index("y"), cols.index("z")

            if timestep < min_step:
                for _ in range(natoms):
                    fh.readline()
                continue

            coords: dict[int, tuple[float, float, float]] = {}
            for _ in range(natoms):
                parts = fh.readline().split()
                aid = int(parts[id_idx])
                if aid in needed_ids:
                    coords[aid] = (
                        float(parts[x_idx]),
                        float(parts[y_idx]),
                        float(parts[z_idx]),
                    )
            yield timestep, box, coords


def min_image(dr: np.ndarray, box: np.ndarray) -> np.ndarray:
    return dr - box * np.round(dr / box)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dump", default="dump.at_prod")
    parser.add_argument("--min-step", type=int, default=70000)
    parser.add_argument("--data", default=DATA_FILE)
    args = parser.parse_args()
    dump_file = args.dump
    min_step = args.min_step

    triples = parse_ring_angles(args.data)
    counts = {t: len(v) for t, v in triples.items()}
    print(f"Ring C-N-C angle instances parsed from {DATA_FILE}: {counts}")
    total_expected = sum(counts.values())
    if total_expected == 0:
        print("ERROR: no ring angle instances found -- check DATA_FILE/type IDs.")
        sys.exit(1)

    needed_ids: set[int] = set()
    for trip_list in triples.values():
        for i, j, k in trip_list:
            needed_ids.update((i, j, k))
    print(f"Unique atom IDs needed: {len(needed_ids)}")

    residuals_by_type: dict[int, list[float]] = {t: [] for t in RING_ANGLE_COEFFS}
    raw_angles_by_type: dict[int, list[float]] = {t: [] for t in RING_ANGLE_COEFFS}
    cc_dist_by_type: dict[int, list[float]] = {t: [] for t in RING_ANGLE_COEFFS}

    n_frames = 0
    for timestep, box, coords in load_needed_atom_positions(dump_file, needed_ids, min_step):
        n_frames += 1
        for atype, trip_list in triples.items():
            _, theta0_deg = RING_ANGLE_COEFFS[atype]
            theta0_rad = math.radians(theta0_deg)
            for i, j, k in trip_list:
                if i not in coords or j not in coords or k not in coords:
                    continue
                ri = np.array(coords[i])
                rj = np.array(coords[j])
                rk = np.array(coords[k])
                a = min_image(ri - rj, box)
                b = min_image(rk - rj, box)
                cosang = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
                cosang = max(-1.0, min(1.0, cosang))
                theta = math.acos(cosang)
                raw_angles_by_type[atype].append(theta)
                residuals_by_type[atype].append(theta - theta0_rad)

                dCC = np.linalg.norm(min_image(ri - rk, box))
                cc_dist_by_type[atype].append(dCC)
        if n_frames % 50 == 0:
            print(f"  ... processed {n_frames} frames (last timestep={timestep})", flush=True)

    print(f"\nTotal frames used (timestep >= {min_step}): {n_frames}")

    kT = KB_KCAL_PER_MOL_K * TEMP_K
    print(f"kT at {TEMP_K} K = {kT:.6f} kcal/mol\n")

    all_residuals = []
    all_cc = []
    print(
        f"{'type':>5} {'K':>10} {'theta0':>8} {'N':>8} {'mean(deg)':>10} "
        f"{'std(deg)':>9} {'Var(rad^2)':>11} {'kT/2K':>9} {'ratio':>7} | "
        f"{'C-C mean(nm)':>13} {'C-C std(nm)':>12}"
    )
    for atype in sorted(RING_ANGLE_COEFFS):
        K, theta0_deg = RING_ANGLE_COEFFS[atype]
        res = np.array(residuals_by_type[atype])
        raw = np.array(raw_angles_by_type[atype])
        cc = np.array(cc_dist_by_type[atype]) / 10.0  # Angstrom -> nm
        if len(res) == 0:
            print(f"{atype:>5}  (no samples collected)")
            continue
        var_meas = float(np.var(res))
        var_pred = kT / (2.0 * K)
        ratio = var_meas / var_pred
        print(
            f"{atype:>5} {K:>10.4f} {theta0_deg:>8.3f} {len(res):>8} "
            f"{math.degrees(np.mean(raw)):>10.3f} {math.degrees(np.std(raw)):>9.3f} "
            f"{var_meas:>11.6f} {var_pred:>9.6f} {ratio:>7.3f} | "
            f"{np.mean(cc):>13.4f} {np.std(cc):>12.4f}"
        )
        all_residuals.append(res)
        all_cc.append(cc)

    pooled_res = np.concatenate(all_residuals)
    pooled_cc = np.concatenate(all_cc)
    K = 40.000044  # identical for all three ring angle types
    var_meas = float(np.var(pooled_res))
    var_pred = kT / (2.0 * K)
    print(
        f"\nPooled across all 3 ring angle types "
        f"(N={len(pooled_res)}, {n_frames} frames x 500 molecules x 3 angles):"
    )
    print(
        f"  measured Var(theta-theta0) = {var_meas:.6f} rad^2 "
        f"(std = {math.degrees(math.sqrt(var_meas)):.3f} deg)"
    )
    print(
        f"  equipartition prediction kT/(2K) = {var_pred:.6f} rad^2 "
        f"(std = {math.degrees(math.sqrt(var_pred)):.3f} deg)"
    )
    print(f"  ratio measured/predicted = {var_meas / var_pred:.4f}")
    print(
        f"  raw C-C ring-flanking distance: mean = {np.mean(pooled_cc):.4f} nm, "
        f"std = {np.std(pooled_cc):.4f} nm  (ref peak position = 0.228 nm)"
    )
    print()
    if var_meas < 0.85 * var_pred:
        print(
            "=> Measured angle variance is SUBSTANTIALLY BELOW the equipartition "
            "prediction: genuine dynamical under-broadening (hypothesis A -- "
            "Langevin damping / fix-bm restraint / charge-asymmetry rigidity) "
            "is supported, independent of any RDF normalization issue."
        )
    elif var_meas > 1.15 * var_pred:
        print(
            "=> Measured angle variance EXCEEDS the equipartition prediction: "
            "the angle itself is not under-broadened. If the RDF peak height "
            "is still ~30% too high, that inflation must come from the RDF "
            "script's normalization/binning (hypothesis B -- box density "
            "never relaxed via NPT), not from genuine narrowing."
        )
    else:
        print(
            "=> Measured angle variance is CONSISTENT with the equipartition "
            "prediction (within ~15%): the underlying angle distribution is "
            "thermally normal. This rules out genuine dynamical narrowing "
            "(hypothesis A) as the primary driver and points at the RDF "
            "script's box-density normalization (hypothesis B) as the likely "
            "cause of the too-tall C-C(ring)/C-N peaks."
        )


if __name__ == "__main__":
    main()
