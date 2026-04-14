# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy"]
# ///
"""Verify backmapping simulation results for dodecane.

Checks:
  1. Thermodynamic stability (from LAMMPS log)
  2. Lambda convergence (from dump file — all atoms reach λ=1)
  3. Bond length distribution (from final data file)
  4. Angle distribution (from final data file)
  5. Density (from final data file box + masses)

Usage:
  uv run verify_backmapping.py [--log log.lammps] [--dump dump.backmap] [--data dodecane_final.data]
"""

from __future__ import annotations

import argparse
import contextlib
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np

PASS = 0
FAIL = 0
WARN = 0

# Reference values for dodecane (GROMOS united-atom)
REF_BOND_LENGTH = 1.53  # Å (C-C bond)
REF_ANGLE = 111.0  # degrees (C-C-C angle)
REF_DENSITY = 0.749  # g/cm³ at 298 K, 1 atm
BOND_TOL = 0.15  # Å tolerance on mean bond length
ANGLE_TOL = 5.0  # degrees tolerance on mean angle
DENSITY_TOL = 0.15  # g/cm³ tolerance (generous for small system + NVT)
TEMP_TOL = 50.0  # K tolerance on mean temperature


def check(condition: bool, msg: str) -> None:
    global PASS, FAIL
    if condition:
        PASS += 1
        print(f"  PASS: {msg}")
    else:
        FAIL += 1
        print(f"  FAIL: {msg}")


def warn(msg: str) -> None:
    global WARN
    WARN += 1
    print(f"  WARN: {msg}")


def parse_log(path: Path) -> dict[str, list[float]]:
    """Parse LAMMPS log file for thermo output."""
    data: dict[str, list[float]] = defaultdict(list)
    in_thermo = False
    headers: list[str] = []

    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith("Step "):
            headers = stripped.split()
            in_thermo = True
            continue
        if in_thermo:
            if stripped.startswith("Loop time") or not stripped:
                in_thermo = False
                continue
            parts = stripped.split()
            if len(parts) == len(headers):
                try:
                    for h, v in zip(headers, parts, strict=False):
                        data[h].append(float(v))
                except ValueError:
                    in_thermo = False
    return data


def verify_thermodynamics(log_path: Path, target_temp: float = 298.0) -> None:
    """Check temperature and energy stability from log file."""
    print(f"\n--- Thermodynamic stability ({log_path.name}) ---")
    data = parse_log(log_path)

    if "Temp" not in data:
        warn("No temperature data found in log file")
        return

    temps = np.array(data["Temp"])
    np.array(data["Step"])

    # Split into backmapping phase (first half) and production (second half)
    mid = len(temps) // 2
    prod_temps = temps[mid:]

    mean_t = np.mean(prod_temps)
    std_t = np.std(prod_temps)

    # If thermo_modify temp at_temp was used, reported T is the AT temperature.
    # Otherwise, CG atoms (with zero velocity) dilute the value:
    #   T_reported ~= (N_at / N_total) * T_actual
    is_correct_temp = abs(mean_t - target_temp) < TEMP_TOL
    if is_correct_temp:
        check(True, f"Production <T> = {mean_t:.1f} ± {std_t:.1f} K (target {target_temp} K)")
    else:
        # Estimate what temperature would be if CG atoms dilute it
        corrected = mean_t * 1.5  # typical N_total/N_at ratio (180/120 for dodecane)
        if abs(corrected - target_temp) < TEMP_TOL:
            check(
                True,
                f"Production <T> = {mean_t:.1f} K (AT-corrected ≈ {corrected:.0f} K; "
                f"add 'compute at_temp at_atoms temp' + 'thermo_modify temp at_temp')",
            )
        else:
            check(False, f"Production <T> = {mean_t:.1f} ± {std_t:.1f} K (target {target_temp} K)")

    # Energy should not diverge
    if "TotEng" in data:
        energies = np.array(data["TotEng"])
        prod_e = energies[mid:]
        mean_e = np.mean(prod_e)
        std_e = np.std(prod_e)
        # No energy explosion: std should be < 50% of |mean|
        ratio = std_e / max(abs(mean_e), 1.0)
        check(
            ratio < 0.5,
            f"Production <E> = {mean_e:.1f} +/- {std_e:.1f} kcal/mol (sigma/|mu|={ratio:.2f})",
        )

    if "PotEng" in data:
        pe = np.array(data["PotEng"])
        # PE should not explode (no values > 10* the production mean)
        prod_pe_mean = np.mean(pe[mid:])
        max_pe = np.max(pe)
        check(
            max_pe < 10 * abs(prod_pe_mean) + 1000,
            f"Max PE = {max_pe:.1f}, production <PE> = {prod_pe_mean:.1f} (no explosion)",
        )


def parse_dump_last_frame(path: Path) -> tuple[list[dict], int]:
    """Parse last frame of a LAMMPS custom dump for lambda values."""
    frames = []
    current: list[str] = []

    for line in path.read_text().splitlines():
        if line.startswith("ITEM: TIMESTEP"):
            if current:
                frames.append(current)
            current = []
        current.append(line)
    if current:
        frames.append(current)

    if not frames:
        return [], 0

    last = frames[-1]
    atoms = []
    in_atoms = False
    for line in last:
        if line.startswith("ITEM: NUMBER OF ATOMS"):
            continue
        if line.startswith("ITEM: ATOMS"):
            headers = line.split()[2:]  # skip "ITEM: ATOMS"
            in_atoms = True
            continue
        if in_atoms and line.strip():
            parts = line.split()
            atom = {}
            for h, v in zip(headers, parts, strict=False):
                try:
                    atom[h] = float(v)
                except ValueError:
                    atom[h] = v
            atoms.append(atom)
        if not in_atoms:
            with contextlib.suppress(ValueError):
                int(line.strip())

    timestep = 0
    if last and last[1].strip().isdigit():
        timestep = int(last[1].strip())

    return atoms, timestep


def verify_lambda(dump_path: Path) -> None:
    """Check that all atoms reached λ=1.0 in the final frame."""
    print(f"\n--- Lambda convergence ({dump_path.name}) ---")
    atoms, timestep = parse_dump_last_frame(dump_path)

    if not atoms:
        warn("No atoms found in dump file")
        return

    if "f_bm" not in atoms[0]:
        warn("No f_bm (lambda) column in dump — check dump format")
        return

    lambdas = np.array([a["f_bm"] for a in atoms])
    at_lambdas = lambdas[lambdas > 0]  # CG atoms have lambda=0

    if len(at_lambdas) == 0:
        warn("No AT atoms with lambda > 0 found")
        return

    min_l = np.min(at_lambdas)
    max_l = np.max(at_lambdas)
    mean_l = np.mean(at_lambdas)

    check(min_l >= 0.99, f"Min AT λ = {min_l:.4f} (should be ≥ 0.99)")
    check(max_l <= 1.01, f"Max AT λ = {max_l:.4f} (should be ≤ 1.01)")
    check(
        abs(mean_l - 1.0) < 0.01,
        f"Mean AT λ = {mean_l:.4f} at step {timestep} (should be ~1.0)",
    )


def parse_data_file(path: Path) -> tuple[list[dict], list[dict], list[dict], np.ndarray, dict]:
    """Parse LAMMPS data file for atoms, bonds, angles, box, masses."""
    text = path.read_text()
    lines = text.splitlines()

    atoms: list[dict] = []
    bonds: list[dict] = []
    angles: list[dict] = []
    box = np.zeros(3)
    masses: dict[int, float] = {}
    section = None

    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        if stripped.startswith("Atoms"):
            section = "atoms"
            continue
        if stripped.startswith("Bonds"):
            section = "bonds"
            continue
        if stripped.startswith("Angles"):
            section = "angles"
            continue
        if stripped.startswith("Velocities"):
            section = "vel"
            continue
        if stripped.startswith("Masses"):
            section = "masses"
            continue

        if "xlo" in stripped:
            p = stripped.split()
            box[0] = float(p[1]) - float(p[0])
        elif "ylo" in stripped:
            p = stripped.split()
            box[1] = float(p[1]) - float(p[0])
        elif "zlo" in stripped:
            p = stripped.split()
            box[2] = float(p[1]) - float(p[0])
        elif section == "masses":
            p = stripped.split()
            if len(p) >= 2:
                with contextlib.suppress(ValueError):
                    masses[int(p[0])] = float(p[1])
        elif section == "atoms":
            p = stripped.split()
            if len(p) >= 7:
                atoms.append(
                    {
                        "id": int(p[0]),
                        "mol": int(p[1]),
                        "type": int(p[2]),
                        "x": float(p[4]),
                        "y": float(p[5]),
                        "z": float(p[6]),
                    }
                )
        elif section == "bonds":
            p = stripped.split()
            if len(p) >= 4:
                bonds.append({"type": int(p[1]), "i": int(p[2]), "j": int(p[3])})
        elif section == "angles":
            p = stripped.split()
            if len(p) >= 5:
                angles.append({"type": int(p[1]), "i": int(p[2]), "j": int(p[3]), "k": int(p[4])})

    return atoms, bonds, angles, box, masses


def min_image_dist(r1: np.ndarray, r2: np.ndarray, box: np.ndarray) -> float:
    """Minimum-image distance between two positions."""
    delta = r1 - r2
    delta -= np.round(delta / box) * box
    return float(np.linalg.norm(delta))


def min_image_angle(r1: np.ndarray, r2: np.ndarray, r3: np.ndarray, box: np.ndarray) -> float:
    """Angle i-j-k using minimum image convention (in degrees)."""
    d1 = r1 - r2
    d1 -= np.round(d1 / box) * box
    d2 = r3 - r2
    d2 -= np.round(d2 / box) * box
    cos_theta = np.dot(d1, d2) / (np.linalg.norm(d1) * np.linalg.norm(d2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return float(np.degrees(np.arccos(cos_theta)))


def verify_structure(data_path: Path) -> None:
    """Check bond lengths, angles, and density from the final data file."""
    print(f"\n--- Structural properties ({data_path.name}) ---")
    atoms, bonds, angles, box, masses = parse_data_file(data_path)

    if not atoms:
        warn("No atoms found in data file")
        return

    id_to_pos = {a["id"]: np.array([a["x"], a["y"], a["z"]]) for a in atoms}

    # Identify AT atoms (types 3, 4 for dodecane: CH3, CH2)
    cg_types = {1, 2}
    at_atoms = [a for a in atoms if a["type"] not in cg_types]
    at_ids = {a["id"] for a in at_atoms}

    # Bond lengths (AT-AT bonds only)
    at_bond_lengths = []
    for b in bonds:
        if b["i"] in at_ids and b["j"] in at_ids:
            d = min_image_dist(id_to_pos[b["i"]], id_to_pos[b["j"]], box)
            at_bond_lengths.append(d)

    if at_bond_lengths:
        bl = np.array(at_bond_lengths)
        mean_bl = np.mean(bl)
        std_bl = np.std(bl)
        check(
            abs(mean_bl - REF_BOND_LENGTH) < BOND_TOL,
            f"Mean AT bond length = {mean_bl:.3f} ± {std_bl:.3f} Å (ref {REF_BOND_LENGTH} Å)",
        )
        check(
            np.max(bl) < 3.0,
            f"Max AT bond length = {np.max(bl):.3f} Å (no broken bonds)",
        )
    else:
        warn("No AT-AT bonds found for bond length analysis")

    # Angles (AT-AT-AT angles only)
    at_angles_val = []
    for ang in angles:
        if ang["i"] in at_ids and ang["j"] in at_ids and ang["k"] in at_ids:
            theta = min_image_angle(
                id_to_pos[ang["i"]],
                id_to_pos[ang["j"]],
                id_to_pos[ang["k"]],
                box,
            )
            at_angles_val.append(theta)

    if at_angles_val:
        av = np.array(at_angles_val)
        mean_a = np.mean(av)
        std_a = np.std(av)
        check(
            abs(mean_a - REF_ANGLE) < ANGLE_TOL,
            f"Mean AT angle = {mean_a:.1f} ± {std_a:.1f}° (ref {REF_ANGLE}°)",
        )
    else:
        warn("No AT-AT-AT angles found for angle analysis")

    # Density from AT atoms only
    total_mass_amu = sum(masses.get(a["type"], 0.0) for a in at_atoms)
    volume_ang3 = box[0] * box[1] * box[2]

    if volume_ang3 > 0 and total_mass_amu > 0:
        # amu / Å³ → g/cm³: multiply by 1.66054
        density = total_mass_amu / volume_ang3 * 1.66054
        check(
            abs(density - REF_DENSITY) < DENSITY_TOL,
            f"AT density = {density:.3f} g/cm³ (ref {REF_DENSITY} g/cm³, NVT — box not relaxed)",
        )
    else:
        warn("Could not compute density")

    # Chain integrity: every AT atom should be bonded
    bonded_at = set()
    for b in bonds:
        if b["i"] in at_ids:
            bonded_at.add(b["i"])
        if b["j"] in at_ids:
            bonded_at.add(b["j"])
    unbonded = at_ids - bonded_at
    check(len(unbonded) == 0, f"All {len(at_ids)} AT atoms are bonded ({len(unbonded)} unbonded)")


@dataclass
class RDFData:
    r: np.ndarray
    gr: list[np.ndarray]


RDF_PAIR_LABELS = ["CH3-CH3", "CH3-CH2", "CH2-CH2"]
RDF_PEAK_POS_TOL = 0.3  # Å
RDF_PEAK_HT_TOL = 0.35  # relative
RDF_RMSD_TOL = 0.25  # absolute


def parse_lammps_rdf(path: Path) -> RDFData:
    """Parse LAMMPS fix ave/time RDF output (mode vector), return last block."""
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


def verify_rdf(backmap_path: Path, reference_path: Path) -> None:
    """Compare backmapped RDF against reference AT RDF."""
    print(f"\n--- RDF comparison ({backmap_path.name} vs {reference_path.name}) ---")

    try:
        bm = parse_lammps_rdf(backmap_path)
        ref = parse_lammps_rdf(reference_path)
    except (ValueError, Exception) as e:
        warn(f"Could not parse RDF files: {e}")
        return

    # Interpolate onto common grid
    r_min = max(bm.r[0], ref.r[0])
    r_max = min(bm.r[-1], ref.r[-1])
    n_pts = min(len(bm.r), len(ref.r))
    r_common = np.linspace(r_min, r_max, n_pts)

    n_pairs = min(len(bm.gr), len(ref.gr))
    for i in range(n_pairs):
        label = RDF_PAIR_LABELS[i] if i < len(RDF_PAIR_LABELS) else f"pair_{i + 1}"

        g_bm = np.interp(r_common, bm.r, bm.gr[i])
        g_ref = np.interp(r_common, ref.r, ref.gr[i])

        # First peak (r > 2 Å to skip the exclusion zone)
        mask = r_common > 2.0
        r_cut, g_bm_cut, g_ref_cut = r_common[mask], g_bm[mask], g_ref[mask]

        pos_bm = r_cut[np.argmax(g_bm_cut)]
        pos_ref = r_cut[np.argmax(g_ref_cut)]
        ht_bm = float(np.max(g_bm_cut))
        ht_ref = float(np.max(g_ref_cut))

        pos_err = abs(pos_bm - pos_ref)
        check(
            pos_err < RDF_PEAK_POS_TOL,
            f"{label} 1st peak position: Δ={pos_err:.2f} Å "
            f"(backmap={pos_bm:.2f}, ref={pos_ref:.2f}, tol={RDF_PEAK_POS_TOL})",
        )

        ht_rel = abs(ht_bm - ht_ref) / max(ht_ref, 1e-6)
        check(
            ht_rel < RDF_PEAK_HT_TOL,
            f"{label} 1st peak height: rel err={ht_rel:.2f} "
            f"(backmap={ht_bm:.3f}, ref={ht_ref:.3f}, tol={RDF_PEAK_HT_TOL})",
        )

        rmsd = float(np.sqrt(np.mean((g_bm - g_ref) ** 2)))
        check(
            rmsd < RDF_RMSD_TOL,
            f"{label} RMSD(g(r)) = {rmsd:.4f} (tol={RDF_RMSD_TOL})",
        )


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify backmapping results")
    parser.add_argument("--log", type=Path, default=Path("log.lammps"))
    parser.add_argument("--dump", type=Path, default=Path("dump.backmap"))
    parser.add_argument("--data", type=Path, default=Path("dodecane_final.data"))
    parser.add_argument("--temp", type=float, default=298.0, help="Target temperature (K)")
    parser.add_argument(
        "--rdf-backmap",
        type=Path,
        default=Path("rdf_backmap.dat"),
        help="RDF from backmapping Phase 3",
    )
    parser.add_argument(
        "--rdf-reference",
        type=Path,
        default=Path("rdf_reference.dat"),
        help="RDF from standalone AT reference simulation",
    )
    args = parser.parse_args()

    if args.log.exists():
        verify_thermodynamics(args.log, args.temp)
    else:
        print(f"\nSKIP: {args.log} not found")

    if args.dump.exists():
        verify_lambda(args.dump)
    else:
        print(f"\nSKIP: {args.dump} not found")

    if args.data.exists():
        verify_structure(args.data)
    else:
        print(f"\nSKIP: {args.data} not found")

    if args.rdf_backmap.exists() and args.rdf_reference.exists():
        verify_rdf(args.rdf_backmap, args.rdf_reference)
    elif args.rdf_backmap.exists() or args.rdf_reference.exists():
        missing = args.rdf_reference if args.rdf_backmap.exists() else args.rdf_backmap
        print(f"\nSKIP: RDF comparison — {missing} not found")
        print("  To generate reference RDF:")
        print("    uv run extract_at_system.py --input dodecane_final.data")
        print("    lmp -in in.dodecane_at_ref")
    else:
        print("\nSKIP: RDF comparison — no RDF files found")

    print(f"\n{'=' * 50}")
    print(f"Results: {PASS} passed, {FAIL} failed, {WARN} warnings")
    if FAIL > 0:
        print("\nNext steps for failures:")
        print("  - Bond/angle issues → check AT placement geometry or increase equilibration")
        print("  - Lambda not converged → increase backmapping steps (reduce alpha)")
        print("  - Temperature drift → check thermostat damping or reduce timestep")
        print("  - Density off → expected for NVT (fixed box); run NPT for density match")
        print("  - RDF mismatch → increase Phase 3 run length or decrease alpha")
    if FAIL == 0 and WARN == 0:
        print("\nAll checks passed — backmapped structure matches reference AT simulation.")
    return 1 if FAIL > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
