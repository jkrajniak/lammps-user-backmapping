#!/usr/bin/env python3
"""Patch melamine_hybrid_final.data: replace the per-atom ring N (type 2)
and ring C (type 3) partial charges with their group means, testing the
polymer-chemistry panel hypothesis that topol_aa.top's single-conformer
ESP/RESP charge fit was never symmetrized across the molecule's 3-fold
arm symmetry (ring N spread 0.36e, ring C spread 0.31e -- confirmed
directly in melamine_hybrid_final.data: molecule 160's ring N charges are
-0.901187/-0.723387/-0.541197, ring C are 0.993526/0.826629/0.686578,
distinct per atom despite sharing a LAMMPS type), and that this asymmetry
is a real contributor to the ~30% too-tall C-C(ring)/C-N RDF peaks
independent of the (now ruled out) Langevin damping and fix-bm mechanisms.

Symmetrizing each group to its own mean exactly preserves each group's
charge sum (hence overall molecular neutrality) -- verified: ring N sum
unchanged at -2.165771e (3 x -0.7219236667), ring C sum unchanged at
2.506733e (3 x 0.8355776667). Only ring N/C are touched; every other atom
(exocyclic N, bridge C, O, H's -- all showing much tighter ~5-9% spread
in the chemistry panel's charge audit) is left exactly as-is, to keep this
a surgical, single-variable test.

Usage: python3 symmetrize_ring_charges.py melamine_hybrid_final.data melamine_hybrid_final_symmetrized.data
"""

import sys

RING_N_TYPE = 2
RING_C_TYPE = 3
RING_N_MEAN = -0.7219236666666666
RING_C_MEAN = 0.8355776666666667


def main():
    src, dst = sys.argv[1], sys.argv[2]
    in_atoms = False
    n_patched_n = 0
    n_patched_c = 0
    out_lines = []
    with open(src) as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("Atoms"):
                in_atoms = True
                out_lines.append(line)
                continue
            if stripped.startswith("Velocities"):
                in_atoms = False
                out_lines.append(line)
                continue
            if in_atoms and stripped:
                parts = stripped.split()
                if len(parts) >= 7 and parts[0].isdigit():
                    atype = int(parts[2])
                    if atype == RING_N_TYPE:
                        parts[3] = f"{RING_N_MEAN:.7f}"
                        n_patched_n += 1
                    elif atype == RING_C_TYPE:
                        parts[3] = f"{RING_C_MEAN:.7f}"
                        n_patched_c += 1
                    out_lines.append(" ".join(parts) + "\n")
                    continue
            out_lines.append(line)

    with open(dst, "w") as fh:
        fh.writelines(out_lines)

    print(f"Patched {n_patched_n} ring-N (type {RING_N_TYPE}) atoms -> charge {RING_N_MEAN:.7f}")
    print(f"Patched {n_patched_c} ring-C (type {RING_C_TYPE}) atoms -> charge {RING_C_MEAN:.7f}")
    print("Expected: 1500 each (500 molecules x 3 ring N/C atoms)")
    print(f"Wrote {dst}")


if __name__ == "__main__":
    main()
