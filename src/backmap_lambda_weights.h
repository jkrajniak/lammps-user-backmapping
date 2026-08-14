/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* Pure lambda-weighting math shared by all backmap styles.

   Deliberately free of any LAMMPS headers/types so it can be compiled and
   unit-tested standalone (see tests/unit/) without building LAMMPS.

   Weighting uses a single global lambda scalar (not a per-particle
   product) and per-atom CG-bead membership, giving three cases:
     - CG-CG (is_cg=true):                       w = 1 - lambda_global
     - AT-AT intra-bead (same CG bead):          w = 1 always
     - AT-AT inter-bead (different CG beads):    w = lambda_global

   The intra-bead case is deliberately NOT gated by lambda: it represents
   real intra-molecular chemistry (bonds/angles/dihedrals/LJ within the
   same fragment) that exists independent of the CG/AT resolution ramp,
   matching the original ESPResSo++ production driver
   (bakery/src/start_backmapping.py), which builds the full AT interaction
   list once, unconditionally, before the resolution ramp is ever
   activated. Only the *inter-bead* (intermolecular) AT-AT term and the
   CG-CG term are actually resolution-dependent. See notebook
   2026-07-12_at-intrabead-always-on-fix.md. */

#ifndef LMP_BACKMAP_LAMBDA_WEIGHTS_H
#define LMP_BACKMAP_LAMBDA_WEIGHTS_H

#include <cmath>

namespace LAMMPS_NS {
namespace BackmapLambda {

static constexpr double WEIGHT_ZERO_THRESHOLD = 1.0e-10;

// Clamp lambda: negative values (e.g. from a negative lambda0) → 0, cap at 1.
inline double clamp_lambda(double val) {
  if (val < 0.0) return 0.0;
  if (val > 1.0) return 1.0;
  return val;
}

// True iff all given local-or-ghost atom indices map to the same CG bead
// in atom2cg. Used to distinguish AT-AT intra-bead from inter-bead
// interactions. A negative bead index (unmapped, e.g. a CG atom itself)
// never counts as "same" via this check.
inline bool same_bead(const int *atom2cg, int i1, int i2) {
  int b = atom2cg[i1];
  return b >= 0 && b == atom2cg[i2];
}

inline bool same_bead(const int *atom2cg, int i1, int i2, int i3) {
  int b = atom2cg[i1];
  return b >= 0 && b == atom2cg[i2] && b == atom2cg[i3];
}

inline bool same_bead(const int *atom2cg, int i1, int i2, int i3, int i4) {
  int b = atom2cg[i1];
  return b >= 0 && b == atom2cg[i2] && b == atom2cg[i3] && b == atom2cg[i4];
}

// Compute the weight for an interaction using a single global lambda
// scalar (not a per-particle product) and CG-bead co-membership:
//   CG term (is_cg=true):                       w = 1 - lambda_global
//   AT intra-bead term (same_bead=true):         w = 1 always
//   AT inter-bead term (same_bead=false):        w = lambda_global
inline double compute_weight3(bool is_same_bead, bool is_cg,
                              double lambda_global) {
  if (is_cg) return 1.0 - lambda_global;
  if (is_same_bead) return 1.0;
  return lambda_global;
}

// Returns true when the weight is negligible and the interaction
// can be skipped for efficiency.
inline bool is_almost_zero(double w) {
  return std::fabs(w) < WEIGHT_ZERO_THRESHOLD;
}

}  // namespace BackmapLambda
}  // namespace LAMMPS_NS

#endif
