/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* Shared lambda-access utilities for all backmap styles.
   Provides: fix lookup, per-atom lambda access, weight computation,
   and negligible-weight check.

   Usage pattern in a style's compute():
     double *lambda = BackmapLambda::extract_lambda(fix_backmap);
     for each interaction:
       double li = BackmapLambda::clamp_lambda(lambda[i]);
       double lj = BackmapLambda::clamp_lambda(lambda[j]);
       double w  = BackmapLambda::compute_weight(li, lj, is_cg);
       if (BackmapLambda::is_almost_zero(w)) continue;
       ... compute force, scale by w ...
*/

#ifndef LMP_BACKMAP_LAMBDA_H
#define LMP_BACKMAP_LAMBDA_H

#include <cmath>
#include <cstring>

#include "error.h"
#include "fix.h"
#include "lmptype.h"
#include "modify.h"

namespace LAMMPS_NS {
namespace BackmapLambda {

static constexpr double WEIGHT_ZERO_THRESHOLD = 1.0e-10;

// Scan the fix list for a fix with style "backmap". Aborts if not found.
inline Fix *find_fix_backmap(LAMMPS *lmp, const char *caller) {
  for (int i = 0; i < lmp->modify->nfix; i++) {
    if (strcmp(lmp->modify->fix[i]->style, "backmap") == 0)
      return lmp->modify->fix[i];
  }
  lmp->error->all(FLERR, "{} requires fix backmap", caller);
  return nullptr;
}

// Extract the per-atom lambda array from fix backmap.
// Call once per compute() invocation; the pointer is valid for
// the current timestep (may move on grow_arrays).
inline double *extract_lambda(Fix *fix_backmap) {
  int dim = 0;
  auto *ptr = static_cast<double *>(fix_backmap->extract("lambda", dim));
  return ptr;
}

// Clamp lambda: negative values (from nonuniform init) → 0, cap at 1.
inline double clamp_lambda(double val) {
  if (val < 0.0) return 0.0;
  if (val > 1.0) return 1.0;
  return val;
}

// Compute the weight for a pair interaction.
//   AT style (is_cg=false):  w = λ_i × λ_j
//   CG style (is_cg=true):   w = 1 − λ_i × λ_j
inline double compute_weight(double lambda_i, double lambda_j, bool is_cg) {
  double w = lambda_i * lambda_j;
  return is_cg ? (1.0 - w) : w;
}

// Returns true when the weight is negligible and the interaction
// can be skipped for efficiency.
inline bool is_almost_zero(double w) {
  return std::fabs(w) < WEIGHT_ZERO_THRESHOLD;
}

}  // namespace BackmapLambda
}  // namespace LAMMPS_NS

#endif
