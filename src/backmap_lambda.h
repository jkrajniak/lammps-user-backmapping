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
   Provides: fix lookup, atom2cg/lambda_global access. The pure weighting
   math (compute_weight3, same_bead, clamp_lambda, is_almost_zero) lives in
   backmap_lambda_weights.h so it can be unit-tested without LAMMPS.

   Usage pattern in a style's compute():
     int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
     double lambda_global =
         BackmapLambda::clamp_lambda(*BackmapLambda::extract_lambda_global(fix_backmap));
     for each interaction (e.g. pair i,j):
       bool same = BackmapLambda::same_bead(atom2cg, i, j);
       double w  = BackmapLambda::compute_weight3(same, is_cg, lambda_global);
       if (BackmapLambda::is_almost_zero(w)) continue;
       ... compute force, scale by w ...
*/

#ifndef LMP_BACKMAP_LAMBDA_H
#define LMP_BACKMAP_LAMBDA_H

#include <cstring>

#include "backmap_lambda_weights.h"
#include "error.h"
#include "fix.h"
#include "lmptype.h"
#include "modify.h"

namespace LAMMPS_NS {
namespace BackmapLambda {

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

// Extract the per-atom CG-bead-membership map (local-or-ghost index of the
// mapped CG bead, -1 if none) from fix backmap.
inline int *extract_atom2cg(Fix *fix_backmap) {
  int dim = 0;
  auto *ptr = static_cast<int *>(fix_backmap->extract("atom2cg", dim));
  return ptr;
}

// Extract a pointer to the single authoritative global lambda scalar from
// fix backmap. Dereference to get the current value.
inline double *extract_lambda_global(Fix *fix_backmap) {
  int dim = 0;
  auto *ptr = static_cast<double *>(fix_backmap->extract("lambda_global", dim));
  return ptr;
}

}  // namespace BackmapLambda
}  // namespace LAMMPS_NS

#endif
