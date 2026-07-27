/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(backmap/capforce,FixBackmapCapforce);
// clang-format on
#else

#ifndef LMP_FIX_BACKMAP_CAPFORCE_H
#define LMP_FIX_BACKMAP_CAPFORCE_H

#include "fix.h"

namespace LAMMPS_NS {

/* Absolute per-atom force clamp (ESPResSo++ CapForce port).
 *
 * Syntax: fix ID group-ID backmap/capforce Fmax [ramp delta]
 *
 * Fmax and delta are in LAMMPS force units (kcal/(mol·Å) for units real).
 * If |F| > Fmax, F is rescaled to Fmax. Optional ramp adds delta to Fmax
 * each timestep (E++ dacron uses +10 kJ/(mol·nm)/step ≈ +0.239 in real units).
 *
 * Define after fix backmap so CG→AT force redistribution is included.
 * f_ID returns the number of atoms capped on the last post_force call.
 */

class FixBackmapCapforce : public Fix {
 public:
  FixBackmapCapforce(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  void end_of_step() override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  double fmax;
  double fmaxsq;
  double ramp;
  bigint ncapped;
  double fmax_observed;

  void apply_cap();
};

}  // namespace LAMMPS_NS

#endif
#endif
