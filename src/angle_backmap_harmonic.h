/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(backmap/harmonic,AngleBackmapHarmonic);
// clang-format on
#else

#ifndef LMP_ANGLE_BACKMAP_HARMONIC_H
#define LMP_ANGLE_BACKMAP_HARMONIC_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleBackmapHarmonic : public Angle {
 public:
  AngleBackmapHarmonic(class LAMMPS *);
  ~AngleBackmapHarmonic() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int) override;

 protected:
  double *k, *theta0;
  int *is_cg;
  class Fix *fix_backmap;

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
