/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS
// clang-format off
BondStyle(backmap/harmonic,BondBackmapHarmonic);
// clang-format on
#else

#ifndef LMP_BOND_BACKMAP_HARMONIC_H
#define LMP_BOND_BACKMAP_HARMONIC_H

#include "bond.h"

namespace LAMMPS_NS {

class BondBackmapHarmonic : public Bond {
 public:
  BondBackmapHarmonic(class LAMMPS *);
  ~BondBackmapHarmonic() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_distance(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, double, int, int, double &) override;

 protected:
  double *k, *r0;
  int *is_cg;
  class Fix *fix_backmap;

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
