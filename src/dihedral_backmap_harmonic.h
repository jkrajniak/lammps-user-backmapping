/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef DIHEDRAL_CLASS
// clang-format off
DihedralStyle(backmap/harmonic,DihedralBackmapHarmonic);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_BACKMAP_HARMONIC_H
#define LMP_DIHEDRAL_BACKMAP_HARMONIC_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralBackmapHarmonic : public Dihedral {
 public:
  DihedralBackmapHarmonic(class LAMMPS *);
  ~DihedralBackmapHarmonic() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;

 protected:
  double *k, *cos_shift, *sin_shift;
  int *sign, *multiplicity, *is_cg;
  class Fix *fix_backmap;

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
