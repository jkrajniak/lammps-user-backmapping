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
DihedralStyle(backmap/fourier,DihedralBackmapFourier);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_BACKMAP_FOURIER_H
#define LMP_DIHEDRAL_BACKMAP_FOURIER_H

#include <vector>

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralBackmapFourier : public Dihedral {
 public:
  DihedralBackmapFourier(class LAMMPS *);
  ~DihedralBackmapFourier() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  struct Term {
    double k;
    int n;
    double shift;  // degrees
    double cos_shift, sin_shift;
  };
  std::vector<std::vector<Term> > terms;
  int *is_cg;
  class Fix *fix_backmap;

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
