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
DihedralStyle(ryckaert,DihedralRyckaert);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_RYCKAERT_H
#define LMP_DIHEDRAL_RYCKAERT_H

#include "dihedral.h"

namespace LAMMPS_NS {

class DihedralRyckaert : public Dihedral {
 public:
  DihedralRyckaert(class LAMMPS *);
  ~DihedralRyckaert() override;
  void compute(int, int) override;
  void coeff(int, char **) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;

 protected:
  double *c0, *c1, *c2, *c3, *c4, *c5;

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
