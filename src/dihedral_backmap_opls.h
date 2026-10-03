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
DihedralStyle(backmap/opls,DihedralBackmapOpls);
// clang-format on
#else

#ifndef LMP_DIHEDRAL_BACKMAP_OPLS_H
#define LMP_DIHEDRAL_BACKMAP_OPLS_H

#include "dihedral_backmap_ryckaert.h"

namespace LAMMPS_NS {

class DihedralBackmapOpls : public DihedralBackmapRyckaert {
 public:
  DihedralBackmapOpls(class LAMMPS *);
  ~DihedralBackmapOpls() override;
  void coeff(int, char **) override;
  void init_style() override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_data(FILE *) override;

 protected:
  double *k1, *k2, *k3, *k4;

  void allocate() override;
  void set_rb(int);
};

}  // namespace LAMMPS_NS

#endif
#endif
