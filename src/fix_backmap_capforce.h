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

class FixBackmapCapforce : public Fix {
 public:
  FixBackmapCapforce(class LAMMPS *, int, char **);
  int setmask() override;
  void post_force(int) override;

 private:
  double fmax;
  double fmaxsq;
};

}  // namespace LAMMPS_NS

#endif
#endif
