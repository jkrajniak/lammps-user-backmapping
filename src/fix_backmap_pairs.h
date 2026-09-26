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
FixStyle(backmap/pairs,FixBackmapPairs);
// clang-format on
#else

#ifndef LMP_FIX_BACKMAP_PAIRS_H
#define LMP_FIX_BACKMAP_PAIRS_H

#include <vector>

#include "fix.h"

namespace LAMMPS_NS {

class FixBackmapPairs : public Fix {
 public:
  FixBackmapPairs(class LAMMPS *, int, char **);
  ~FixBackmapPairs() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  double memory_usage() override;

 private:
  struct PairEntry {
    tagint id1, id2;
    double sigma, epsilon;
    double qq_scale;  // 1-4 Coulomb scale (GROMACS fudgeQQ); 0 = LJ only
    int is_cg;
  };

  std::vector<PairEntry> pairs;
  double cut;
  double cutsq;
  bool has_coulomb;
  double energy_local[2];  // LJ-14, Coulomb-14 on this rank
  class Fix *fix_backmap;

  void read_file(const char *filename);
};

}  // namespace LAMMPS_NS

#endif
#endif
