/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(backmap,PairBackmap);
// clang-format on
#else

#ifndef LMP_PAIR_BACKMAP_H
#define LMP_PAIR_BACKMAP_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBackmap : public Pair {
 public:
  PairBackmap(class LAMMPS *);
  ~PairBackmap() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  // Sub-styles
  class Pair *pair_at;
  class Pair *pair_cg;
  char *style_at;
  char *style_cg;

  double cut_at, cut_cg, cut_global;
  double **cut;

  // Per-type-pair classification: 0=none, 1=atomistic, 2=cg
  int **pair_kind;

  // Fix backmap for lambda access
  class Fix *fix_backmap;

  // Temporary force/energy arrays for sub-style computation
  double **f_at, **f_cg;
  double eng_vdwl_at, eng_vdwl_cg;
  double eng_coul_at, eng_coul_cg;
  double virial_at[6], virial_cg[6];
  int nmax_force;

  enum { NONE = 0, ATOMISTIC = 1, CG = 2 };

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
