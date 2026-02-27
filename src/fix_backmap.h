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
FixStyle(backmap,FixBackmap);
// clang-format on
#else

#ifndef LMP_FIX_BACKMAP_H
#define LMP_FIX_BACKMAP_H

#include <map>
#include <vector>

#include "fix.h"

namespace LAMMPS_NS {

class FixBackmap : public Fix {
 public:
  FixBackmap(class LAMMPS *, int, char **);
  ~FixBackmap() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void post_force(int) override;
  void end_of_step() override;
  int modify_param(int, char **) override;

  // Per-atom array grow/copy/exchange for domain decomposition
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  // Ghost atom communication for lambda values
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  // Restart support
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  void *extract(const char *, int &) override;
  double memory_usage() override;

 private:
  struct MolMap {
    int cg_local;
    std::vector<int> at_local;
    double cg_mass;
  };

  int cg_type;
  double alpha;
  double lambda0;
  int nonuniform;
  int ramp_active;

  double *lambda;
  int maxatom;

  std::map<tagint, MolMap> mol_map;

  void build_molecule_map();
  void validate_masses();
};

}  // namespace LAMMPS_NS

#endif
#endif
