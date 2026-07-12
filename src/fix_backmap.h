/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under the
   GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(backmap,FixBackmap);
// clang-format on
#else

#ifndef LMP_FIX_BACKMAP_H
#define LMP_FIX_BACKMAP_H

#include <set>
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
  void pre_force(int) override;
  void post_force(int) override;
  void end_of_step() override;
  int modify_param(int, char **) override;
  double compute_scalar() override;

  // Per-atom array grow/copy/exchange for domain decomposition
  void grow_arrays(int) override;
  void copy_arrays(int, int, int) override;
  int pack_exchange(int, double *) override;
  int unpack_exchange(int, double *) override;

  // Ghost atom communication for lambda + CG force (forward) and COM (reverse)
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;

  // Restart support
  int pack_restart(int, double *) override;
  void unpack_restart(int, int) override;
  int size_restart(int) override;
  int maxsize_restart() override;

  void *extract(const char *, int &) override;
  double memory_usage() override;

 private:
  struct BeadMap {
    int cg_local;               // local index of the CG atom
    std::vector<int> at_local;  // local indices of AT atoms for this bead
    double at_mass_sum;         // sum of AT atom type masses
  };

  std::set<int> cg_types;
  double alpha;
  double lambda0;
  int nonuniform;
  int ramp_active;

  double *lambda;
  double lambda_global;  // single authoritative global lambda scalar
  int maxatom;

  // MPI-correct partner map and communication scratch
  int *atom2cg;  // per atom: local-or-ghost index of mapped CG bead, -1 if none
  double
      *com_buf;  // per CG atom (4): mass, mdx, mdy, mdz (reverse-comm scratch)
  double *cg_fwd;  // per CG atom (4): fx, fy, fz, mass (forward-comm scratch
                   // for ghosts)

  std::vector<BeadMap>
      bead_map;  // used for setup logging + validate_masses only

  bool is_cg_type(int t) const { return cg_types.count(t) > 0; }
  void build_bead_map();
  void validate_masses();
  void zero_com_buf();
};

}  // namespace LAMMPS_NS

#endif
#endif
