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
FixStyle(backmap/ringguard,FixBackmapRingGuard);
// clang-format on
#else

#ifndef LMP_FIX_BACKMAP_RINGGUARD_H
#define LMP_FIX_BACKMAP_RINGGUARD_H

#include <vector>

#include "fix.h"

namespace LAMMPS_NS {

/* Repulsive wall around ring centroids that keeps atoms from passing through
 * rings while the atomistic nonbonded repulsion is still scaled down.
 *
 * Syntax: fix ID group-ID backmap/ringguard K rc file
 *
 * Each non-comment line of file lists the atom IDs of one ring, optionally
 * followed by | and the IDs of atoms the ring does not repel (its own
 * hydrogens and substituents, which sit close to the centroid). Any other atom
 * of the group at distance d < rc from the ring centroid (unweighted mean of
 * the ring atoms) adds E = K (rc - d)^2; the reaction force is shared equally
 * by the ring atoms.
 *
 * Each ring is evaluated once, by the processor that owns its first atom. The
 * atoms near the centroid are found with a cell grid over owned and ghost
 * atoms, not a neighbor list, because neighbor lists leave out special pairs
 * with weight 0 (an atom two bonds from a ring atom is one). Forces go into a
 * per-atom buffer on owned and ghost atoms, which a reverse communication of
 * the fix sums onto the owners, so the result does not depend on the
 * decomposition.
 *
 * f_ID is the energy, f_ID[1] the energy, f_ID[2] the number of atom-ring
 * contacts inside rc. The virial is not tallied.
 */

class FixBackmapRingGuard : public Fix {
 public:
  FixBackmapRingGuard(class LAMMPS *, int, char **);
  ~FixBackmapRingGuard() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;
  int pack_reverse_comm(int, int, double *) override;
  void unpack_reverse_comm(int, int *, double *) override;
  double memory_usage() override;

 private:
  double kwall, rc;
  std::vector<tagint> ring_atoms;  // all rings, concatenated
  std::vector<int> ring_start;     // ring m is [ring_start[m], ring_start[m+1])
  std::vector<tagint> excl_atoms;  // atoms not repelled, per ring, concatenated
  std::vector<int> excl_start;     // same layout as ring_start
  double **fbuf;
  int nmax;
  double energy_local, contacts_local;
  double sums[2];
  bool reduced;
  std::vector<int> cell_head, cell_next;
  double cell_lo[3];
  int cell_n[3];

  void read_rings(const char *);
  void build_cells();
};

}  // namespace LAMMPS_NS

#endif
#endif
