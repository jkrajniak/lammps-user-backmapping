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

/* Repulsive sphere around ring centroids that keeps atoms from passing through
 * rings while the atomistic nonbonded repulsion is still scaled down.
 *
 * Syntax: fix ID group-ID backmap/ringguard K rc rcb file
 *
 * Each non-comment line of file lists the atom IDs of one ring, optionally
 * followed by | and the IDs of the atoms bonded to the ring (its hydrogens and
 * the first atom of each substituent). An atom of the group that is not in
 * the ring and lies at distance d < rc from the ring centroid (unweighted mean
 * of the ring atoms) adds E = K (rc - d)^2, an atom bonded to the ring
 * E = K (rcb - d)^2 for d < rcb; the reaction force is shared equally by the
 * ring atoms.
 *
 * A bond can sweep through a ring while both its atoms stay outside the
 * sphere, so each bond of the group (neither atom in the ring) also adds
 * w(s) - w(d_near), w(x) = K (rw - x)^2 for x < rw, where s is the distance of
 * the bond segment from the centroid and d_near that of the nearer bond atom.
 * The term is zero when the nearest point of the segment is an end, so the
 * energy is continuous; rw = rcb when a bond atom is bonded to the ring.
 *
 * Each atom term is evaluated by the owner of the atom, each bond term by the
 * owner of the atom that stores the bond, against the centroids of the rings
 * whose atoms are all present as owned or ghost atoms (a periodic cell grid).
 * A ring is built around the owned image of its first atom, or the ghost image
 * nearest the subdomain, and skipped if its atoms are not all minimum images
 * of that anchor (a ghost at the edge of the shell can see the other ring
 * atoms only across the box).
 * Neighbor lists are not used: they leave out special pairs with weight 0. The
 * forces go into a per-atom buffer on owned and ghost atoms, which a reverse
 * communication of the fix sums onto the owners, so the result does not
 * depend on the decomposition.
 *
 * f_ID is the energy, f_ID[1] the energy, f_ID[2] the number of atom-ring and
 * bond-ring contacts inside the sphere. The virial is not tallied.
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
  double kwall, rc, rcb;
  std::vector<tagint> ring_atoms;  // all rings, concatenated
  std::vector<int> ring_start;     // ring m is [ring_start[m], ring_start[m+1])
  std::vector<tagint> excl_atoms;  // atoms bonded to each ring, concatenated
  std::vector<int> excl_start;     // same layout as ring_start
  double **fbuf;
  int nmax;
  double energy_local, contacts_local;
  double sums[2];
  bool reduced;
  std::vector<int> ring_idx;     // local index of each ring atom, as ring_atoms
  std::vector<double> ring_cen;  // centroid of ring m at [3m, 3m+3)
  std::vector<int> cell_head;    // first ring in each cell of the centroid grid
  std::vector<int> cell_next;    // next ring in the same cell
  double cell_lo[3], cell_w[3];
  int cell_n[3];

  void read_rings(const char *);
  void build_rings();
  int ring_anchor(tagint) const;
  int cell_of(double, int) const;
  void rings_near(const double *, const double *, std::vector<int> &) const;
  int relation(int, tagint) const;
  void push_ring(int, const double *);
};

}  // namespace LAMMPS_NS

#endif
#endif
