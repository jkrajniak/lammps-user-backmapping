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
FixStyle(backmap/stereo,FixBackmapStereo);
// clang-format on
#else

#ifndef LMP_FIX_BACKMAP_STEREO_H
#define LMP_FIX_BACKMAP_STEREO_H

#include <vector>

#include "fix.h"

namespace LAMMPS_NS {

/* Improper-dihedral restraints that keep the handedness of stereocenters.
 *
 * Syntax: fix ID group-ID backmap/stereo K file
 *
 * Each non-comment line of file holds "atom1 atom2 atom3 atom4 phi_t" (atom
 * IDs, phi_t in degrees with the IUPAC sign of the dihedral). Every restraint
 * adds E = K [1 - cos(phi - phi_t)] which is zero at phi = phi_t. K is in
 * energy units.
 *
 * fix restrain dihedral has the same functional form, but in post_force it adds
 * forces to ghost atoms with newton_bond on, after the reverse communication,
 * so restraints that span two processors lose part of their force. Here every
 * processor that owns one of the four atoms computes the restraint, applies
 * forces to its own atoms only, and tallies its share of the energy,
 * independent of newton_bond.
 *
 * f_ID is the total restraint energy (use fix_modify energy yes to add it to
 * the potential energy, needed for minimization); f_ID[m] is the energy of
 * restraint m (1-based).
 */

class FixBackmapStereo : public Fix {
 public:
  FixBackmapStereo(class LAMMPS *, int, char **);
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void min_post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  double kspring;
  std::vector<tagint> ids;  // 4 per restraint
  std::vector<double> cos_shift, sin_shift;
  std::vector<double> energy_local, energy_all;
  bool vector_reduced;

  void read_restraints(const char *);
  double restrain(int);
};

}  // namespace LAMMPS_NS

#endif
#endif
