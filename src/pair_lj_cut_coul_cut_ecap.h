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
PairStyle(lj/cut/coul/cut/ecap,PairLJCutCoulCutEcap);
// clang-format on
#else

#ifndef LMP_PAIR_LJ_CUT_COUL_CUT_ECAP_H
#define LMP_PAIR_LJ_CUT_COUL_CUT_ECAP_H

#include "pair_lj_cut_coul_cut.h"

namespace LAMMPS_NS {

/* ESPResSo++ LennardJonesEnergyCapped port for AT LJ inside backmap.

   For r <= caprad (= cap_factor * sigma): LJ force = 0, LJ energy = V(caprad).
   Coulomb is unchanged by default (matches bakery's tools_sim.py port).
   Default cap_factor = 0.5 (bakery tools_sim.py).

   Optional 4th arg coul_cap_radius (default 0.0 = disabled) is NOT part of
   the bakery/E++ reference: it caps same-pair Coulomb the same way (force=0,
   energy held flat at V(coul_cap_radius)) for r <= coul_cap_radius. Added to
   absorb unrelaxed same-fragment 1-5+ electrostatic strain that plain LJ
   capping cannot reach (see
   research/experiments/20260802_melamine-bakery-rerun.md). Leave at 0.0 to
   reproduce the original bakery-faithful behavior exactly.

   Syntax: pair_style lj/cut/coul/cut/ecap cut_lj [cut_coul [cap_factor
   [coul_cap_radius]]]
*/

class PairLJCutCoulCutEcap : public PairLJCutCoulCut {
 public:
  PairLJCutCoulCutEcap(class LAMMPS *);
  ~PairLJCutCoulCutEcap() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  double init_one(int, int) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  double cap_factor;
  double **capradsq;
  double **lj_cap_eng;
  double coul_cap_radius;
  double coul_capradsq;

  void allocate() override;
};

}  // namespace LAMMPS_NS

#endif
#endif
