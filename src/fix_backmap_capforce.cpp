/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_backmap_capforce.h"

#include <cmath>

#include "atom.h"
#include "error.h"
#include "group.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBackmapCapforce::FixBackmapCapforce(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), fmax(0.0), fmaxsq(0.0) {
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix backmap/capforce", error);

  fmax = utils::numeric(FLERR, arg[3], false, lmp);
  if (fmax <= 0.0) error->all(FLERR, "fix backmap/capforce Fmax must be > 0");
  fmaxsq = fmax * fmax;
}

/* ---------------------------------------------------------------------- */

int FixBackmapCapforce::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::post_force(int /*vflag*/) {
  double **f = atom->f;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int groupbit = group->bitmask[igroup];

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    const double fx = f[i][0];
    const double fy = f[i][1];
    const double fz = f[i][2];
    const double fsq = fx * fx + fy * fy + fz * fz;
    if (fsq > fmaxsq) {
      const double scaling = fmax / sqrt(fsq);
      f[i][0] *= scaling;
      f[i][1] *= scaling;
      f[i][2] *= scaling;
    }
  }
}
