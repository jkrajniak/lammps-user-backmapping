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
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBackmapCapforce::FixBackmapCapforce(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg),
      fmax(0.0),
      fmaxsq(0.0),
      ramp(0.0),
      ncapped(0),
      fmax_observed(0.0) {
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix backmap/capforce", error);

  fmax = utils::numeric(FLERR, arg[3], false, lmp);
  if (fmax <= 0.0) error->all(FLERR, "fix backmap/capforce Fmax must be > 0");
  fmaxsq = fmax * fmax;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "ramp") == 0) {
      if (iarg + 1 >= narg)
        error->all(FLERR, "Illegal fix backmap/capforce command");
      ramp = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix backmap/capforce command");
    }
  }

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  extscalar = 0;
  extvector = 0;
  global_freq = 1;

  if (ramp != 0.0) {
    // Advance Fmax once per timestep after integration.
    nevery = 1;
  }
}

/* ---------------------------------------------------------------------- */

int FixBackmapCapforce::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  if (ramp != 0.0) mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::init() {}

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::min_setup(int vflag) { min_post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::post_force(int /*vflag*/) { apply_cap(); }

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::min_post_force(int /*vflag*/) { apply_cap(); }

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::apply_cap() {
  double **f = atom->f;
  int *mask = atom->mask;
  const int nlocal = atom->nlocal;
  const int groupbit = group->bitmask[igroup];

  bigint ncapped_local = 0;
  double fmax_local = 0.0;

  for (int i = 0; i < nlocal; i++) {
    if (!(mask[i] & groupbit)) continue;

    const double fx = f[i][0];
    const double fy = f[i][1];
    const double fz = f[i][2];
    const double fsq = fx * fx + fy * fy + fz * fz;
    if (fsq > fmax_local * fmax_local) fmax_local = sqrt(fsq);

    if (fsq > fmaxsq) {
      const double scaling = fmax / sqrt(fsq);
      f[i][0] *= scaling;
      f[i][1] *= scaling;
      f[i][2] *= scaling;
      ncapped_local++;
    }
  }

  bigint ncapped_all = 0;
  double fmax_all = 0.0;
  MPI_Allreduce(&ncapped_local, &ncapped_all, 1, MPI_LMP_BIGINT, MPI_SUM,
                world);
  MPI_Allreduce(&fmax_local, &fmax_all, 1, MPI_DOUBLE, MPI_MAX, world);
  ncapped = ncapped_all;
  fmax_observed = fmax_all;
}

/* ---------------------------------------------------------------------- */

void FixBackmapCapforce::end_of_step() {
  if (ramp == 0.0) return;
  fmax += ramp;
  if (fmax <= 0.0) {
    fmax = 0.0;
    fmaxsq = 0.0;
    return;
  }
  fmaxsq = fmax * fmax;
}

/* ---------------------------------------------------------------------- */

double FixBackmapCapforce::compute_scalar() {
  return static_cast<double>(ncapped);
}

/* ---------------------------------------------------------------------- */

double FixBackmapCapforce::compute_vector(int n) {
  // 0: current Fmax threshold; 1: max |F| observed before capping this step
  if (n == 0) return fmax;
  return fmax_observed;
}
