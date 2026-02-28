/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* angle_style backmap/harmonic — lambda-weighted harmonic angle for
   backmapping.

   E = w × 0.5 × K × (θ - θ₀)²

   Weight w uses lambda of the first and last atoms (i and k in i-j-k).

   Syntax:
     angle_style backmap/harmonic
     angle_coeff N at/cg K theta0 */

#include "angle_backmap_harmonic.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AngleBackmapHarmonic::AngleBackmapHarmonic(LAMMPS *lmp)
    : Angle(lmp),
      k(nullptr),
      theta0(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

AngleBackmapHarmonic::~AngleBackmapHarmonic() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void AngleBackmapHarmonic::allocate() {
  allocated = 1;
  int n = atom->nangletypes + 1;

  memory->create(setflag, n, "angle:setflag");
  memory->create(k, n, "angle:k");
  memory->create(theta0, n, "angle:theta0");
  memory->create(is_cg, n, "angle:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapHarmonic::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for angle_coeff backmap/harmonic");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  int cg_flag;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "angle_coeff backmap/harmonic: 2nd arg must be 'at' or 'cg'");

  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  double theta0_one = utils::numeric(FLERR, arg[3], false, lmp) * MY_PI / 180.0;

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = theta0_one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void AngleBackmapHarmonic::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "angle_style backmap/harmonic");
}

/* ---------------------------------------------------------------------- */

void AngleBackmapHarmonic::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  double *lam = BackmapLambda::extract_lambda(fix_backmap);
  if (!lam)
    error->all(FLERR, "angle_style backmap/harmonic: cannot extract lambda");

  double **x = atom->x;
  double **f = atom->f;
  int **anglelist = neighbor->anglelist;
  int nanglelist = neighbor->nanglelist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (int n = 0; n < nanglelist; n++) {
    int i1 = anglelist[n][0];
    int i2 = anglelist[n][1];
    int i3 = anglelist[n][2];
    int atype = anglelist[n][3];

    // Weight from first and last atoms of the angle triple
    double li = BackmapLambda::clamp_lambda(lam[i1]);
    double lk = BackmapLambda::clamp_lambda(lam[i3]);
    double w = BackmapLambda::compute_weight(li, lk, is_cg[atype]);

    if (BackmapLambda::is_almost_zero(w)) continue;

    // Vectors from central atom j to end atoms i and k
    double delx1 = x[i1][0] - x[i2][0];
    double dely1 = x[i1][1] - x[i2][1];
    double delz1 = x[i1][2] - x[i2][2];

    double delx2 = x[i3][0] - x[i2][0];
    double dely2 = x[i3][1] - x[i2][1];
    double delz2 = x[i3][2] - x[i2][2];

    double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    double rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    double r1 = sqrt(rsq1);
    double r2 = sqrt(rsq2);

    // cos(theta) and theta
    double c = (delx1 * delx2 + dely1 * dely2 + delz1 * delz2) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    double s = sqrt(1.0 - c * c);
    if (s < 0.001) s = 0.001;
    s = 1.0 / s;

    double dtheta = acos(c) - theta0[atype];
    double tk = w * k[atype] * dtheta;

    // Force on each atom
    double a = -tk * s;
    double a11 = a * c / rsq1;
    double a12 = -a / (r1 * r2);
    double a22 = a * c / rsq2;

    double f1[3], f3[3];
    f1[0] = a11 * delx1 + a12 * delx2;
    f1[1] = a11 * dely1 + a12 * dely2;
    f1[2] = a11 * delz1 + a12 * delz2;
    f3[0] = a22 * delx2 + a12 * delx1;
    f3[1] = a22 * dely2 + a12 * dely1;
    f3[2] = a22 * delz2 + a12 * delz1;

    f[i1][0] += f1[0];
    f[i1][1] += f1[1];
    f[i1][2] += f1[2];
    f[i2][0] -= f1[0] + f3[0];
    f[i2][1] -= f1[1] + f3[1];
    f[i2][2] -= f1[2] + f3[2];
    f[i3][0] += f3[0];
    f[i3][1] += f3[1];
    f[i3][2] += f3[2];

    double eangle = 0.0;
    if (eflag) eangle = w * k[atype] * dtheta * dtheta;

    if (evflag)
      ev_tally(i1, i2, i3, nlocal, newton_bond, eangle, f1, f3, delx1, dely1,
               delz1, delx2, dely2, delz2);
  }
}

/* ---------------------------------------------------------------------- */

double AngleBackmapHarmonic::equilibrium_angle(int i) { return theta0[i]; }

/* ---------------------------------------------------------------------- */

double AngleBackmapHarmonic::single(int atype, int i1, int i2, int i3) {
  double **x = atom->x;

  double delx1 = x[i1][0] - x[i2][0];
  double dely1 = x[i1][1] - x[i2][1];
  double delz1 = x[i1][2] - x[i2][2];
  domain->minimum_image(FLERR, delx1, dely1, delz1);
  double r1 = sqrt(delx1 * delx1 + dely1 * dely1 + delz1 * delz1);

  double delx2 = x[i3][0] - x[i2][0];
  double dely2 = x[i3][1] - x[i2][1];
  double delz2 = x[i3][2] - x[i2][2];
  domain->minimum_image(FLERR, delx2, dely2, delz2);
  double r2 = sqrt(delx2 * delx2 + dely2 * dely2 + delz2 * delz2);

  double c = (delx1 * delx2 + dely1 * dely2 + delz1 * delz2) / (r1 * r2);
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;
  double dtheta = acos(c) - theta0[atype];

  double w = 1.0;
  if (fix_backmap) {
    double *lam = BackmapLambda::extract_lambda(fix_backmap);
    if (lam) {
      double li = BackmapLambda::clamp_lambda(lam[i1]);
      double lk = BackmapLambda::clamp_lambda(lam[i3]);
      w = BackmapLambda::compute_weight(li, lk, is_cg[atype]);
    }
  }

  return w * k[atype] * dtheta * dtheta;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapHarmonic::write_restart(FILE *fp) {
  fwrite(k + 1, sizeof(double), atom->nangletypes, fp);
  fwrite(theta0 + 1, sizeof(double), atom->nangletypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nangletypes, fp);
}

void AngleBackmapHarmonic::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, k + 1, sizeof(double), atom->nangletypes, fp, nullptr,
                  error);
    utils::sfread(FLERR, theta0 + 1, sizeof(double), atom->nangletypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->nangletypes, fp, nullptr,
                  error);
  }
  MPI_Bcast(k + 1, atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(theta0 + 1, atom->nangletypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(is_cg + 1, atom->nangletypes, MPI_INT, 0, world);

  for (int i = 1; i <= atom->nangletypes; i++) setflag[i] = 1;
}
