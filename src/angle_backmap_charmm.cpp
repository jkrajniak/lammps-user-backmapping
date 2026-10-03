/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* angle_style backmap/charmm -- lambda-weighted CHARMM angle (harmonic +
   Urey-Bradley) for backmapping.

   E = w [K (theta - theta0)^2 + K_ub (r13 - r_ub)^2]

   Same form and kernel as LAMMPS angle_style charmm (MOLECULE); w is the
   compute_weight3 weight of the angle (1 intra-bead, lambda for `at` across
   beads, 1 - lambda for `cg`). GROMACS angle func 5 (theta0, k_theta, r13,
   k_UB with 1/2 k conventions) maps to K = k_theta/2, K_ub = k_UB/2.

   Syntax:
     angle_style backmap/charmm
     angle_coeff N at/cg K theta0 K_ub r_ub */

#include "angle_backmap_charmm.h"

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
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

AngleBackmapCharmm::AngleBackmapCharmm(LAMMPS *lmp)
    : Angle(lmp),
      k(nullptr),
      theta0(nullptr),
      k_ub(nullptr),
      r_ub(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

AngleBackmapCharmm::~AngleBackmapCharmm() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(theta0);
    memory->destroy(k_ub);
    memory->destroy(r_ub);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void AngleBackmapCharmm::allocate() {
  allocated = 1;
  int n = atom->nangletypes + 1;

  memory->create(setflag, n, "angle:setflag");
  memory->create(k, n, "angle:k");
  memory->create(theta0, n, "angle:theta0");
  memory->create(k_ub, n, "angle:k_ub");
  memory->create(r_ub, n, "angle:r_ub");
  memory->create(is_cg, n, "angle:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapCharmm::coeff(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR, "Incorrect args for angle_coeff backmap/charmm");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nangletypes, ilo, ihi, error);

  int cg_flag = 0;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "angle_coeff backmap/charmm: 2nd arg must be 'at' or 'cg'");

  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  double theta0_one = utils::numeric(FLERR, arg[3], false, lmp);
  double k_ub_one = utils::numeric(FLERR, arg[4], false, lmp);
  double r_ub_one = utils::numeric(FLERR, arg[5], false, lmp);

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    theta0[i] = theta0_one / 180.0 * MY_PI;
    k_ub[i] = k_ub_one;
    r_ub[i] = r_ub_one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void AngleBackmapCharmm::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "angle_style backmap/charmm");
}

/* ----------------------------------------------------------------------
   weight of one angle; 1 if no fix backmap is present
------------------------------------------------------------------------- */

double AngleBackmapCharmm::weight(int type, int i1, int i2, int i3) {
  if (!fix_backmap) return 1.0;
  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(
        FLERR,
        "angle_style backmap/charmm: cannot extract atom2cg/lambda_global");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);
  bool same_bead = BackmapLambda::same_bead(atom2cg, i1, i2, i3);
  return BackmapLambda::compute_weight3(same_bead, is_cg[type], lambda_global);
}

/* ---------------------------------------------------------------------- */

void AngleBackmapCharmm::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

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
    int type = anglelist[n][3];

    double w = weight(type, i1, i2, i3);
    if (BackmapLambda::is_almost_zero(w)) continue;

    double delx1 = x[i1][0] - x[i2][0];
    double dely1 = x[i1][1] - x[i2][1];
    double delz1 = x[i1][2] - x[i2][2];
    double rsq1 = delx1 * delx1 + dely1 * dely1 + delz1 * delz1;
    double r1 = sqrt(rsq1);

    double delx2 = x[i3][0] - x[i2][0];
    double dely2 = x[i3][1] - x[i2][1];
    double delz2 = x[i3][2] - x[i2][2];
    double rsq2 = delx2 * delx2 + dely2 * dely2 + delz2 * delz2;
    double r2 = sqrt(rsq2);

    // Urey-Bradley 1-3 distance
    double delxUB = x[i3][0] - x[i1][0];
    double delyUB = x[i3][1] - x[i1][1];
    double delzUB = x[i3][2] - x[i1][2];
    double rUB = sqrt(delxUB * delxUB + delyUB * delyUB + delzUB * delzUB);

    double dr = rUB - r_ub[type];
    double rk = w * k_ub[type] * dr;
    double forceUB = (rUB > 0.0) ? -2.0 * rk / rUB : 0.0;

    double eangle = 0.0;
    if (eflag) eangle = rk * dr;

    double c = (delx1 * delx2 + dely1 * dely2 + delz1 * delz2) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    double s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;
    s = 1.0 / s;

    double dtheta = acos(c) - theta0[type];
    double tk = w * k[type] * dtheta;
    if (eflag) eangle += tk * dtheta;

    double a = -2.0 * tk * s;
    double a11 = a * c / rsq1;
    double a12 = -a / (r1 * r2);
    double a22 = a * c / rsq2;

    double f1[3], f3[3];
    f1[0] = a11 * delx1 + a12 * delx2 - delxUB * forceUB;
    f1[1] = a11 * dely1 + a12 * dely2 - delyUB * forceUB;
    f1[2] = a11 * delz1 + a12 * delz2 - delzUB * forceUB;
    f3[0] = a22 * delx2 + a12 * delx1 + delxUB * forceUB;
    f3[1] = a22 * dely2 + a12 * dely1 + delyUB * forceUB;
    f3[2] = a22 * delz2 + a12 * delz1 + delzUB * forceUB;

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }
    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= f1[0] + f3[0];
      f[i2][1] -= f1[1] + f3[1];
      f[i2][2] -= f1[2] + f3[2];
    }
    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }

    if (evflag)
      ev_tally(i1, i2, i3, nlocal, newton_bond, eangle, f1, f3, delx1, dely1,
               delz1, delx2, dely2, delz2);
  }
}

/* ---------------------------------------------------------------------- */

double AngleBackmapCharmm::equilibrium_angle(int i) { return theta0[i]; }

/* ---------------------------------------------------------------------- */

double AngleBackmapCharmm::single(int type, int i1, int i2, int i3) {
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

  double delxUB = x[i3][0] - x[i1][0];
  double delyUB = x[i3][1] - x[i1][1];
  double delzUB = x[i3][2] - x[i1][2];
  domain->minimum_image(FLERR, delxUB, delyUB, delzUB);
  double rUB = sqrt(delxUB * delxUB + delyUB * delyUB + delzUB * delzUB);

  double c = (delx1 * delx2 + dely1 * dely2 + delz1 * delz2) / (r1 * r2);
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  double dtheta = acos(c) - theta0[type];
  double dr = rUB - r_ub[type];
  double w = weight(type, i1, i2, i3);
  return w * (k[type] * dtheta * dtheta + k_ub[type] * dr * dr);
}

/* ---------------------------------------------------------------------- */

void AngleBackmapCharmm::write_restart(FILE *fp) {
  fwrite(k + 1, sizeof(double), atom->nangletypes, fp);
  fwrite(theta0 + 1, sizeof(double), atom->nangletypes, fp);
  fwrite(k_ub + 1, sizeof(double), atom->nangletypes, fp);
  fwrite(r_ub + 1, sizeof(double), atom->nangletypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nangletypes, fp);
}

void AngleBackmapCharmm::read_restart(FILE *fp) {
  allocate();
  int n = atom->nangletypes;
  if (comm->me == 0) {
    utils::sfread(FLERR, k + 1, sizeof(double), n, fp, nullptr, error);
    utils::sfread(FLERR, theta0 + 1, sizeof(double), n, fp, nullptr, error);
    utils::sfread(FLERR, k_ub + 1, sizeof(double), n, fp, nullptr, error);
    utils::sfread(FLERR, r_ub + 1, sizeof(double), n, fp, nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), n, fp, nullptr, error);
  }
  MPI_Bcast(k + 1, n, MPI_DOUBLE, 0, world);
  MPI_Bcast(theta0 + 1, n, MPI_DOUBLE, 0, world);
  MPI_Bcast(k_ub + 1, n, MPI_DOUBLE, 0, world);
  MPI_Bcast(r_ub + 1, n, MPI_DOUBLE, 0, world);
  MPI_Bcast(is_cg + 1, n, MPI_INT, 0, world);

  for (int i = 1; i <= n; i++) setflag[i] = 1;
}

/* proc 0 writes to data file, one line per type in coeff() argument order */

void AngleBackmapCharmm::write_data(FILE *fp) {
  for (int i = 1; i <= atom->nangletypes; i++)
    fprintf(fp, "%d %s %.15g %.15g %.15g %.15g\n", i, is_cg[i] ? "cg" : "at",
            k[i], theta0[i] * 180.0 / MY_PI, k_ub[i], r_ub[i]);
}
