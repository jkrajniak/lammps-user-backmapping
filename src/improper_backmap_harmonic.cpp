/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* improper_style backmap/harmonic -- lambda-weighted harmonic improper for
   backmapping.

   E = w K (chi - chi0)^2

   Same form and kernel as LAMMPS improper_style harmonic (MOLECULE); chi is
   the angle between the planes (I,J,K) and (J,K,L). GROMACS improper func 2
   (xi0, k_xi with E = 1/2 k (xi - xi0)^2) maps to K = k_xi/2, chi0 = xi0.
   w is the compute_weight3 weight (1 intra-bead, lambda for `at` across
   beads, 1 - lambda for `cg`).

   Syntax:
     improper_style backmap/harmonic
     improper_coeff N at/cg K chi0 */

#include "improper_backmap_harmonic.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double TOLERANCE = 0.05;
static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

ImproperBackmapHarmonic::ImproperBackmapHarmonic(LAMMPS *lmp)
    : Improper(lmp),
      k(nullptr),
      chi(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

ImproperBackmapHarmonic::~ImproperBackmapHarmonic() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(chi);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperBackmapHarmonic::allocate() {
  allocated = 1;
  int n = atom->nimpropertypes + 1;

  memory->create(setflag, n, "improper:setflag");
  memory->create(k, n, "improper:k");
  memory->create(chi, n, "improper:chi");
  memory->create(is_cg, n, "improper:is_cg");
  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void ImproperBackmapHarmonic::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for improper_coeff backmap/harmonic");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nimpropertypes, ilo, ihi, error);

  int cg_flag = 0;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "improper_coeff backmap/harmonic: 2nd arg must be 'at' or 'cg'");

  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  double chi_one = utils::numeric(FLERR, arg[3], false, lmp);

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    chi[i] = DEG2RAD * chi_one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void ImproperBackmapHarmonic::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "improper_style backmap/harmonic");
}

/* ---------------------------------------------------------------------- */

void ImproperBackmapHarmonic::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "improper_style backmap/harmonic: cannot extract "
               "atom2cg/lambda_global");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

  double **x = atom->x;
  double **f = atom->f;
  int **improperlist = neighbor->improperlist;
  int nimproperlist = neighbor->nimproperlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (int n = 0; n < nimproperlist; n++) {
    int i1 = improperlist[n][0];
    int i2 = improperlist[n][1];
    int i3 = improperlist[n][2];
    int i4 = improperlist[n][3];
    int type = improperlist[n][4];

    bool same_bead = BackmapLambda::same_bead(atom2cg, i1, i2, i3, i4);
    double w =
        BackmapLambda::compute_weight3(same_bead, is_cg[type], lambda_global);
    if (BackmapLambda::is_almost_zero(w)) continue;

    double vb1x = x[i1][0] - x[i2][0];
    double vb1y = x[i1][1] - x[i2][1];
    double vb1z = x[i1][2] - x[i2][2];

    double vb2x = x[i3][0] - x[i2][0];
    double vb2y = x[i3][1] - x[i2][1];
    double vb2z = x[i3][2] - x[i2][2];

    double vb3x = x[i4][0] - x[i3][0];
    double vb3y = x[i4][1] - x[i3][1];
    double vb3z = x[i4][2] - x[i3][2];

    double ss1 = 1.0 / (vb1x * vb1x + vb1y * vb1y + vb1z * vb1z);
    double ss2 = 1.0 / (vb2x * vb2x + vb2y * vb2y + vb2z * vb2z);
    double ss3 = 1.0 / (vb3x * vb3x + vb3y * vb3y + vb3z * vb3z);

    double r1 = sqrt(ss1);
    double r2 = sqrt(ss2);
    double r3 = sqrt(ss3);

    double c0 = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * r1 * r3;
    double c1 = (vb1x * vb2x + vb1y * vb2y + vb1z * vb2z) * r1 * r2;
    double c2 = -(vb3x * vb2x + vb3y * vb2y + vb3z * vb2z) * r3 * r2;

    double s1 = 1.0 - c1 * c1;
    if (s1 < SMALL) s1 = SMALL;
    s1 = 1.0 / s1;

    double s2 = 1.0 - c2 * c2;
    if (s2 < SMALL) s2 = SMALL;
    s2 = 1.0 / s2;

    double s12 = sqrt(s1 * s2);
    double c = (c1 * c2 + c0) * s12;

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE))
      problem(FLERR, i1, i2, i3, i4);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    double s = sqrt(1.0 - c * c);
    if (s < SMALL) s = SMALL;

    double domega = acos(c) - chi[type];
    double a = w * k[type] * domega;

    double eimproper = 0.0;
    if (eflag) eimproper = a * domega;

    a = -a * 2.0 / s;
    c = c * a;
    s12 = s12 * a;
    double a11 = c * ss1 * s1;
    double a22 = -ss2 * (2.0 * c0 * s12 - c * (s1 + s2));
    double a33 = c * ss3 * s2;
    double a12 = -r1 * r2 * (c1 * c * s1 + c2 * s12);
    double a13 = -r1 * r3 * s12;
    double a23 = r2 * r3 * (c2 * c * s2 + c1 * s12);

    double sx2 = a22 * vb2x + a23 * vb3x + a12 * vb1x;
    double sy2 = a22 * vb2y + a23 * vb3y + a12 * vb1y;
    double sz2 = a22 * vb2z + a23 * vb3z + a12 * vb1z;

    double f1[3], f2[3], f3[3], f4[3];
    f1[0] = a12 * vb2x + a13 * vb3x + a11 * vb1x;
    f1[1] = a12 * vb2y + a13 * vb3y + a11 * vb1y;
    f1[2] = a12 * vb2z + a13 * vb3z + a11 * vb1z;

    f2[0] = -sx2 - f1[0];
    f2[1] = -sy2 - f1[1];
    f2[2] = -sz2 - f1[2];

    f4[0] = a23 * vb2x + a33 * vb3x + a13 * vb1x;
    f4[1] = a23 * vb2y + a33 * vb3y + a13 * vb1y;
    f4[2] = a23 * vb2z + a33 * vb3z + a13 * vb1z;

    f3[0] = sx2 - f4[0];
    f3[1] = sy2 - f4[1];
    f3[2] = sz2 - f4[2];

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += f1[0];
      f[i1][1] += f1[1];
      f[i1][2] += f1[2];
    }
    if (newton_bond || i2 < nlocal) {
      f[i2][0] += f2[0];
      f[i2][1] += f2[1];
      f[i2][2] += f2[2];
    }
    if (newton_bond || i3 < nlocal) {
      f[i3][0] += f3[0];
      f[i3][1] += f3[1];
      f[i3][2] += f3[2];
    }
    if (newton_bond || i4 < nlocal) {
      f[i4][0] += f4[0];
      f[i4][1] += f4[1];
      f[i4][2] += f4[2];
    }

    if (evflag)
      ev_tally(i1, i2, i3, i4, nlocal, newton_bond, eimproper, f1, f3, f4, vb1x,
               vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void ImproperBackmapHarmonic::write_restart(FILE *fp) {
  fwrite(k + 1, sizeof(double), atom->nimpropertypes, fp);
  fwrite(chi + 1, sizeof(double), atom->nimpropertypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nimpropertypes, fp);
}

void ImproperBackmapHarmonic::read_restart(FILE *fp) {
  allocate();
  int n = atom->nimpropertypes;
  if (comm->me == 0) {
    utils::sfread(FLERR, k + 1, sizeof(double), n, fp, nullptr, error);
    utils::sfread(FLERR, chi + 1, sizeof(double), n, fp, nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), n, fp, nullptr, error);
  }
  MPI_Bcast(k + 1, n, MPI_DOUBLE, 0, world);
  MPI_Bcast(chi + 1, n, MPI_DOUBLE, 0, world);
  MPI_Bcast(is_cg + 1, n, MPI_INT, 0, world);
  for (int i = 1; i <= n; i++) setflag[i] = 1;
}

/* proc 0 writes to data file, one line per type in coeff() argument order */

void ImproperBackmapHarmonic::write_data(FILE *fp) {
  for (int i = 1; i <= atom->nimpropertypes; i++)
    fprintf(fp, "%d %s %.15g %.15g\n", i, is_cg[i] ? "cg" : "at", k[i],
            RAD2DEG * chi[i]);
}
