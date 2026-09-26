/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* dihedral_style backmap/fourier -- lambda-weighted multi-term periodic
   dihedral for backmapping.

   E = w sum_i K_i [1 + cos(n_i phi - d_i)]

   Same form and kernel as LAMMPS dihedral_style fourier (EXTRA-MOLECULE),
   with an arbitrary phase d_i in degrees. GROMACS dihedral func 1/4/9
   (phi_s, k, n) maps term by term to d = phi_s, K = k, n = n. w is the
   compute_weight3 weight of the dihedral (1 intra-bead, lambda for `at`
   across beads, 1 - lambda for `cg`).

   Syntax:
     dihedral_style backmap/fourier
     dihedral_coeff N at/cg m K1 n1 d1 ... Km nm dm */

#include "dihedral_backmap_fourier.h"

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

/* ---------------------------------------------------------------------- */

DihedralBackmapFourier::DihedralBackmapFourier(LAMMPS *lmp)
    : Dihedral(lmp), is_cg(nullptr), fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

DihedralBackmapFourier::~DihedralBackmapFourier() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapFourier::allocate() {
  allocated = 1;
  int n = atom->ndihedraltypes + 1;

  memory->create(setflag, n, "dihedral:setflag");
  memory->create(is_cg, n, "dihedral:is_cg");
  terms.assign(n, {});
  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapFourier::coeff(int narg, char **arg) {
  if (narg < 6)
    error->all(FLERR, "Incorrect args for dihedral_coeff backmap/fourier");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  int cg_flag = 0;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "dihedral_coeff backmap/fourier: 2nd arg must be 'at' or 'cg'");

  int m = utils::inumeric(FLERR, arg[2], false, lmp);
  if (m < 1)
    error->all(FLERR,
               "dihedral_coeff backmap/fourier: number of terms must be >= 1");
  if (narg != 3 + 3 * m)
    error->all(
        FLERR,
        "dihedral_coeff backmap/fourier: expected {} args for {} terms, got {}",
        3 + 3 * m, m, narg);

  std::vector<Term> one(m);
  for (int j = 0; j < m; j++) {
    one[j].k = utils::numeric(FLERR, arg[3 + 3 * j], false, lmp);
    one[j].n = utils::inumeric(FLERR, arg[4 + 3 * j], false, lmp);
    one[j].shift = utils::numeric(FLERR, arg[5 + 3 * j], false, lmp);
    if (one[j].n < 0)
      error->all(FLERR,
                 "dihedral_coeff backmap/fourier: multiplicity must be >= 0");
    one[j].cos_shift = cos(MY_PI * one[j].shift / 180.0);
    one[j].sin_shift = sin(MY_PI * one[j].shift / 180.0);
  }

  for (int i = ilo; i <= ihi; i++) {
    terms[i] = one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapFourier::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "dihedral_style backmap/fourier");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapFourier::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(
        FLERR,
        "dihedral_style backmap/fourier: cannot extract atom2cg/lambda_global");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (int n = 0; n < ndihedrallist; n++) {
    int i1 = dihedrallist[n][0];
    int i2 = dihedrallist[n][1];
    int i3 = dihedrallist[n][2];
    int i4 = dihedrallist[n][3];
    int type = dihedrallist[n][4];

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
    double vb2xm = -vb2x;
    double vb2ym = -vb2y;
    double vb2zm = -vb2z;

    double vb3x = x[i4][0] - x[i3][0];
    double vb3y = x[i4][1] - x[i3][1];
    double vb3z = x[i4][2] - x[i3][2];

    double ax = vb1y * vb2zm - vb1z * vb2ym;
    double ay = vb1z * vb2xm - vb1x * vb2zm;
    double az = vb1x * vb2ym - vb1y * vb2xm;
    double bx = vb3y * vb2zm - vb3z * vb2ym;
    double by = vb3z * vb2xm - vb3x * vb2zm;
    double bz = vb3x * vb2ym - vb3y * vb2xm;

    double rasq = ax * ax + ay * ay + az * az;
    double rbsq = bx * bx + by * by + bz * bz;
    double rgsq = vb2xm * vb2xm + vb2ym * vb2ym + vb2zm * vb2zm;
    double rg = sqrt(rgsq);

    double rginv = 0.0, ra2inv = 0.0, rb2inv = 0.0;
    if (rg > 0) rginv = 1.0 / rg;
    if (rasq > 0) ra2inv = 1.0 / rasq;
    if (rbsq > 0) rb2inv = 1.0 / rbsq;
    double rabinv = sqrt(ra2inv * rb2inv);

    double c = (ax * bx + ay * by + az * bz) * rabinv;
    double s = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE))
      problem(FLERR, i1, i2, i3, i4);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    // p = sum_j k_j (1 + cos(n_j phi - d_j)); df = -dp/dphi (as in
    // dihedral_fourier)
    double edihedral = 0.0;
    double df = 0.0;
    for (const auto &t : terms[type]) {
      int m = t.n;
      double p_ = 1.0;
      double df1_ = 0.0, ddf1_ = 0.0;
      for (int i = 0; i < m; i++) {
        ddf1_ = p_ * c - df1_ * s;
        df1_ = p_ * s + df1_ * c;
        p_ = ddf1_;
      }
      p_ = p_ * t.cos_shift + df1_ * t.sin_shift;
      df1_ = df1_ * t.cos_shift - ddf1_ * t.sin_shift;
      df1_ *= -m;
      p_ += 1.0;
      if (m == 0) {
        p_ = 1.0 + t.cos_shift;
        df1_ = 0.0;
      }
      if (eflag) edihedral += w * t.k * p_;
      df += -w * t.k * df1_;
    }

    double fg = vb1x * vb2xm + vb1y * vb2ym + vb1z * vb2zm;
    double hg = vb3x * vb2xm + vb3y * vb2ym + vb3z * vb2zm;
    double fga = fg * ra2inv * rginv;
    double hgb = hg * rb2inv * rginv;
    double gaa = -ra2inv * rg;
    double gbb = rb2inv * rg;

    double dtfx = gaa * ax, dtfy = gaa * ay, dtfz = gaa * az;
    double dtgx = fga * ax - hgb * bx;
    double dtgy = fga * ay - hgb * by;
    double dtgz = fga * az - hgb * bz;
    double dthx = gbb * bx, dthy = gbb * by, dthz = gbb * bz;

    double sx2 = df * dtgx, sy2 = df * dtgy, sz2 = df * dtgz;

    double f1[3], f2[3], f3[3], f4[3];
    f1[0] = df * dtfx;
    f1[1] = df * dtfy;
    f1[2] = df * dtfz;
    f2[0] = sx2 - f1[0];
    f2[1] = sy2 - f1[1];
    f2[2] = sz2 - f1[2];
    f4[0] = df * dthx;
    f4[1] = df * dthy;
    f4[2] = df * dthz;
    f3[0] = -sx2 - f4[0];
    f3[1] = -sy2 - f4[1];
    f3[2] = -sz2 - f4[2];

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
      ev_tally(i1, i2, i3, i4, nlocal, newton_bond, edihedral, f1, f3, f4, vb1x,
               vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapFourier::write_restart(FILE *fp) {
  int n = atom->ndihedraltypes;
  fwrite(is_cg + 1, sizeof(int), n, fp);
  for (int i = 1; i <= n; i++) {
    int m = static_cast<int>(terms[i].size());
    fwrite(&m, sizeof(int), 1, fp);
    for (const auto &t : terms[i]) {
      fwrite(&t.k, sizeof(double), 1, fp);
      fwrite(&t.n, sizeof(int), 1, fp);
      fwrite(&t.shift, sizeof(double), 1, fp);
    }
  }
}

void DihedralBackmapFourier::read_restart(FILE *fp) {
  allocate();
  int n = atom->ndihedraltypes;
  if (comm->me == 0)
    utils::sfread(FLERR, is_cg + 1, sizeof(int), n, fp, nullptr, error);
  MPI_Bcast(is_cg + 1, n, MPI_INT, 0, world);

  for (int i = 1; i <= n; i++) {
    int m = 0;
    if (comm->me == 0)
      utils::sfread(FLERR, &m, sizeof(int), 1, fp, nullptr, error);
    MPI_Bcast(&m, 1, MPI_INT, 0, world);
    terms[i].resize(m);
    for (auto &t : terms[i]) {
      if (comm->me == 0) {
        utils::sfread(FLERR, &t.k, sizeof(double), 1, fp, nullptr, error);
        utils::sfread(FLERR, &t.n, sizeof(int), 1, fp, nullptr, error);
        utils::sfread(FLERR, &t.shift, sizeof(double), 1, fp, nullptr, error);
      }
      MPI_Bcast(&t.k, 1, MPI_DOUBLE, 0, world);
      MPI_Bcast(&t.n, 1, MPI_INT, 0, world);
      MPI_Bcast(&t.shift, 1, MPI_DOUBLE, 0, world);
      t.cos_shift = cos(MY_PI * t.shift / 180.0);
      t.sin_shift = sin(MY_PI * t.shift / 180.0);
    }
    setflag[i] = 1;
  }
}

/* proc 0 writes to data file, one line per type in coeff() argument order */

void DihedralBackmapFourier::write_data(FILE *fp) {
  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    fprintf(fp, "%d %s %d", i, is_cg[i] ? "cg" : "at",
            static_cast<int>(terms[i].size()));
    for (const auto &t : terms[i])
      fprintf(fp, " %.15g %d %.15g", t.k, t.n, t.shift);
    fprintf(fp, "\n");
  }
}
