/* ----------------------------------------------------------------------
   dihedral_style backmap/harmonic — lambda-weighted OPLS improper dihedral.

   E = w × K × (1 + cos(n φ - shift))
   w uses global lambda and CG-bead co-membership across all four atoms.

   Syntax:
     dihedral_style backmap/harmonic
     dihedral_coeff N at/cg K sign multiplicity
------------------------------------------------------------------------- */

#include "dihedral_backmap_harmonic.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "utils.h"

using namespace LAMMPS_NS;

static constexpr double TOLERANCE = 0.05;

/* ---------------------------------------------------------------------- */

DihedralBackmapHarmonic::DihedralBackmapHarmonic(LAMMPS *lmp)
    : Dihedral(lmp),
      k(nullptr),
      cos_shift(nullptr),
      sin_shift(nullptr),
      sign(nullptr),
      multiplicity(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralBackmapHarmonic::~DihedralBackmapHarmonic() {
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(sign);
    memory->destroy(multiplicity);
    memory->destroy(cos_shift);
    memory->destroy(sin_shift);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapHarmonic::allocate() {
  allocated = 1;
  int n = atom->ndihedraltypes + 1;

  memory->create(setflag, n, "dihedral:setflag");
  memory->create(k, n, "dihedral:k");
  memory->create(sign, n, "dihedral:sign");
  memory->create(multiplicity, n, "dihedral:multiplicity");
  memory->create(cos_shift, n, "dihedral:cos_shift");
  memory->create(sin_shift, n, "dihedral:sin_shift");
  memory->create(is_cg, n, "dihedral:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapHarmonic::coeff(int narg, char **arg) {
  if (narg != 5)
    error->all(FLERR, "Incorrect args for dihedral_coeff backmap/harmonic");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  int cg_flag;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "dihedral_coeff backmap/harmonic: 2nd arg must be 'at' or 'cg'");

  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  int sign_one = utils::inumeric(FLERR, arg[3], false, lmp);
  int multiplicity_one = utils::inumeric(FLERR, arg[4], false, lmp);

  if (sign_one != -1 && sign_one != 1)
    error->all(FLERR, "Incorrect sign arg for dihedral_coeff backmap/harmonic");
  if (multiplicity_one < 0)
    error->all(
        FLERR,
        "Incorrect multiplicity arg for dihedral_coeff backmap/harmonic");

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    sign[i] = sign_one;
    if (sign[i] == 1) {
      cos_shift[i] = 1;
      sin_shift[i] = 0;
    } else {
      cos_shift[i] = -1;
      sin_shift[i] = 0;
    }
    multiplicity[i] = multiplicity_one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapHarmonic::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "dihedral_style backmap/harmonic");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapHarmonic::compute(int eflag, int vflag) {
  int i1, i2, i3, i4, i, m, n, type;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z, vb2xm, vb2ym,
      vb2zm;
  double edihedral, f1[3], f2[3], f3[3], f4[3];
  double ax, ay, az, bx, by, bz, rasq, rbsq, rgsq, rg, rginv, ra2inv, rb2inv,
      rabinv;
  double df, df1, ddf1, fg, hg, fga, hgb, gaa, gbb;
  double dtfx, dtfy, dtfz, dtgx, dtgy, dtgz, dthx, dthy, dthz;
  double c, s, p, sx2, sy2, sz2;

  edihedral = 0.0;
  ev_init(eflag, vflag);

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "dihedral_style backmap/harmonic: cannot extract atom2cg/"
               "lambda_global");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

  double **x = atom->x;
  double **f = atom->f;
  int **dihedrallist = neighbor->dihedrallist;
  int ndihedrallist = neighbor->ndihedrallist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (n = 0; n < ndihedrallist; n++) {
    i1 = dihedrallist[n][0];
    i2 = dihedrallist[n][1];
    i3 = dihedrallist[n][2];
    i4 = dihedrallist[n][3];
    type = dihedrallist[n][4];

    bool same_bead = BackmapLambda::same_bead(atom2cg, i1, i2, i3, i4);
    double w =
        BackmapLambda::compute_weight3(same_bead, is_cg[type], lambda_global);

    if (BackmapLambda::is_almost_zero(w)) continue;

    vb1x = x[i1][0] - x[i2][0];
    vb1y = x[i1][1] - x[i2][1];
    vb1z = x[i1][2] - x[i2][2];

    vb2x = x[i3][0] - x[i2][0];
    vb2y = x[i3][1] - x[i2][1];
    vb2z = x[i3][2] - x[i2][2];

    vb2xm = -vb2x;
    vb2ym = -vb2y;
    vb2zm = -vb2z;

    vb3x = x[i4][0] - x[i3][0];
    vb3y = x[i4][1] - x[i3][1];
    vb3z = x[i4][2] - x[i3][2];

    ax = vb1y * vb2zm - vb1z * vb2ym;
    ay = vb1z * vb2xm - vb1x * vb2zm;
    az = vb1x * vb2ym - vb1y * vb2xm;
    bx = vb3y * vb2zm - vb3z * vb2ym;
    by = vb3z * vb2xm - vb3x * vb2zm;
    bz = vb3x * vb2ym - vb3y * vb2xm;

    rasq = ax * ax + ay * ay + az * az;
    rbsq = bx * bx + by * by + bz * bz;
    rgsq = vb2xm * vb2xm + vb2ym * vb2ym + vb2zm * vb2zm;
    rg = sqrt(rgsq);

    rginv = ra2inv = rb2inv = 0.0;
    if (rg > 0) rginv = 1.0 / rg;
    if (rasq > 0) ra2inv = 1.0 / rasq;
    if (rbsq > 0) rb2inv = 1.0 / rbsq;
    rabinv = sqrt(ra2inv * rb2inv);

    c = (ax * bx + ay * by + az * bz) * rabinv;
    s = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE))
      error->one(FLERR, "Dihedral problem: {} {} {} {}", i1 + 1, i2 + 1, i3 + 1,
                 i4 + 1);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    m = multiplicity[type];
    p = 1.0;
    ddf1 = df1 = 0.0;

    for (i = 0; i < m; i++) {
      ddf1 = p * c - df1 * s;
      df1 = p * s + df1 * c;
      p = ddf1;
    }

    p = p * cos_shift[type] + df1 * sin_shift[type];
    df1 = df1 * cos_shift[type] - ddf1 * sin_shift[type];
    df1 *= -m;
    p += 1.0;

    if (m == 0) {
      p = 1.0 + cos_shift[type];
      df1 = 0.0;
    }

    if (eflag) edihedral = w * k[type] * p;

    fg = vb1x * vb2xm + vb1y * vb2ym + vb1z * vb2zm;
    hg = vb3x * vb2xm + vb3y * vb2ym + vb3z * vb2zm;
    fga = fg * ra2inv * rginv;
    hgb = hg * rb2inv * rginv;
    gaa = -ra2inv * rg;
    gbb = rb2inv * rg;

    dtfx = gaa * ax;
    dtfy = gaa * ay;
    dtfz = gaa * az;
    dtgx = fga * ax - hgb * bx;
    dtgy = fga * ay - hgb * by;
    dtgz = fga * az - hgb * bz;
    dthx = gbb * bx;
    dthy = gbb * by;
    dthz = gbb * bz;

    df = -w * k[type] * df1;

    sx2 = df * dtgx;
    sy2 = df * dtgy;
    sz2 = df * dtgz;

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

void DihedralBackmapHarmonic::write_restart(FILE *fp) {
  fwrite(k + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(sign + 1, sizeof(int), atom->ndihedraltypes, fp);
  fwrite(multiplicity + 1, sizeof(int), atom->ndihedraltypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->ndihedraltypes, fp);
}

void DihedralBackmapHarmonic::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, k + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, sign + 1, sizeof(int), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, multiplicity + 1, sizeof(int), atom->ndihedraltypes,
                  fp, nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->ndihedraltypes, fp,
                  nullptr, error);
  }
  MPI_Bcast(k + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(sign + 1, atom->ndihedraltypes, MPI_INT, 0, world);
  MPI_Bcast(multiplicity + 1, atom->ndihedraltypes, MPI_INT, 0, world);
  MPI_Bcast(is_cg + 1, atom->ndihedraltypes, MPI_INT, 0, world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    setflag[i] = 1;
    if (sign[i] == 1) {
      cos_shift[i] = 1;
      sin_shift[i] = 0;
    } else {
      cos_shift[i] = -1;
      sin_shift[i] = 0;
    }
  }
}
