/* ----------------------------------------------------------------------
   dihedral_style backmap/ryckaert — lambda-weighted RB dihedral.

   E = w × Σ Cn cos^n(φ)
   Weight w uses lambda of endpoints i and l.

   Syntax:
     dihedral_style backmap/ryckaert
     dihedral_coeff N at/cg C0 C1 C2 C3 C4 C5
------------------------------------------------------------------------- */

#include "dihedral_backmap_ryckaert.h"

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
static constexpr double SMALL = 0.001;
static constexpr double SMALLER = 0.00001;

/* ---------------------------------------------------------------------- */

DihedralBackmapRyckaert::DihedralBackmapRyckaert(LAMMPS *lmp)
    : Dihedral(lmp),
      c0(nullptr),
      c1(nullptr),
      c2(nullptr),
      c3(nullptr),
      c4(nullptr),
      c5(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {
  writedata = 1;
}

/* ---------------------------------------------------------------------- */

DihedralBackmapRyckaert::~DihedralBackmapRyckaert() {
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(c0);
    memory->destroy(c1);
    memory->destroy(c2);
    memory->destroy(c3);
    memory->destroy(c4);
    memory->destroy(c5);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapRyckaert::allocate() {
  allocated = 1;
  int n = atom->ndihedraltypes + 1;

  memory->create(setflag, n, "dihedral:setflag");
  memory->create(c0, n, "dihedral:c0");
  memory->create(c1, n, "dihedral:c1");
  memory->create(c2, n, "dihedral:c2");
  memory->create(c3, n, "dihedral:c3");
  memory->create(c4, n, "dihedral:c4");
  memory->create(c5, n, "dihedral:c5");
  memory->create(is_cg, n, "dihedral:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapRyckaert::coeff(int narg, char **arg) {
  if (narg != 8)
    error->all(FLERR, "Incorrect args for dihedral_coeff backmap/ryckaert");
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
               "dihedral_coeff backmap/ryckaert: 2nd arg must be 'at' or 'cg'");

  double c0_one = utils::numeric(FLERR, arg[2], false, lmp);
  double c1_one = utils::numeric(FLERR, arg[3], false, lmp);
  double c2_one = utils::numeric(FLERR, arg[4], false, lmp);
  double c3_one = utils::numeric(FLERR, arg[5], false, lmp);
  double c4_one = utils::numeric(FLERR, arg[6], false, lmp);
  double c5_one = utils::numeric(FLERR, arg[7], false, lmp);

  for (int i = ilo; i <= ihi; i++) {
    c0[i] = c0_one;
    c1[i] = c1_one;
    c2[i] = c2_one;
    c3[i] = c3_one;
    c4[i] = c4_one;
    c5[i] = c5_one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapRyckaert::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "dihedral_style backmap/ryckaert");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapRyckaert::compute(int eflag, int vflag) {
  int i1, i2, i3, i4, n, type;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z, vb2xm, vb2ym,
      vb2zm;
  double edihedral, f1[3], f2[3], f3[3], f4[3];
  double sb1, sb2, sb3, rb1, rb3, c0val, b1mag2, b1mag, b2mag2;
  double b2mag, b3mag2, b3mag, ctmp, r12c1, c1mag, r12c2;
  double c2mag, sc1, sc2, s1, s12, c, p, pd, a, a11, a22;
  double a33, a12, a13, a23, sx2, sy2, sz2;
  double s2, cx, cy, cz, cmag, dx, phi, si, siinv;

  edihedral = 0.0;
  ev_init(eflag, vflag);

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "dihedral_style backmap/ryckaert: cannot extract atom2cg/"
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

    sb1 = 1.0 / (vb1x * vb1x + vb1y * vb1y + vb1z * vb1z);
    sb2 = 1.0 / (vb2x * vb2x + vb2y * vb2y + vb2z * vb2z);
    sb3 = 1.0 / (vb3x * vb3x + vb3y * vb3y + vb3z * vb3z);

    rb1 = sqrt(sb1);
    rb3 = sqrt(sb3);

    c0val = (vb1x * vb3x + vb1y * vb3y + vb1z * vb3z) * rb1 * rb3;

    b1mag2 = vb1x * vb1x + vb1y * vb1y + vb1z * vb1z;
    b1mag = sqrt(b1mag2);
    b2mag2 = vb2x * vb2x + vb2y * vb2y + vb2z * vb2z;
    b2mag = sqrt(b2mag2);
    b3mag2 = vb3x * vb3x + vb3y * vb3y + vb3z * vb3z;
    b3mag = sqrt(b3mag2);

    ctmp = vb1x * vb2x + vb1y * vb2y + vb1z * vb2z;
    r12c1 = 1.0 / (b1mag * b2mag);
    c1mag = ctmp * r12c1;

    ctmp = vb2xm * vb3x + vb2ym * vb3y + vb2zm * vb3z;
    r12c2 = 1.0 / (b2mag * b3mag);
    c2mag = ctmp * r12c2;

    sc1 = sqrt(MAX(1.0 - c1mag * c1mag, 0.0));
    if (sc1 < SMALL) sc1 = SMALL;
    sc1 = 1.0 / sc1;

    sc2 = sqrt(MAX(1.0 - c2mag * c2mag, 0.0));
    if (sc2 < SMALL) sc2 = SMALL;
    sc2 = 1.0 / sc2;

    s1 = sc1 * sc1;
    s2 = sc2 * sc2;
    s12 = sc1 * sc2;
    c = (c0val + c1mag * c2mag) * s12;

    cx = vb1y * vb2z - vb1z * vb2y;
    cy = vb1z * vb2x - vb1x * vb2z;
    cz = vb1x * vb2y - vb1y * vb2x;
    cmag = sqrt(cx * cx + cy * cy + cz * cz);
    dx = (cx * vb3x + cy * vb3y + cz * vb3z) / cmag / b3mag;

    if (c > 1.0 + TOLERANCE || c < (-1.0 - TOLERANCE))
      error->one(FLERR,
                 "Dihedral backmap/ryckaert problem: i={} j={} k={} l={}", i1,
                 i2, i3, i4);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    phi = acos(c);
    if (dx < 0.0) phi *= -1.0;
    si = sin(phi);
    if (fabs(si) < SMALLER) si = SMALLER;
    siinv = 1.0 / si;

    double cc = c;
    double cc2 = cc * cc;
    double cc3 = cc2 * cc;
    double cc4 = cc3 * cc;
    double cc5 = cc4 * cc;

    p = c0[type] + c1[type] * cc + c2[type] * cc2 + c3[type] * cc3 +
        c4[type] * cc4 + c5[type] * cc5;
    pd = c1[type] + 2.0 * c2[type] * cc + 3.0 * c3[type] * cc2 +
         4.0 * c4[type] * cc3 + 5.0 * c5[type] * cc4;

    p *= w;
    pd *= w;

    if (eflag) edihedral = p;

    a = pd;
    c = c * a;
    s12 = s12 * a;
    a11 = c * sb1 * s1;
    a22 = -sb2 * (2.0 * c0val * s12 - c * (s1 + s2));
    a33 = c * sb3 * s2;
    a12 = -r12c1 * (c1mag * c * s1 + c2mag * s12);
    a13 = -rb1 * rb3 * s12;
    a23 = r12c2 * (c2mag * c * s2 + c1mag * s12);

    sx2 = a12 * vb1x + a22 * vb2x + a23 * vb3x;
    sy2 = a12 * vb1y + a22 * vb2y + a23 * vb3y;
    sz2 = a12 * vb1z + a22 * vb2z + a23 * vb3z;

    f1[0] = a11 * vb1x + a12 * vb2x + a13 * vb3x;
    f1[1] = a11 * vb1y + a12 * vb2y + a13 * vb3y;
    f1[2] = a11 * vb1z + a12 * vb2z + a13 * vb3z;

    f2[0] = -sx2 - f1[0];
    f2[1] = -sy2 - f1[1];
    f2[2] = -sz2 - f1[2];

    f4[0] = a13 * vb1x + a23 * vb2x + a33 * vb3x;
    f4[1] = a13 * vb1y + a23 * vb2y + a33 * vb3y;
    f4[2] = a13 * vb1z + a23 * vb2z + a33 * vb3z;

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
      ev_tally(i1, i2, i3, i4, nlocal, newton_bond, edihedral, f1, f3, f4, vb1x,
               vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapRyckaert::write_restart(FILE *fp) {
  fwrite(c0 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(c1 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(c2 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(c3 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(c4 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(c5 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->ndihedraltypes, fp);
}

void DihedralBackmapRyckaert::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, c0 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, c1 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, c2 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, c3 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, c4 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, c5 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->ndihedraltypes, fp,
                  nullptr, error);
  }
  MPI_Bcast(c0 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(c1 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(c2 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(c3 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(c4 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(c5 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(is_cg + 1, atom->ndihedraltypes, MPI_INT, 0, world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) setflag[i] = 1;
}
