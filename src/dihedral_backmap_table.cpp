/* ----------------------------------------------------------------------
   dihedral_style backmap/table — lambda-weighted tabulated dihedral.

   F = w × F_table(φ)
   E = w × E_table(φ)

   Syntax:
     dihedral_style backmap/table linear N
     dihedral_coeff M at/cg filename keyword
------------------------------------------------------------------------- */

#include "dihedral_backmap_table.h"

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
#include "text_file_reader.h"
#include "tokenizer.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

static constexpr double TOLERANCE = 0.05;
static constexpr double SMALL = 0.001;

/* ---------------------------------------------------------------------- */

DihedralBackmapTable::DihedralBackmapTable(LAMMPS *lmp)
    : Dihedral(lmp),
      tabstyle(0),
      tablength(0),
      ntables(0),
      tables(nullptr),
      tabindex(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

DihedralBackmapTable::~DihedralBackmapTable() {
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::allocate() {
  allocated = 1;
  int n = atom->ndihedraltypes + 1;

  memory->create(setflag, n, "dihedral:setflag");
  memory->create(tabindex, n, "dihedral:tabindex");
  memory->create(is_cg, n, "dihedral:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::settings(int narg, char **arg) {
  if (narg != 2)
    error->all(FLERR, "Illegal dihedral_style backmap/table command");

  if (strcmp(arg[0], "linear") == 0)
    tabstyle = 0;
  else if (strcmp(arg[0], "spline") == 0)
    tabstyle = 1;
  else
    error->all(FLERR, "Unknown table style in dihedral_style backmap/table");

  tablength = utils::inumeric(FLERR, arg[1], false, lmp);
  if (tablength < 2)
    error->all(FLERR, "Illegal table length in dihedral_style backmap/table");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for dihedral_coeff backmap/table");
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
               "dihedral_coeff backmap/table: 2nd arg must be 'at' or 'cg'");

  tables = static_cast<Table *>(memory->srealloc(
      tables, (ntables + 1) * sizeof(Table), "dihedral:tables"));
  Table *tb = &tables[ntables];
  null_table(tb);

  if (comm->me == 0) read_table(tb, arg[2], arg[3]);
  bcast_table(tb);

  tb->delta = (tb->hi - tb->lo) / (tablength - 1);
  tb->invdelta = 1.0 / tb->delta;

  compute_table(tb);

  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
  ntables++;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "dihedral_style backmap/table");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::compute(int eflag, int vflag) {
  int i1, i2, i3, i4, n, type;
  double vb1x, vb1y, vb1z, vb2x, vb2y, vb2z, vb3x, vb3y, vb3z, vb2xm, vb2ym,
      vb2zm;
  double edihedral, f1[3], f2[3], f3[3], f4[3];
  double sb1, sb2, sb3, rb1, rb3, c0val, b1mag2, b1mag, b2mag2;
  double b2mag, b3mag2, b3mag, ctmp, r12c1, c1mag, r12c2;
  double c2mag, sc1, sc2, s1, s12, c, pd, a, a11, a22;
  double a33, a12, a13, a23, sx2, sy2, sz2;
  double s2, cx, cy, cz, cmag, dx, phi;

  edihedral = 0.0;
  ev_init(eflag, vflag);

  double *lam = BackmapLambda::extract_lambda(fix_backmap);
  if (!lam)
    error->all(FLERR, "dihedral_style backmap/table: cannot extract lambda");

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

    double li = BackmapLambda::clamp_lambda(lam[i1]);
    double ll = BackmapLambda::clamp_lambda(lam[i4]);
    double w = BackmapLambda::compute_weight(li, ll, is_cg[type]);

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
      error->one(FLERR, "Dihedral backmap/table problem: i={} j={} k={} l={}",
                 i1, i2, i3, i4);

    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;

    phi = acos(c);
    if (dx < 0.0) phi *= -1.0;

    double phi_deg = phi * RAD2DEG;
    double fforce;
    double p = uf_lookup(tabindex[type], phi_deg, fforce);

    p *= w;
    pd = -fforce * (MY_PI / 180.0) * w;

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

double DihedralBackmapTable::uf_lookup(int tindex, double phi_deg,
                                       double &fforce) {
  Table *tb = &tables[tindex];

  if (phi_deg < tb->lo) phi_deg = tb->lo;
  if (phi_deg > tb->hi) phi_deg = tb->hi;

  double fraction = (phi_deg - tb->lo) * tb->invdelta;
  int itable = static_cast<int>(fraction);
  if (itable >= tablength - 1) itable = tablength - 2;
  if (itable < 0) itable = 0;
  fraction -= itable;

  double eng = tb->e[itable] + fraction * tb->de[itable];
  fforce = tb->f[itable] + fraction * tb->df[itable];

  return eng;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::null_table(Table *tb) {
  tb->ninput = 0;
  tb->lo = tb->hi = 0.0;
  tb->phifile = tb->efile = tb->ffile = nullptr;
  tb->delta = tb->invdelta = 0.0;
  tb->phi = tb->e = tb->f = nullptr;
  tb->de = tb->df = nullptr;
}

void DihedralBackmapTable::free_table(Table *tb) {
  memory->destroy(tb->phifile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->phi);
  memory->destroy(tb->e);
  memory->destroy(tb->f);
  memory->destroy(tb->de);
  memory->destroy(tb->df);
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::read_table(Table *tb, const char *file,
                                      const char *keyword) {
  std::string filecontent = utils::get_potential_file_path(file);
  if (filecontent.empty())
    error->one(FLERR, "Cannot open dihedral table file {}", file);

  auto reader = TextFileReader(filecontent, "dihedral table");
  reader.ignore_comments = true;

  while (true) {
    auto line = reader.next_line();
    if (!line)
      error->one(FLERR, "Did not find keyword {} in table file {}", keyword,
                 file);
    ValueTokenizer values(line);
    std::string word = values.next_string();
    if (word == keyword) break;
  }

  {
    auto line = reader.next_line();
    if (!line)
      error->one(FLERR, "Missing parameter line after keyword {} in {}",
                 keyword, file);
    ValueTokenizer values(line);
    std::string key = values.next_string();
    if (key == "N" && values.has_next()) tb->ninput = values.next_int();
  }

  if (tb->ninput <= 1)
    error->one(FLERR, "Invalid table length in file {}", file);

  memory->create(tb->phifile, tb->ninput, "dihedral:phifile");
  memory->create(tb->efile, tb->ninput, "dihedral:efile");
  memory->create(tb->ffile, tb->ninput, "dihedral:ffile");

  for (int i = 0; i < tb->ninput; i++) {
    auto line = reader.next_line();
    if (!line)
      error->one(FLERR, "Premature end of dihedral table file {}", file);
    ValueTokenizer values(line);
    values.next_int();
    tb->phifile[i] = values.next_double();
    tb->efile[i] = values.next_double();
    tb->ffile[i] = values.next_double();
  }

  tb->lo = tb->phifile[0];
  tb->hi = tb->phifile[tb->ninput - 1];
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::bcast_table(Table *tb) {
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
  MPI_Bcast(&tb->lo, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb->hi, 1, MPI_DOUBLE, 0, world);

  if (comm->me != 0) {
    memory->create(tb->phifile, tb->ninput, "dihedral:phifile");
    memory->create(tb->efile, tb->ninput, "dihedral:efile");
    memory->create(tb->ffile, tb->ninput, "dihedral:ffile");
  }

  MPI_Bcast(tb->phifile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->efile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::compute_table(Table *tb) {
  memory->create(tb->phi, tablength, "dihedral:phi");
  memory->create(tb->e, tablength, "dihedral:e");
  memory->create(tb->f, tablength, "dihedral:f");
  memory->create(tb->de, tablength, "dihedral:de");
  memory->create(tb->df, tablength, "dihedral:df");

  for (int i = 0; i < tablength; i++) {
    tb->phi[i] = tb->lo + i * tb->delta;

    double span = tb->phifile[tb->ninput - 1] - tb->phifile[0];
    if (span <= 0.0) span = 1.0;
    double fraction = (tb->phi[i] - tb->phifile[0]) / span;
    int idx = static_cast<int>(fraction * (tb->ninput - 1));
    if (idx >= tb->ninput - 1) idx = tb->ninput - 2;
    if (idx < 0) idx = 0;
    double frac = fraction * (tb->ninput - 1) - idx;

    tb->e[i] = tb->efile[idx] + frac * (tb->efile[idx + 1] - tb->efile[idx]);
    tb->f[i] = tb->ffile[idx] + frac * (tb->ffile[idx + 1] - tb->ffile[idx]);
  }

  for (int i = 0; i < tablength - 1; i++) {
    tb->de[i] = tb->e[i + 1] - tb->e[i];
    tb->df[i] = tb->f[i + 1] - tb->f[i];
  }
  tb->de[tablength - 1] = 0.0;
  tb->df[tablength - 1] = 0.0;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapTable::write_restart(FILE *fp) {
  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
  fwrite(is_cg + 1, sizeof(int), atom->ndihedraltypes, fp);
}

void DihedralBackmapTable::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->ndihedraltypes, fp,
                  nullptr, error);
  }
  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);
  MPI_Bcast(is_cg + 1, atom->ndihedraltypes, MPI_INT, 0, world);
}
