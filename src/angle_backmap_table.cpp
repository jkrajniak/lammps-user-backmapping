/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* angle_style backmap/table — lambda-weighted tabulated angle for backmapping.

   F = w × F_table(θ)
   E = w × E_table(θ)

   Syntax:
     angle_style backmap/table linear N
     angle_coeff M at/cg filename keyword */

#include "angle_backmap_table.h"

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
#include "text_file_reader.h"
#include "tokenizer.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */

AngleBackmapTable::AngleBackmapTable(LAMMPS *lmp)
    : Angle(lmp),
      tabstyle(0),
      tablength(0),
      ntables(0),
      tables(nullptr),
      tabindex(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

AngleBackmapTable::~AngleBackmapTable() {
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::allocate() {
  allocated = 1;
  int n = atom->nangletypes + 1;

  memory->create(setflag, n, "angle:setflag");
  memory->create(tabindex, n, "angle:tabindex");
  memory->create(is_cg, n, "angle:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::settings(int narg, char **arg) {
  if (narg != 2) error->all(FLERR, "Illegal angle_style backmap/table command");

  if (strcmp(arg[0], "linear") == 0)
    tabstyle = 0;
  else if (strcmp(arg[0], "spline") == 0)
    tabstyle = 1;
  else
    error->all(FLERR, "Unknown table style in angle_style backmap/table");

  tablength = utils::inumeric(FLERR, arg[1], false, lmp);
  if (tablength < 2)
    error->all(FLERR, "Illegal table length in angle_style backmap/table");
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for angle_coeff backmap/table");
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
               "angle_coeff backmap/table: 2nd arg must be 'at' or 'cg'");

  tables = static_cast<Table *>(
      memory->srealloc(tables, (ntables + 1) * sizeof(Table), "angle:tables"));
  Table *tb = &tables[ntables];
  null_table(tb);

  if (comm->me == 0) read_table(tb, arg[2], arg[3]);
  bcast_table(tb);

  tb->delta = (tb->hi - tb->lo) / (tablength - 1);
  tb->invdelta = 1.0 / tb->delta;

  if (tabstyle == 1) spline_table(tb);
  compute_table(tb);

  for (int i = ilo; i <= ihi; i++) {
    tabindex[i] = ntables;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
  ntables++;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "angle_style backmap/table");
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "angle_style backmap/table: cannot extract atom2cg/"
               "lambda_global");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

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

    bool same_bead = BackmapLambda::same_bead(atom2cg, i1, i2, i3);
    double w =
        BackmapLambda::compute_weight3(same_bead, is_cg[atype], lambda_global);

    if (BackmapLambda::is_almost_zero(w)) continue;

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

    double c = (delx1 * delx2 + dely1 * dely2 + delz1 * delz2) / (r1 * r2);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    double s = sqrt(1.0 - c * c);
    if (s < 0.001) s = 0.001;
    s = 1.0 / s;

    double theta_deg = acos(c) * RAD2DEG;
    double fforce;
    double eangle = uf_lookup(tabindex[atype], theta_deg, fforce);

    eangle *= w;
    fforce *= w;
    fforce *= -1.0;

    double a = -fforce * s;
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

    if (evflag)
      ev_tally(i1, i2, i3, nlocal, newton_bond, eangle, f1, f3, delx1, dely1,
               delz1, delx2, dely2, delz2);
  }
}

/* ---------------------------------------------------------------------- */

double AngleBackmapTable::equilibrium_angle(int atype) {
  Table *tb = &tables[tabindex[atype]];
  return 0.5 * (tb->lo + tb->hi) * MY_PI / 180.0;
}

/* ---------------------------------------------------------------------- */

double AngleBackmapTable::single(int atype, int i1, int i2, int i3) {
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
  double theta_deg = acos(c) * RAD2DEG;

  double fforce;
  double eng = uf_lookup(tabindex[atype], theta_deg, fforce);

  double w = 1.0;
  if (fix_backmap) {
    int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
    double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
    if (atom2cg && lam_global_ptr) {
      double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);
      bool same_bead = BackmapLambda::same_bead(atom2cg, i1, i2, i3);
      w = BackmapLambda::compute_weight3(same_bead, is_cg[atype],
                                         lambda_global);
    }
  }
  return eng * w;
}

/* ---------------------------------------------------------------------- */

double AngleBackmapTable::uf_lookup(int tindex, double theta, double &fforce) {
  Table *tb = &tables[tindex];

  if (theta < tb->lo) theta = tb->lo;
  if (theta > tb->hi) theta = tb->hi;

  double fraction = (theta - tb->lo) * tb->invdelta;
  int itable = static_cast<int>(fraction);
  if (itable >= tablength - 1) itable = tablength - 2;
  if (itable < 0) itable = 0;
  fraction -= itable;

  double eng = tb->e[itable] + fraction * tb->de[itable];
  fforce = tb->f[itable] + fraction * tb->df[itable];

  return eng;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::null_table(Table *tb) {
  tb->ninput = 0;
  tb->lo = tb->hi = 0.0;
  tb->rfile = tb->efile = tb->ffile = nullptr;
  tb->e2file = tb->f2file = nullptr;
  tb->delta = tb->invdelta = 0.0;
  tb->r = tb->e = tb->f = nullptr;
  tb->de = tb->df = nullptr;
}

void AngleBackmapTable::free_table(Table *tb) {
  memory->destroy(tb->rfile);
  memory->destroy(tb->efile);
  memory->destroy(tb->ffile);
  memory->destroy(tb->e2file);
  memory->destroy(tb->f2file);
  memory->destroy(tb->r);
  memory->destroy(tb->e);
  memory->destroy(tb->f);
  memory->destroy(tb->de);
  memory->destroy(tb->df);
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::read_table(Table *tb, const char *file,
                                   const char *keyword) {
  std::string filecontent = utils::get_potential_file_path(file);
  if (filecontent.empty())
    error->one(FLERR, "Cannot open angle table file {}", file);

  auto reader = TextFileReader(filecontent, "angle table");
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

  memory->create(tb->rfile, tb->ninput, "angle:rfile");
  memory->create(tb->efile, tb->ninput, "angle:efile");
  memory->create(tb->ffile, tb->ninput, "angle:ffile");

  for (int i = 0; i < tb->ninput; i++) {
    auto line = reader.next_line();
    if (!line) error->one(FLERR, "Premature end of angle table file {}", file);
    ValueTokenizer values(line);
    values.next_int();
    tb->rfile[i] = values.next_double();
    tb->efile[i] = values.next_double();
    tb->ffile[i] = values.next_double();
  }

  tb->lo = tb->rfile[0];
  tb->hi = tb->rfile[tb->ninput - 1];
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::bcast_table(Table *tb) {
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
  MPI_Bcast(&tb->lo, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb->hi, 1, MPI_DOUBLE, 0, world);

  if (comm->me != 0) {
    memory->create(tb->rfile, tb->ninput, "angle:rfile");
    memory->create(tb->efile, tb->ninput, "angle:efile");
    memory->create(tb->ffile, tb->ninput, "angle:ffile");
  }

  MPI_Bcast(tb->rfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->efile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::spline_table(Table *tb) {
  memory->create(tb->e2file, tb->ninput, "angle:e2file");
  memory->create(tb->f2file, tb->ninput, "angle:f2file");

  auto *u = new double[tb->ninput];

  tb->e2file[0] = 0.0;
  u[0] = 0.0;
  for (int i = 1; i < tb->ninput - 1; i++) {
    double h1 = tb->rfile[i] - tb->rfile[i - 1];
    double h2 = tb->rfile[i + 1] - tb->rfile[i];
    double sig = h1 / (h1 + h2);
    double p = sig * tb->e2file[i - 1] + 2.0;
    tb->e2file[i] = (sig - 1.0) / p;
    u[i] = (tb->efile[i + 1] - tb->efile[i]) / h2 -
           (tb->efile[i] - tb->efile[i - 1]) / h1;
    u[i] = (6.0 * u[i] / (h1 + h2) - sig * u[i - 1]) / p;
  }
  tb->e2file[tb->ninput - 1] = 0.0;
  for (int i = tb->ninput - 2; i >= 0; i--)
    tb->e2file[i] = tb->e2file[i] * tb->e2file[i + 1] + u[i];

  tb->f2file[0] = 0.0;
  u[0] = 0.0;
  for (int i = 1; i < tb->ninput - 1; i++) {
    double h1 = tb->rfile[i] - tb->rfile[i - 1];
    double h2 = tb->rfile[i + 1] - tb->rfile[i];
    double sig = h1 / (h1 + h2);
    double p = sig * tb->f2file[i - 1] + 2.0;
    tb->f2file[i] = (sig - 1.0) / p;
    u[i] = (tb->ffile[i + 1] - tb->ffile[i]) / h2 -
           (tb->ffile[i] - tb->ffile[i - 1]) / h1;
    u[i] = (6.0 * u[i] / (h1 + h2) - sig * u[i - 1]) / p;
  }
  tb->f2file[tb->ninput - 1] = 0.0;
  for (int i = tb->ninput - 2; i >= 0; i--)
    tb->f2file[i] = tb->f2file[i] * tb->f2file[i + 1] + u[i];

  delete[] u;
}

/* ---------------------------------------------------------------------- */

void AngleBackmapTable::compute_table(Table *tb) {
  memory->create(tb->r, tablength, "angle:r");
  memory->create(tb->e, tablength, "angle:e");
  memory->create(tb->f, tablength, "angle:f");
  memory->create(tb->de, tablength, "angle:de");
  memory->create(tb->df, tablength, "angle:df");

  for (int i = 0; i < tablength; i++) {
    tb->r[i] = tb->lo + i * tb->delta;

    double fraction =
        (tb->r[i] - tb->rfile[0]) / (tb->rfile[tb->ninput - 1] - tb->rfile[0]);
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

void AngleBackmapTable::write_restart(FILE *fp) {
  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nangletypes, fp);
}

void AngleBackmapTable::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->nangletypes, fp, nullptr,
                  error);
  }
  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);
  MPI_Bcast(is_cg + 1, atom->nangletypes, MPI_INT, 0, world);
}
