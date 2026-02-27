/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* bond_style backmap/table — lambda-weighted tabulated bond for backmapping.

   F = w × F_table(r)
   E = w × E_table(r)

   Syntax:
     bond_style backmap/table linear N
     bond_coeff M at/cg filename keyword */

#include "bond_backmap_table.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "tokenizer.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

static constexpr int MAXTABLE = 10000;

/* ---------------------------------------------------------------------- */

BondBackmapTable::BondBackmapTable(LAMMPS *lmp)
    : Bond(lmp),
      tabstyle(0),
      tablength(0),
      ntables(0),
      tables(nullptr),
      tabindex(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

BondBackmapTable::~BondBackmapTable() {
  for (int m = 0; m < ntables; m++) free_table(&tables[m]);
  memory->sfree(tables);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(tabindex);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::allocate() {
  allocated = 1;
  int n = atom->nbondtypes + 1;

  memory->create(setflag, n, "bond:setflag");
  memory->create(tabindex, n, "bond:tabindex");
  memory->create(is_cg, n, "bond:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::settings(int narg, char **arg) {
  if (narg != 2) error->all(FLERR, "Illegal bond_style backmap/table command");

  if (strcmp(arg[0], "linear") == 0)
    tabstyle = 0;
  else if (strcmp(arg[0], "spline") == 0)
    tabstyle = 1;
  else
    error->all(FLERR, "Unknown table style in bond_style backmap/table");

  tablength = utils::inumeric(FLERR, arg[1], false, lmp);
  if (tablength < 2)
    error->all(FLERR, "Illegal table length in bond_style backmap/table");
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for bond_coeff backmap/table");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  int cg_flag;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR, "bond_coeff backmap/table: 2nd arg must be 'at' or 'cg'");

  // Allocate new table
  tables = static_cast<Table *>(
      memory->srealloc(tables, (ntables + 1) * sizeof(Table), "bond:tables"));
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

void BondBackmapTable::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "bond_style backmap/table");
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  double *lam = BackmapLambda::extract_lambda(fix_backmap);
  if (!lam)
    error->all(FLERR, "bond_style backmap/table: cannot extract lambda");

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  for (int n = 0; n < nbondlist; n++) {
    int i1 = bondlist[n][0];
    int i2 = bondlist[n][1];
    int btype = bondlist[n][2];

    double delx = x[i1][0] - x[i2][0];
    double dely = x[i1][1] - x[i2][1];
    double delz = x[i1][2] - x[i2][2];

    double rsq = delx * delx + dely * dely + delz * delz;
    double r = sqrt(rsq);

    // Lambda weighting
    double li = BackmapLambda::clamp_lambda(lam[i1]);
    double lj = BackmapLambda::clamp_lambda(lam[i2]);
    double w = BackmapLambda::compute_weight(li, lj, is_cg[btype]);

    if (BackmapLambda::is_almost_zero(w)) continue;

    double fforce;
    double ebond = uf_lookup(tabindex[btype], r, fforce);

    fforce *= w;
    ebond *= w;

    double fbond;
    if (r > 0.0)
      fbond = fforce / r;
    else
      fbond = 0.0;

    f[i1][0] += delx * fbond;
    f[i1][1] += dely * fbond;
    f[i1][2] += delz * fbond;
    f[i2][0] -= delx * fbond;
    f[i2][1] -= dely * fbond;
    f[i2][2] -= delz * fbond;

    if (evflag)
      ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

double BondBackmapTable::equilibrium_distance(int btype) {
  Table *tb = &tables[tabindex[btype]];
  return 0.5 * (tb->lo + tb->hi);
}

/* ---------------------------------------------------------------------- */

double BondBackmapTable::single(int btype, double rsq, int i, int j,
                                double &fforce) {
  double r = sqrt(rsq);
  double eng = uf_lookup(tabindex[btype], r, fforce);

  double w = 1.0;
  if (fix_backmap) {
    double *lam = BackmapLambda::extract_lambda(fix_backmap);
    if (lam) {
      double li = BackmapLambda::clamp_lambda(lam[i]);
      double lj_val = BackmapLambda::clamp_lambda(lam[j]);
      w = BackmapLambda::compute_weight(li, lj_val, is_cg[btype]);
    }
  }
  fforce *= w;
  return eng * w;
}

/* ---------------------------------------------------------------------- */

double BondBackmapTable::uf_lookup(int tindex, double r, double &fforce) {
  Table *tb = &tables[tindex];

  if (r < tb->lo || r > tb->hi)
    error->one(FLERR, "Bond distance {} out of range [{},{}] in backmap/table",
               r, tb->lo, tb->hi);

  double fraction = (r - tb->lo) * tb->invdelta;
  int itable = static_cast<int>(fraction);
  if (itable >= tablength - 1) itable = tablength - 2;
  fraction -= itable;

  double eng = tb->e[itable] + fraction * tb->de[itable];
  fforce = tb->f[itable] + fraction * tb->df[itable];

  return eng;
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::null_table(Table *tb) {
  tb->ninput = 0;
  tb->lo = tb->hi = 0.0;
  tb->rfile = tb->efile = tb->ffile = nullptr;
  tb->e2file = tb->f2file = nullptr;
  tb->delta = tb->invdelta = 0.0;
  tb->r = tb->e = tb->f = nullptr;
  tb->de = tb->df = nullptr;
}

void BondBackmapTable::free_table(Table *tb) {
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

void BondBackmapTable::read_table(Table *tb, const char *file,
                                  const char *keyword) {
  std::string filecontent = utils::get_potential_file_path(file);
  if (filecontent.empty())
    error->one(FLERR, "Cannot open bond table file {}", file);

  auto reader = TextFileReader(filecontent, "bond table");
  reader.ignore_comments = true;

  // Find keyword section
  while (true) {
    auto line = reader.next_line();
    if (!line)
      error->one(FLERR, "Did not find keyword {} in table file {}", keyword,
                 file);
    ValueTokenizer values(line);
    std::string word = values.next_string();
    if (word == keyword) {
      if (values.has_next()) tb->ninput = values.next_int();
      break;
    }
  }

  if (tb->ninput <= 1)
    error->one(FLERR, "Invalid table length in file {}", file);

  memory->create(tb->rfile, tb->ninput, "bond:rfile");
  memory->create(tb->efile, tb->ninput, "bond:efile");
  memory->create(tb->ffile, tb->ninput, "bond:ffile");

  for (int i = 0; i < tb->ninput; i++) {
    auto line = reader.next_line();
    if (!line) error->one(FLERR, "Premature end of bond table file {}", file);
    ValueTokenizer values(line);
    values.next_int();  // index
    tb->rfile[i] = values.next_double();
    tb->efile[i] = values.next_double();
    tb->ffile[i] = values.next_double();
  }

  tb->lo = tb->rfile[0];
  tb->hi = tb->rfile[tb->ninput - 1];
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::bcast_table(Table *tb) {
  MPI_Bcast(&tb->ninput, 1, MPI_INT, 0, world);
  MPI_Bcast(&tb->lo, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&tb->hi, 1, MPI_DOUBLE, 0, world);

  if (comm->me != 0) {
    memory->create(tb->rfile, tb->ninput, "bond:rfile");
    memory->create(tb->efile, tb->ninput, "bond:efile");
    memory->create(tb->ffile, tb->ninput, "bond:ffile");
  }

  MPI_Bcast(tb->rfile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->efile, tb->ninput, MPI_DOUBLE, 0, world);
  MPI_Bcast(tb->ffile, tb->ninput, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

void BondBackmapTable::spline_table(Table *tb) {
  memory->create(tb->e2file, tb->ninput, "bond:e2file");
  memory->create(tb->f2file, tb->ninput, "bond:f2file");

  // Natural cubic spline (second derivatives at endpoints = 0)
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

void BondBackmapTable::compute_table(Table *tb) {
  memory->create(tb->r, tablength, "bond:r");
  memory->create(tb->e, tablength, "bond:e");
  memory->create(tb->f, tablength, "bond:f");
  memory->create(tb->de, tablength, "bond:de");
  memory->create(tb->df, tablength, "bond:df");

  for (int i = 0; i < tablength; i++) {
    tb->r[i] = tb->lo + i * tb->delta;

    // Linear interpolation from input data
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

void BondBackmapTable::write_restart(FILE *fp) {
  fwrite(&tabstyle, sizeof(int), 1, fp);
  fwrite(&tablength, sizeof(int), 1, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nbondtypes, fp);
}

void BondBackmapTable::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, &tabstyle, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, &tablength, sizeof(int), 1, fp, nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->nbondtypes, fp, nullptr,
                  error);
  }
  MPI_Bcast(&tabstyle, 1, MPI_INT, 0, world);
  MPI_Bcast(&tablength, 1, MPI_INT, 0, world);
  MPI_Bcast(is_cg + 1, atom->nbondtypes, MPI_INT, 0, world);
}
