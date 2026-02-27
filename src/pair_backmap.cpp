/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* pair_style backmap — lambda-weighted non-bonded pair forces for backmapping.

   Delegates force computation to AT and CG sub-styles, then weights by:
     AT: w = λ_i × λ_j
     CG: w = 1 − λ_i × λ_j

   Syntax:
     pair_style backmap cut_at at_style at_args ... cut_cg cg_style cg_args ...
     pair_coeff I J atomistic at_args ...
     pair_coeff I J cg cg_args ...
     pair_coeff I J none

   Reference: Krajniak et al., JCTC 2016, DOI: 10.1021/acs.jctc.6b00595 */

#include "pair_backmap.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairBackmap::PairBackmap(LAMMPS *lmp)
    : Pair(lmp),
      pair_at(nullptr),
      pair_cg(nullptr),
      style_at(nullptr),
      style_cg(nullptr),
      cut_at(0.0),
      cut_cg(0.0),
      cut_global(0.0),
      cut(nullptr),
      pair_kind(nullptr),
      fix_backmap(nullptr),
      f_at(nullptr),
      f_cg(nullptr),
      eng_vdwl_at(0.0),
      eng_vdwl_cg(0.0),
      eng_coul_at(0.0),
      eng_coul_cg(0.0),
      nmax_force(0) {
  restartinfo = 1;
  writedata = 0;
  allocated = 0;
  memset(virial_at, 0, sizeof(virial_at));
  memset(virial_cg, 0, sizeof(virial_cg));
}

/* ---------------------------------------------------------------------- */

PairBackmap::~PairBackmap() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(pair_kind);
  }
  delete[] style_at;
  delete[] style_cg;
  delete pair_at;
  delete pair_cg;
  memory->destroy(f_at);
  memory->destroy(f_cg);
}

/* ---------------------------------------------------------------------- */

void PairBackmap::allocate() {
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(cut, n, n, "pair:cut");
  memory->create(pair_kind, n, n, "pair:pair_kind");

  for (int i = 1; i < n; i++) {
    for (int j = i; j < n; j++) {
      setflag[i][j] = 0;
      pair_kind[i][j] = NONE;
    }
  }
}

/* ---------------------------------------------------------------------- */

/* C1: settings — parse pair_style command, create sub-styles.

   pair_style backmap cut_at at_style at_args ... cut_cg cg_style cg_args ... */

void PairBackmap::settings(int narg, char **arg) {
  if (narg < 4) error->all(FLERR, "Illegal pair_style backmap command");

  // First argument: AT cutoff
  cut_at = utils::numeric(FLERR, arg[0], false, lmp);
  if (cut_at <= 0.0)
    error->all(FLERR, "pair_style backmap: AT cutoff must be positive");

  // Second argument: AT sub-style name
  style_at = utils::strdup(arg[1]);
  int dummy = 0;
  pair_at = force->new_pair(arg[1], 1, dummy);

  // Find CG cutoff: scan for the next numeric argument that could be a cutoff
  // after the AT sub-style arguments
  int cg_start = -1;
  for (int i = 2; i < narg - 1; i++) {
    char *endptr;
    double val = strtod(arg[i], &endptr);
    if (*endptr == '\0' && val > 0.0) {
      // Could be a cutoff — check if next arg is a pair style name
      if (i + 1 < narg && force->pair_map->count(arg[i + 1])) {
        cg_start = i;
        break;
      }
    }
  }

  if (cg_start < 0)
    error->all(FLERR, "pair_style backmap: cannot find CG cutoff/style");

  // Pass AT sub-style settings (args between style name and CG cutoff)
  int narg_at = cg_start - 2;
  if (narg_at > 0) pair_at->settings(narg_at, &arg[2]);

  // CG cutoff and style
  cut_cg = utils::numeric(FLERR, arg[cg_start], false, lmp);
  if (cut_cg <= 0.0)
    error->all(FLERR, "pair_style backmap: CG cutoff must be positive");

  style_cg = utils::strdup(arg[cg_start + 1]);
  pair_cg = force->new_pair(arg[cg_start + 1], 1, dummy);

  int narg_cg = narg - cg_start - 2;
  if (narg_cg > 0) pair_cg->settings(narg_cg, &arg[cg_start + 2]);

  cut_global = MAX(cut_at, cut_cg);
}

/* ---------------------------------------------------------------------- */

/* C2: coeff — parse pair_coeff commands.

   pair_coeff I J atomistic at_args ...
   pair_coeff I J cg cg_args ...
   pair_coeff I J none */

void PairBackmap::coeff(int narg, char **arg) {
  if (narg < 3) error->all(FLERR, "Incorrect args for pair_coeff backmap");
  if (!allocated) allocate();

  int ilo, ihi, jlo, jhi;
  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);

  if (strcmp(arg[2], "atomistic") == 0) {
    // Forward to AT sub-style
    if (pair_at) pair_at->coeff(narg - 1, &arg[1]);
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo, i); j <= jhi; j++) {
        pair_kind[i][j] = ATOMISTIC;
        setflag[i][j] = 1;
        cut[i][j] = cut_at;
      }
    }
  } else if (strcmp(arg[2], "cg") == 0) {
    // Forward to CG sub-style
    if (pair_cg) pair_cg->coeff(narg - 1, &arg[1]);
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo, i); j <= jhi; j++) {
        pair_kind[i][j] = CG;
        setflag[i][j] = 1;
        cut[i][j] = cut_cg;
      }
    }
  } else if (strcmp(arg[2], "none") == 0) {
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo, i); j <= jhi; j++) {
        pair_kind[i][j] = NONE;
        setflag[i][j] = 1;
        cut[i][j] = 0.0;
      }
    }
  } else {
    error->all(FLERR,
               "pair_coeff backmap: 3rd arg must be atomistic, cg, or none");
  }
}

/* ---------------------------------------------------------------------- */

void PairBackmap::init_style() {
  // Locate fix backmap
  fix_backmap = BackmapLambda::find_fix_backmap(lmp, "pair_style backmap");

  // Initialize sub-styles
  if (pair_at) pair_at->init_style();
  if (pair_cg) pair_cg->init_style();

  // Request a neighbor list
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

double PairBackmap::init_one(int i, int j) {
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

  // Also init sub-style pair parameters for this type pair
  if (pair_kind[i][j] == ATOMISTIC && pair_at)
    pair_at->init_one(i, j);
  else if (pair_kind[i][j] == CG && pair_cg)
    pair_cg->init_one(i, j);

  cut[j][i] = cut[i][j];
  pair_kind[j][i] = pair_kind[i][j];

  return cut[i][j];
}

/* ---------------------------------------------------------------------- */

/* C3 + C4: compute — lambda-weighted pair force and energy computation.

   For each pair:
   - AT type pairs: compute AT force, weight by w_AT = λ_i × λ_j
   - CG type pairs: compute CG force, weight by w_CG = 1 − λ_i × λ_j
   - Cross-type or NONE pairs: skip */

void PairBackmap::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  if (!fix_backmap)
    fix_backmap = BackmapLambda::find_fix_backmap(lmp, "pair_style backmap");

  double *lam = BackmapLambda::extract_lambda(fix_backmap);
  if (!lam)
    error->all(FLERR,
               "pair_style backmap: cannot extract lambda from fix backmap");

  double **x = atom->x;
  double **f = atom->f;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int inum = list->inum;

  for (int ii = 0; ii < inum; ii++) {
    int i = ilist[ii];
    double li = BackmapLambda::clamp_lambda(lam[i]);
    double xi = x[i][0];
    double yi = x[i][1];
    double zi = x[i][2];
    int itype = type[i];

    int *jlist = firstneigh[i];
    int jnum = numneigh[i];

    for (int jj = 0; jj < jnum; jj++) {
      int j = jlist[jj] & NEIGHMASK;
      int jtype = type[j];

      int kind = pair_kind[itype][jtype];
      if (kind == NONE) continue;

      double delx = xi - x[j][0];
      double dely = yi - x[j][1];
      double delz = zi - x[j][2];
      double rsq = delx * delx + dely * dely + delz * delz;

      if (rsq >= cutsq[itype][jtype]) continue;

      double lj = BackmapLambda::clamp_lambda(lam[j]);
      bool is_cg = (kind == CG);
      double w = BackmapLambda::compute_weight(li, lj, is_cg);

      if (BackmapLambda::is_almost_zero(w)) continue;

      // Use the appropriate sub-style's single() to get fpair and eng
      Pair *sub = is_cg ? pair_cg : pair_at;
      double fforce = 0.0;
      double eng = sub->single(i, j, itype, jtype, rsq, 1.0, 1.0, fforce);

      // Scale by lambda weight
      fforce *= w;
      eng *= w;

      // Apply forces
      f[i][0] += delx * fforce;
      f[i][1] += dely * fforce;
      f[i][2] += delz * fforce;
      if (newton_pair || j < nlocal) {
        f[j][0] -= delx * fforce;
        f[j][1] -= dely * fforce;
        f[j][2] -= delz * fforce;
      }

      // Tally energy and virial
      if (eflag) {
        if (eflag_global) eng_vdwl += eng;
        if (eflag_atom) {
          double ehalf = 0.5 * eng;
          if (newton_pair || i < nlocal) eatom[i] += ehalf;
          if (newton_pair || j < nlocal) eatom[j] += ehalf;
        }
      }

      if (evflag) {
        ev_tally(i, j, nlocal, newton_pair, eng, 0.0, fforce, delx, dely, delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

/* C5: single — compute force/energy for a single pair (for diagnostics). */

double PairBackmap::single(int i, int j, int itype, int jtype, double rsq,
                           double factor_coul, double factor_lj,
                           double &fforce) {
  fforce = 0.0;
  int kind = pair_kind[itype][jtype];
  if (kind == NONE) return 0.0;

  Pair *sub = (kind == CG) ? pair_cg : pair_at;
  double eng =
      sub->single(i, j, itype, jtype, rsq, factor_coul, factor_lj, fforce);

  // Apply lambda weighting if fix is available
  if (fix_backmap) {
    double *lam = BackmapLambda::extract_lambda(fix_backmap);
    if (lam) {
      double li = BackmapLambda::clamp_lambda(lam[i]);
      double lj_val = BackmapLambda::clamp_lambda(lam[j]);
      double w = BackmapLambda::compute_weight(li, lj_val, kind == CG);
      fforce *= w;
      eng *= w;
    }
  }
  return eng;
}

/* ---------------------------------------------------------------------- */

void PairBackmap::write_restart(FILE *fp) { write_restart_settings(fp); }

void PairBackmap::read_restart(FILE *fp) {
  read_restart_settings(fp);
  allocate();
}

void PairBackmap::write_restart_settings(FILE *fp) {
  fwrite(&cut_at, sizeof(double), 1, fp);
  fwrite(&cut_cg, sizeof(double), 1, fp);
}

void PairBackmap::read_restart_settings(FILE *fp) {
  if (comm->me == 0) {
    utils::sfread(FLERR, &cut_at, sizeof(double), 1, fp, nullptr, error);
    utils::sfread(FLERR, &cut_cg, sizeof(double), 1, fp, nullptr, error);
  }
  MPI_Bcast(&cut_at, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_cg, 1, MPI_DOUBLE, 0, world);
  cut_global = MAX(cut_at, cut_cg);
}
