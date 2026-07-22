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

   Delegates force computation to AT and CG sub-styles, then weights by a
   single global lambda scalar and CG-bead co-membership:
     CG:                                   w = 1 - lambda_global
     AT intra-bead (same CG bead):        w = 1 always
     AT inter-bead (different CG beads):  w = lambda_global

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

// Inter-bead (intermolecular) AT-AT pair interactions are lambda-weighted
// from lambda_global=0 via compute_weight3(), which already gives a smooth
// onset (w=0 at lambda=0, ramping linearly to w=1). A nonzero hard
// threshold here would keep this term fully OFF (not just weighted near
// zero) until lambda_global crosses it, then suddenly evaluate the real,
// unscaled sub-potential at whatever separation the un-relaxed fragment
// placement happens to have -- a genuine force discontinuity, not a
// smoothing device. This caused PET/Dacron to blow up the instant
// lambda_global crossed a nonzero threshold (T: 300 -> 71,000+ K within a
// few hundred steps), independent of ramp dt/thermostat tuning -- see
// notebook 2026-07-13_at-intrabead-always-on-fix.md. This exact regression
// was already diagnosed and fixed once before, on the unmerged
// feat/pete-example branch (commit eec9dc9, 2026-07-11): "Smooth
// compute_weight() onset suffices; ESPResSo++ has no threshold." Keep it
// at 0.0 -- do not reintroduce a nonzero value without re-deriving why.
static constexpr double LAMBDA_AT_ONSET = 0.0;

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
    // Forward to AT sub-style: rebuild args as [I, J, params...]
    if (pair_at) {
      int sub_narg = narg - 1;
      auto sub_arg = new char *[sub_narg];
      sub_arg[0] = arg[0];
      sub_arg[1] = arg[1];
      for (int k = 3; k < narg; k++) sub_arg[k - 1] = arg[k];
      pair_at->coeff(sub_narg, sub_arg);
      delete[] sub_arg;
    }
    for (int i = ilo; i <= ihi; i++) {
      for (int j = MAX(jlo, i); j <= jhi; j++) {
        pair_kind[i][j] = ATOMISTIC;
        setflag[i][j] = 1;
        cut[i][j] = cut_at;
      }
    }
  } else if (strcmp(arg[2], "cg") == 0) {
    // Forward to CG sub-style: rebuild args as [I, J, params...]
    if (pair_cg) {
      int sub_narg = narg - 1;
      auto sub_arg = new char *[sub_narg];
      sub_arg[0] = arg[0];
      sub_arg[1] = arg[1];
      for (int k = 3; k < narg; k++) sub_arg[k - 1] = arg[k];
      pair_cg->coeff(sub_narg, sub_arg);
      delete[] sub_arg;
    }
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

  // Sub-styles are only used via single() — they don't need their own
  // neighbor lists.  Calling their init_style() creates requests that
  // LAMMPS may merge/copy into ours, causing type-filtered lists.
  // We request a single full list for pair_backmap.
  neighbor->add_request(this, NeighConst::REQ_DEFAULT);
}

/* ---------------------------------------------------------------------- */

double PairBackmap::init_one(int i, int j) {
  if (setflag[i][j] == 0) error->all(FLERR, "All pair coeffs are not set");

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
   - CG type pairs: compute CG force, weight by w_CG = 1 - lambda_global
   - AT type pairs, same CG bead: always full strength, unconditionally
     (real intra-molecular chemistry, independent of lambda_global)
   - AT type pairs, different CG beads: weight by w_AT = lambda_global,
     smoothly from lambda_global=0 (LAMBDA_AT_ONSET=0.0 -- do not
     reintroduce a nonzero hard threshold, see comment above)
   - Cross-type or NONE pairs: skip */

void PairBackmap::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  if (!fix_backmap)
    fix_backmap = BackmapLambda::find_fix_backmap(lmp, "pair_style backmap");

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "pair_style backmap: cannot extract atom2cg/lambda_global "
               "from fix backmap");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

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

      bool is_cg = (kind == CG);
      bool same_bead = BackmapLambda::same_bead(atom2cg, i, j);
      double w =
          BackmapLambda::compute_weight3(same_bead, is_cg, lambda_global);

      if (BackmapLambda::is_almost_zero(w)) continue;

      if (!is_cg && !same_bead && lambda_global < LAMBDA_AT_ONSET) continue;

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
    int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
    double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
    if (atom2cg && lam_global_ptr) {
      double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);
      bool same_bead = BackmapLambda::same_bead(atom2cg, i, j);
      double w =
          BackmapLambda::compute_weight3(same_bead, kind == CG, lambda_global);
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
