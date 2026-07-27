/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "pair_lj_cut_coul_cut_ecap.h"

#include <cmath>

#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutEcap::PairLJCutCoulCutEcap(LAMMPS *lmp)
    : PairLJCutCoulCut(lmp),
      cap_factor(0.5),
      capradsq(nullptr),
      lj_cap_eng(nullptr) {}

/* ---------------------------------------------------------------------- */

PairLJCutCoulCutEcap::~PairLJCutCoulCutEcap() {
  if (copymode) return;
  if (allocated) {
    memory->destroy(capradsq);
    memory->destroy(lj_cap_eng);
  }
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutEcap::allocate() {
  PairLJCutCoulCut::allocate();
  int np1 = atom->ntypes + 1;
  memory->create(capradsq, np1, np1, "pair:capradsq");
  memory->create(lj_cap_eng, np1, np1, "pair:lj_cap_eng");
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutEcap::settings(int narg, char **arg) {
  if (narg < 1 || narg > 3) error->all(FLERR, "Illegal pair_style command");

  cut_lj_global = utils::numeric(FLERR, arg[0], false, lmp);
  if (narg == 1) {
    cut_coul_global = cut_lj_global;
  } else {
    cut_coul_global = utils::numeric(FLERR, arg[1], false, lmp);
  }
  cap_factor = 0.5;
  if (narg == 3) {
    cap_factor = utils::numeric(FLERR, arg[2], false, lmp);
    if (cap_factor <= 0.0)
      error->all(FLERR,
                 "pair_style lj/cut/coul/cut/ecap cap_factor must be > 0");
  }

  if (allocated) {
    for (int i = 1; i <= atom->ntypes; i++)
      for (int j = i; j <= atom->ntypes; j++)
        if (setflag[i][j]) {
          cut_lj[i][j] = cut_lj_global;
          cut_coul[i][j] = cut_coul_global;
        }
  }
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulCutEcap::init_one(int i, int j) {
  double cut = PairLJCutCoulCut::init_one(i, j);

  double caprad = cap_factor * sigma[i][j];
  double csq = caprad * caprad;
  double eng = 0.0;
  if (csq > 0.0 && epsilon[i][j] != 0.0) {
    double r2inv = 1.0 / csq;
    double r6inv = r2inv * r2inv * r2inv;
    eng = r6inv * (lj3[i][j] * r6inv - lj4[i][j]);
  } else {
    csq = 0.0;
  }

  capradsq[i][j] = capradsq[j][i] = csq;
  lj_cap_eng[i][j] = lj_cap_eng[j][i] = eng;
  return cut;
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulCutEcap::compute(int eflag, int vflag) {
  int i, j, ii, jj, inum, jnum, itype, jtype;
  double qtmp, xtmp, ytmp, ztmp, delx, dely, delz, evdwl, ecoul, fpair;
  double rsq, r2inv, r6inv, forcecoul, forcelj, factor_coul, factor_lj;
  int *ilist, *jlist, *numneigh, **firstneigh;

  evdwl = ecoul = 0.0;
  ev_init(eflag, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      factor_coul = special_coul[sbmask(j)];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      jtype = type[j];

      if (rsq < cutsq[itype][jtype]) {
        r2inv = 1.0 / rsq;

        if (rsq < cut_coulsq[itype][jtype])
          forcecoul = qqrd2e * qtmp * q[j] * sqrt(r2inv);
        else
          forcecoul = 0.0;

        if (rsq < cut_ljsq[itype][jtype]) {
          if (rsq > capradsq[itype][jtype]) {
            r6inv = r2inv * r2inv * r2inv;
            forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
          } else {
            forcelj = 0.0;
          }
        } else
          forcelj = 0.0;

        fpair = (factor_coul * forcecoul + factor_lj * forcelj) * r2inv;

        f[i][0] += delx * fpair;
        f[i][1] += dely * fpair;
        f[i][2] += delz * fpair;
        if (newton_pair || j < nlocal) {
          f[j][0] -= delx * fpair;
          f[j][1] -= dely * fpair;
          f[j][2] -= delz * fpair;
        }

        if (eflag) {
          if (rsq < cut_coulsq[itype][jtype])
            ecoul = factor_coul * qqrd2e * qtmp * q[j] * sqrt(r2inv);
          else
            ecoul = 0.0;
          if (rsq < cut_ljsq[itype][jtype]) {
            if (rsq > capradsq[itype][jtype]) {
              r6inv = r2inv * r2inv * r2inv;
              evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) -
                      offset[itype][jtype];
            } else {
              evdwl = lj_cap_eng[itype][jtype] - offset[itype][jtype];
            }
            evdwl *= factor_lj;
          } else
            evdwl = 0.0;
        }

        if (evflag)
          ev_tally(i, j, nlocal, newton_pair, evdwl, ecoul, fpair, delx, dely,
                   delz);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

double PairLJCutCoulCutEcap::single(int i, int j, int itype, int jtype,
                                    double rsq, double factor_coul,
                                    double factor_lj, double &fforce) {
  double r2inv, r6inv, forcecoul, forcelj, phicoul, philj;

  r2inv = 1.0 / rsq;
  if (rsq < cut_coulsq[itype][jtype])
    forcecoul = force->qqrd2e * atom->q[i] * atom->q[j] * sqrt(r2inv);
  else
    forcecoul = 0.0;

  if (rsq < cut_ljsq[itype][jtype]) {
    if (rsq > capradsq[itype][jtype]) {
      r6inv = r2inv * r2inv * r2inv;
      forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
      philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) -
              offset[itype][jtype];
    } else {
      forcelj = 0.0;
      philj = lj_cap_eng[itype][jtype] - offset[itype][jtype];
    }
  } else {
    forcelj = 0.0;
    philj = 0.0;
  }
  fforce = (factor_coul * forcecoul + factor_lj * forcelj) * r2inv;

  double eng = 0.0;
  if (rsq < cut_coulsq[itype][jtype]) {
    phicoul = force->qqrd2e * atom->q[i] * atom->q[j] * sqrt(r2inv);
    eng += factor_coul * phicoul;
  }
  if (rsq < cut_ljsq[itype][jtype]) eng += factor_lj * philj;
  return eng;
}
