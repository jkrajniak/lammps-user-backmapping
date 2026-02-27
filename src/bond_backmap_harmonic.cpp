/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* bond_style backmap/harmonic — lambda-weighted harmonic bond for backmapping.

   F = -w × k × (r - r0)
   E = w × 0.5 × k × (r - r0)²

   where w = λ_i × λ_j (at) or 1 − λ_i × λ_j (cg).

   Syntax:
     bond_style backmap/harmonic
     bond_coeff N at/cg K r0 */

#include "bond_backmap_harmonic.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondBackmapHarmonic::BondBackmapHarmonic(LAMMPS *lmp)
    : Bond(lmp),
      k(nullptr),
      r0(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

BondBackmapHarmonic::~BondBackmapHarmonic() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void BondBackmapHarmonic::allocate() {
  allocated = 1;
  int n = atom->nbondtypes + 1;

  memory->create(setflag, n, "bond:setflag");
  memory->create(k, n, "bond:k");
  memory->create(r0, n, "bond:r0");
  memory->create(is_cg, n, "bond:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void BondBackmapHarmonic::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for bond_coeff backmap/harmonic");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  int cg_flag;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "bond_coeff backmap/harmonic: 2nd arg must be 'at' or 'cg'");

  double k_one = utils::numeric(FLERR, arg[2], false, lmp);
  double r0_one = utils::numeric(FLERR, arg[3], false, lmp);

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void BondBackmapHarmonic::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "bond_style backmap/harmonic");
}

/* ---------------------------------------------------------------------- */

void BondBackmapHarmonic::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  double *lam = BackmapLambda::extract_lambda(fix_backmap);
  if (!lam)
    error->all(FLERR, "bond_style backmap/harmonic: cannot extract lambda");

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
    double dr = r - r0[btype];
    double rk = k[btype] * dr;

    // Lambda weighting
    double li = BackmapLambda::clamp_lambda(lam[i1]);
    double lj = BackmapLambda::clamp_lambda(lam[i2]);
    double w = BackmapLambda::compute_weight(li, lj, is_cg[btype]);

    if (BackmapLambda::is_almost_zero(w)) continue;

    double fbond;
    if (r > 0.0)
      fbond = -2.0 * w * rk / r;
    else
      fbond = 0.0;

    f[i1][0] += delx * fbond;
    f[i1][1] += dely * fbond;
    f[i1][2] += delz * fbond;
    f[i2][0] -= delx * fbond;
    f[i2][1] -= dely * fbond;
    f[i2][2] -= delz * fbond;

    double ebond = 0.0;
    if (eflag) ebond = w * rk * dr;

    if (evflag)
      ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

double BondBackmapHarmonic::equilibrium_distance(int i) { return r0[i]; }

/* ---------------------------------------------------------------------- */

double BondBackmapHarmonic::single(int btype, double rsq, int i, int j,
                                   double &fforce) {
  double r = sqrt(rsq);
  double dr = r - r0[btype];
  double rk = k[btype] * dr;

  double w = 1.0;
  if (fix_backmap) {
    double *lam = BackmapLambda::extract_lambda(fix_backmap);
    if (lam) {
      double li = BackmapLambda::clamp_lambda(lam[i]);
      double lj = BackmapLambda::clamp_lambda(lam[j]);
      w = BackmapLambda::compute_weight(li, lj, is_cg[btype]);
    }
  }

  fforce = 0.0;
  if (r > 0.0) fforce = -2.0 * w * rk / r;
  return w * rk * dr;
}

/* ---------------------------------------------------------------------- */

void BondBackmapHarmonic::write_restart(FILE *fp) {
  fwrite(k + 1, sizeof(double), atom->nbondtypes, fp);
  fwrite(r0 + 1, sizeof(double), atom->nbondtypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nbondtypes, fp);
}

void BondBackmapHarmonic::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, k + 1, sizeof(double), atom->nbondtypes, fp, nullptr,
                  error);
    utils::sfread(FLERR, r0 + 1, sizeof(double), atom->nbondtypes, fp, nullptr,
                  error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->nbondtypes, fp, nullptr,
                  error);
  }
  MPI_Bcast(k + 1, atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(r0 + 1, atom->nbondtypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(is_cg + 1, atom->nbondtypes, MPI_INT, 0, world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}
