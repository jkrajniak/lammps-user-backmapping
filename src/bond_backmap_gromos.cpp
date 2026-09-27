/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* bond_style backmap/gromos — lambda-weighted GROMOS-96 quartic bond.

   Same functional form as bond_style gromos (GROMACS bond func 2):
     E = w × ¼ K (r² - r0²)²
     F = -w × K r (r² - r0²)   (fbond = -w·K·(r² - r0²))
   with K the GROMACS kb converted to energy/distance^4.

   w comes from the single global lambda and CG-bead co-membership:
   1 - lambda (cg), 1 (at, same bead), lambda (at, different beads).

   Syntax:
     bond_style backmap/gromos
     bond_coeff N at/cg K r0 */

#include "bond_backmap_gromos.h"

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

BondBackmapGromos::BondBackmapGromos(LAMMPS *lmp)
    : Bond(lmp),
      k(nullptr),
      r0(nullptr),
      is_cg(nullptr),
      fix_backmap(nullptr) {}

/* ---------------------------------------------------------------------- */

BondBackmapGromos::~BondBackmapGromos() {
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(r0);
    memory->destroy(is_cg);
  }
}

/* ---------------------------------------------------------------------- */

void BondBackmapGromos::allocate() {
  allocated = 1;
  int n = atom->nbondtypes + 1;

  memory->create(setflag, n, "bond:setflag");
  memory->create(k, n, "bond:k");
  memory->create(r0, n, "bond:r0");
  memory->create(is_cg, n, "bond:is_cg");

  for (int i = 1; i < n; i++) setflag[i] = 0;
}

/* ---------------------------------------------------------------------- */

void BondBackmapGromos::coeff(int narg, char **arg) {
  if (narg != 4)
    error->all(FLERR, "Incorrect args for bond_coeff backmap/gromos");
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
               "bond_coeff backmap/gromos: 2nd arg must be 'at' or 'cg'");

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

void BondBackmapGromos::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "bond_style backmap/gromos");
}

/* ---------------------------------------------------------------------- */

void BondBackmapGromos::compute(int eflag, int vflag) {
  ev_init(eflag, vflag);

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "bond_style backmap/gromos: cannot extract atom2cg/"
               "lambda_global");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

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
    double dr2 = rsq - r0[btype] * r0[btype];
    double kdr2 = k[btype] * dr2;

    // Lambda weighting
    bool same_bead = BackmapLambda::same_bead(atom2cg, i1, i2);
    double w =
        BackmapLambda::compute_weight3(same_bead, is_cg[btype], lambda_global);

    if (BackmapLambda::is_almost_zero(w)) continue;

    double fbond = -w * kdr2;

    f[i1][0] += delx * fbond;
    f[i1][1] += dely * fbond;
    f[i1][2] += delz * fbond;
    f[i2][0] -= delx * fbond;
    f[i2][1] -= dely * fbond;
    f[i2][2] -= delz * fbond;

    double ebond = 0.0;
    if (eflag) ebond = 0.25 * w * kdr2 * dr2;

    if (evflag)
      ev_tally(i1, i2, nlocal, newton_bond, ebond, fbond, delx, dely, delz);
  }
}

/* ---------------------------------------------------------------------- */

double BondBackmapGromos::equilibrium_distance(int i) { return r0[i]; }

/* ---------------------------------------------------------------------- */

double BondBackmapGromos::single(int btype, double rsq, int i, int j,
                                 double &fforce) {
  double dr2 = rsq - r0[btype] * r0[btype];
  double kdr2 = k[btype] * dr2;

  double w = 1.0;
  if (fix_backmap) {
    int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
    double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
    if (atom2cg && lam_global_ptr) {
      double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);
      bool same_bead = BackmapLambda::same_bead(atom2cg, i, j);
      w = BackmapLambda::compute_weight3(same_bead, is_cg[btype],
                                         lambda_global);
    }
  }

  fforce = -w * kdr2;
  return 0.25 * w * kdr2 * dr2;
}

/* ---------------------------------------------------------------------- */

void BondBackmapGromos::write_restart(FILE *fp) {
  fwrite(k + 1, sizeof(double), atom->nbondtypes, fp);
  fwrite(r0 + 1, sizeof(double), atom->nbondtypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->nbondtypes, fp);
}

void BondBackmapGromos::read_restart(FILE *fp) {
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

/* ---------------------------------------------------------------------- */

/* proc 0 writes to data file, one line per type in coeff() argument order */

void BondBackmapGromos::write_data(FILE *fp) {
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp, "%d %s %.15g %.15g\n", i, is_cg[i] ? "cg" : "at", k[i], r0[i]);
}
