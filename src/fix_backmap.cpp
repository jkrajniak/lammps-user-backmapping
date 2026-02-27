/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* fix backmap — time-dependent backmapping from CG to AT resolution.

   Manages: lambda ramp, CG-AT molecule mapping, COM position tracking,
   CG force distribution to AT atoms.

   Syntax:
     fix ID group-ID backmap cg_type T alpha A lambda0 L0 [nonuniform yes/no]

   Reference: Krajniak et al., JCTC 2016, DOI: 10.1021/acs.jctc.6b00595 */

#include "fix_backmap.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "random_mars.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBackmap::FixBackmap(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg),
      cg_type(0),
      alpha(0.0),
      lambda0(0.0),
      nonuniform(0),
      ramp_active(1),
      lambda(nullptr),
      maxatom(0) {
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix backmap", error);

  restart_peratom = 1;
  restart_global = 0;
  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  comm_forward = 1;
  create_attribute = 1;

  // Parse required arguments: cg_type T alpha A
  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cg_type") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix backmap cg_type", error);
      cg_type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (cg_type < 1 || cg_type > atom->ntypes)
        error->all(FLERR, "fix backmap cg_type {} out of range [1,{}]", cg_type,
                   atom->ntypes);
      iarg += 2;
    } else if (strcmp(arg[iarg], "alpha") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix backmap alpha", error);
      alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (alpha <= 0.0) error->all(FLERR, "fix backmap alpha must be positive");
      iarg += 2;
    } else if (strcmp(arg[iarg], "lambda0") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix backmap lambda0", error);
      lambda0 = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "nonuniform") == 0) {
      if (iarg + 1 >= narg)
        utils::missing_cmd_args(FLERR, "fix backmap nonuniform", error);
      nonuniform = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, "Illegal fix backmap argument: {}", arg[iarg]);
    }
  }

  if (cg_type == 0) error->all(FLERR, "fix backmap requires cg_type");
  if (alpha == 0.0) error->all(FLERR, "fix backmap requires alpha");

  // Allocate per-atom lambda array
  maxatom = atom->nmax;
  memory->create(lambda, maxatom, "backmap:lambda");

  // Initialize lambda values
  if (nonuniform) {
    auto *rng = new RanMars(lmp, 12345 + comm->me);
    for (int i = 0; i < atom->nlocal; i++) {
      lambda[i] = lambda0 + rng->uniform() * (-10000.0 * alpha);
    }
    delete rng;
  } else {
    for (int i = 0; i < atom->nlocal; i++) lambda[i] = lambda0;
  }
  for (int i = atom->nlocal; i < maxatom; i++) lambda[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

FixBackmap::~FixBackmap() { memory->destroy(lambda); }

/* ---------------------------------------------------------------------- */

int FixBackmap::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::init() {
  build_molecule_map();
  validate_masses();

  if (comm->me == 0)
    utils::logmesg(
        lmp,
        "fix backmap: mapped {} molecules, cg_type={}, alpha={}, lambda0={}\n",
        mol_map.size(), cg_type, alpha, lambda0);
}

/* ---------------------------------------------------------------------- */

void FixBackmap::setup(int /*vflag*/) {
  // Communicate lambda to ghost atoms
  comm->forward_comm(this);

  // Ensure CG positions are at COM at the start
  end_of_step();
}

/* ---------------------------------------------------------------------- */

/* B2: initial_integrate — zero CG velocities and forces as safety net */

void FixBackmap::initial_integrate(int /*vflag*/) {
  int *type = atom->type;
  double **v = atom->v;
  double **f = atom->f;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (type[i] == cg_type) {
      v[i][0] = v[i][1] = v[i][2] = 0.0;
      f[i][0] = f[i][1] = f[i][2] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

/* B3: post_force — distribute CG forces to AT atoms proportional to mass ratio,
   then zero CG forces.
   CG forces are already lambda-weighted by the interaction styles. */

void FixBackmap::post_force(int /*vflag*/) {
  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;

  for (auto &[mol_id, mm] : mol_map) {
    int cg = mm.cg_local;
    if (cg < 0 || cg >= atom->nlocal) continue;

    double fx_cg = f[cg][0];
    double fy_cg = f[cg][1];
    double fz_cg = f[cg][2];
    double m_cg = mm.cg_mass;

    if (m_cg <= 0.0) continue;

    for (int at : mm.at_local) {
      if (at < 0) continue;
      double ratio = mass[type[at]] / m_cg;
      f[at][0] += ratio * fx_cg;
      f[at][1] += ratio * fy_cg;
      f[at][2] += ratio * fz_cg;
    }

    f[cg][0] = f[cg][1] = f[cg][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

/* B4: end_of_step — update CG positions to COM of AT atoms,
   then increment lambda. */

void FixBackmap::end_of_step() {
  double **x = atom->x;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // Update CG positions to COM of AT atoms
  for (auto &[mol_id, mm] : mol_map) {
    int cg = mm.cg_local;
    if (cg < 0 || cg >= nlocal) continue;

    double com_dx = 0.0, com_dy = 0.0, com_dz = 0.0;
    double m_total = 0.0;

    for (int at : mm.at_local) {
      if (at < 0) continue;
      double m_i = mass[type[at]];

      double dx = x[at][0] - x[cg][0];
      double dy = x[at][1] - x[cg][1];
      double dz = x[at][2] - x[cg][2];
      domain->minimum_image(dx, dy, dz);

      com_dx += m_i * dx;
      com_dy += m_i * dy;
      com_dz += m_i * dz;
      m_total += m_i;
    }

    if (m_total > 0.0) {
      x[cg][0] += com_dx / m_total;
      x[cg][1] += com_dy / m_total;
      x[cg][2] += com_dz / m_total;
    }
  }

  // Increment lambda for all local atoms (if ramp is active)
  if (ramp_active) {
    for (int i = 0; i < nlocal; i++) {
      lambda[i] += alpha;
      if (lambda[i] > 1.0) lambda[i] = 1.0;
    }
  }

  // Communicate lambda to ghost atoms
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

/* B5: fix_modify support — activate/deactivate lambda ramp */

int FixBackmap::modify_param(int narg, char **arg) {
  if (narg < 2) return 0;

  if (strcmp(arg[0], "active") == 0) {
    ramp_active = utils::logical(FLERR, arg[1], false, lmp);
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

/* B5: extract — return per-atom lambda array for other styles */

void *FixBackmap::extract(const char *str, int &dim) {
  if (strcmp(str, "lambda") == 0) {
    dim = 1;
    return static_cast<void *>(lambda);
  }
  return nullptr;
}

/* ---------------------------------------------------------------------- */

/* B6: grow_arrays — handle atom array reallocation */

void FixBackmap::grow_arrays(int nmax) {
  maxatom = nmax;
  memory->grow(lambda, maxatom, "backmap:lambda");
  vector_atom = lambda;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::copy_arrays(int i, int j, int /*delflag*/) {
  lambda[j] = lambda[i];
}

/* ---------------------------------------------------------------------- */

/* B6: pack/unpack exchange — atom migration between processors */

int FixBackmap::pack_exchange(int i, double *buf) {
  buf[0] = lambda[i];
  return 1;
}

int FixBackmap::unpack_exchange(int nlocal, double *buf) {
  lambda[nlocal] = buf[0];
  return 1;
}

/* ---------------------------------------------------------------------- */

/* B6: forward comm — communicate lambda to ghost atoms */

int FixBackmap::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/) {
  for (int i = 0; i < n; i++) buf[i] = lambda[list[i]];
  return n;
}

void FixBackmap::unpack_forward_comm(int n, int first, double *buf) {
  for (int i = 0; i < n; i++) lambda[first + i] = buf[i];
}

/* ---------------------------------------------------------------------- */

/* B7: restart support */

int FixBackmap::pack_restart(int i, double *buf) {
  buf[0] = 2;  // size of packed data (including this count)
  buf[1] = lambda[i];
  return 2;
}

void FixBackmap::unpack_restart(int nlocal, int nth) {
  double **extra = atom->extra;
  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int>(extra[nlocal][m]);
  lambda[nlocal] = extra[nlocal][m + 1];
}

int FixBackmap::size_restart(int /*nlocal*/) { return 2; }

int FixBackmap::maxsize_restart() { return 2; }

/* ---------------------------------------------------------------------- */

double FixBackmap::memory_usage() {
  return static_cast<double>(maxatom) * sizeof(double);
}

/* ---------------------------------------------------------------------- */

/* Build the molecule map: for each molecule, find the CG atom and AT atoms.
   Called during init() and should be called after domain decomposition
   (atoms migrating between processors). */

void FixBackmap::build_molecule_map() {
  mol_map.clear();

  int *type = atom->type;
  tagint *molecule = atom->molecule;
  double *mass_type = atom->mass;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;

  // First pass: identify CG atoms (local only own the molecule)
  for (int i = 0; i < nlocal; i++) {
    if (type[i] == cg_type) {
      tagint mol_id = molecule[i];
      mol_map[mol_id].cg_local = i;
      mol_map[mol_id].cg_mass = 0.0;
    }
  }

  // Second pass: assign AT atoms to their molecule's CG atom
  for (int i = 0; i < ntotal; i++) {
    if (type[i] == cg_type) continue;
    tagint mol_id = molecule[i];
    auto it = mol_map.find(mol_id);
    if (it == mol_map.end()) continue;
    it->second.at_local.push_back(i);
    if (i < nlocal) it->second.cg_mass += mass_type[type[i]];
  }

  // For ghost AT atoms, we need the mass from the type array
  // (mass_type is per-type, always available)
  for (auto &[mol_id, mm] : mol_map) {
    double total = 0.0;
    for (int at : mm.at_local) total += mass_type[type[at]];
    mm.cg_mass = total;
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmap::validate_masses() {
  double *mass_type = atom->mass;
  int *type = atom->type;
  constexpr double tol = 1.0e-4;
  int nwarn = 0;

  for (auto &[mol_id, mm] : mol_map) {
    int cg = mm.cg_local;
    if (cg < 0) continue;
    double cg_mass = mass_type[type[cg]];
    double at_mass_sum = mm.cg_mass;
    if (std::fabs(cg_mass - at_mass_sum) > tol) {
      nwarn++;
      if (nwarn <= 5)
        error->warning(FLERR,
                       "fix backmap: molecule {} CG mass ({:.4f}) != "
                       "AT mass sum ({:.4f})",
                       mol_id, cg_mass, at_mass_sum);
    }
  }
  if (nwarn > 5)
    error->warning(FLERR,
                   "fix backmap: {} additional mass mismatches not shown",
                   nwarn - 5);
}
