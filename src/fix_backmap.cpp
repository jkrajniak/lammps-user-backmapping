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

   Manages: lambda ramp, CG-AT bead mapping, COM position tracking,
   CG force distribution to AT atoms.

   Syntax:
     fix ID group-ID backmap cg_type T1 [T2 ...] alpha A lambda0 L0
         [nonuniform yes/no]

   Reference: Krajniak et al., JCTC 2016, DOI: 10.1021/acs.jctc.6b00595 */

#include "fix_backmap.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "random_mars.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static const std::set<std::string> KNOWN_KEYWORDS = {"alpha", "lambda0",
                                                     "nonuniform", "phase"};

/* ---------------------------------------------------------------------- */

FixBackmap::FixBackmap(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg),
      alpha(0.0),
      lambda0(0.0),
      nonuniform(0),
      ramp_active(0),
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

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "cg_type") == 0) {
      iarg++;
      if (iarg >= narg)
        utils::missing_cmd_args(FLERR, "fix backmap cg_type", error);
      while (iarg < narg &&
             KNOWN_KEYWORDS.find(arg[iarg]) == KNOWN_KEYWORDS.end()) {
        char *endptr;
        long val = strtol(arg[iarg], &endptr, 10);
        if (*endptr != '\0' || endptr == arg[iarg])
          error->all(FLERR,
                     "fix backmap cg_type: '{}' is not a valid integer type",
                     arg[iarg]);
        int t = static_cast<int>(val);
        if (t < 1 || t > atom->ntypes)
          error->all(FLERR, "fix backmap cg_type {} out of range [1,{}]", t,
                     atom->ntypes);
        cg_types.insert(t);
        iarg++;
      }
      if (cg_types.empty())
        error->all(FLERR, "fix backmap cg_type requires at least one type");
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

  if (cg_types.empty()) error->all(FLERR, "fix backmap requires cg_type");
  if (alpha == 0.0) error->all(FLERR, "fix backmap requires alpha");

  maxatom = atom->nmax;
  memory->create(lambda, maxatom, "backmap:lambda");
  vector_atom = lambda;

  atom->add_callback(Atom::GROW);
  if (restart_peratom) atom->add_callback(Atom::RESTART);

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

FixBackmap::~FixBackmap() {
  atom->delete_callback(id, Atom::GROW);
  atom->delete_callback(id, Atom::RESTART);
  memory->destroy(lambda);
}

/* ---------------------------------------------------------------------- */

int FixBackmap::setmask() {
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= PRE_FORCE;
  mask |= POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::init() {
  // Bead map built in setup() after neighbor list establishes ghost atoms
}

/* ---------------------------------------------------------------------- */

void FixBackmap::setup(int /*vflag*/) {
  comm->forward_comm(this);

  build_bead_map();
  validate_masses();

  if (comm->me == 0) {
    std::string type_str;
    for (int t : cg_types) {
      if (!type_str.empty()) type_str += ",";
      type_str += std::to_string(t);
    }
    utils::logmesg(lmp,
                   "fix backmap: {} bead mappings, cg_types={{{}}}, alpha={}, "
                   "lambda0={}\n",
                   bead_map.size(), type_str, alpha, lambda0);
  }

  // CG position update deferred to the first timestep to avoid
  // moving atoms outside communication range during setup.
}

/* ---------------------------------------------------------------------- */

void FixBackmap::initial_integrate(int /*vflag*/) {
  if (!ramp_active) return;

  // Rebuild bead map if empty (e.g., after unfix/fix or at start of
  // a new run).  setup() builds it once, but subsequent runs need
  // to rebuild before initial_integrate uses it.
  if (bead_map.empty()) build_bead_map();

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double *mass_type = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  // Update CG bead positions to AT COM.  fix nve has already updated
  // AT positions in its initial_integrate (created before this fix),
  // so CG positions match current AT configuration when forces are
  // computed later in this timestep.
  for (auto &bm : bead_map) {
    int cg = bm.cg_local;
    if (cg < 0 || cg >= nlocal) continue;

    double com_dx = 0.0, com_dy = 0.0, com_dz = 0.0;
    double m_total = 0.0;

    for (int at : bm.at_local) {
      if (at < 0 || at >= nlocal) continue;
      double m_i = mass_type[type[at]];
      double dx = x[at][0] - x[cg][0];
      double dy = x[at][1] - x[cg][1];
      double dz = x[at][2] - x[cg][2];
      // Use round-based wrapping to handle ghost atoms that may be
      // more than one box length from the local CG bead.
      double *prd = domain->prd;
      if (domain->xperiodic) dx -= prd[0] * round(dx / prd[0]);
      if (domain->yperiodic) dy -= prd[1] * round(dy / prd[1]);
      if (domain->zperiodic) dz -= prd[2] * round(dz / prd[2]);
      com_dx += m_i * dx;
      com_dy += m_i * dy;
      com_dz += m_i * dz;
      m_total += m_i;
    }

    if (m_total > 0.0) {
      double frac = lambda[cg];
      x[cg][0] += frac * com_dx / m_total;
      x[cg][1] += frac * com_dy / m_total;
      x[cg][2] += frac * com_dz / m_total;
    }
  }

  // Zero CG velocities and forces so they don't get integrated by NVE
  nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) {
    if (is_cg_type(type[i])) {
      v[i][0] = v[i][1] = v[i][2] = 0.0;
      f[i][0] = f[i][1] = f[i][2] = 0.0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmap::pre_force(int /*vflag*/) {
  if (neighbor->ago == 0) {
    // After a neighbor rebuild, new ghost atoms may appear with
    // uninitialized lambda values. Communicate lambda to ghosts
    // before force computation reads them.
    comm->forward_comm(this);
    build_bead_map();
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmap::post_force(int /*vflag*/) {
  if (!ramp_active) return;

  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (auto &bm : bead_map) {
    int cg = bm.cg_local;
    if (cg < 0 || cg >= nlocal) continue;

    double fx_cg = f[cg][0];
    double fy_cg = f[cg][1];
    double fz_cg = f[cg][2];
    double m_cg = mass[type[cg]];

    if (m_cg <= 0.0) continue;

    for (int at : bm.at_local) {
      if (at < 0 || at >= nlocal) continue;
      double ratio = mass[type[at]] / m_cg;
      f[at][0] += ratio * fx_cg;
      f[at][1] += ratio * fy_cg;
      f[at][2] += ratio * fz_cg;
    }

    f[cg][0] = f[cg][1] = f[cg][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmap::end_of_step() {
  int nlocal = atom->nlocal;

  if (ramp_active) {
    for (int i = 0; i < nlocal; i++) {
      lambda[i] += alpha;
      if (lambda[i] > 1.0) lambda[i] = 1.0;
    }
  }

  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

int FixBackmap::modify_param(int narg, char **arg) {
  if (narg < 2) return 0;

  if (strcmp(arg[0], "active") == 0) {
    ramp_active = utils::logical(FLERR, arg[1], false, lmp);
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void *FixBackmap::extract(const char *str, int &dim) {
  if (strcmp(str, "lambda") == 0) {
    dim = 1;
    return static_cast<void *>(lambda);
  }
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::grow_arrays(int nmax) {
  int old_maxatom = maxatom;
  maxatom = nmax;
  memory->grow(lambda, maxatom, "backmap:lambda");
  for (int i = old_maxatom; i < maxatom; i++) lambda[i] = 0.0;
  vector_atom = lambda;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::copy_arrays(int i, int j, int /*delflag*/) {
  lambda[j] = lambda[i];
}

/* ---------------------------------------------------------------------- */

int FixBackmap::pack_exchange(int i, double *buf) {
  buf[0] = lambda[i];
  return 1;
}

int FixBackmap::unpack_exchange(int nlocal, double *buf) {
  lambda[nlocal] = buf[0];
  return 1;
}

/* ---------------------------------------------------------------------- */

int FixBackmap::pack_forward_comm(int n, int *list, double *buf,
                                  int /*pbc_flag*/, int * /*pbc*/) {
  for (int i = 0; i < n; i++) buf[i] = lambda[list[i]];
  return n;
}

void FixBackmap::unpack_forward_comm(int n, int first, double *buf) {
  for (int i = 0; i < n; i++) lambda[first + i] = buf[i];
}

/* ---------------------------------------------------------------------- */

int FixBackmap::pack_restart(int i, double *buf) {
  buf[0] = 2;
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

/* Build the per-bead map: for each CG bead, find its mapped AT atoms.

   Convention: within each molecule (by molecule ID), CG atoms are listed
   first (sorted by tag), AT atoms next (sorted by tag). The AT atoms are
   distributed evenly among CG beads in tag order:
     CG bead i -> AT atoms [i*apb, (i+1)*apb)
   where apb = num_AT / num_CG.

   Only beads whose CG atom is local are stored (AT atoms may be ghost). */

void FixBackmap::build_bead_map() {
  bead_map.clear();

  int *type = atom->type;
  tagint *tag = atom->tag;
  tagint *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;

  // Group atoms by molecule ID, separating CG and AT.
  // Deduplicate by tag: a local+ghost pair for the same physical atom
  // should only appear once. Prefer local index over ghost.
  struct AtomRef {
    tagint gtag;
    int local_idx;
    bool is_local;
  };
  struct MolInfo {
    std::map<tagint, AtomRef> cg_by_tag;
    std::map<tagint, AtomRef> at_by_tag;
    bool has_local_cg = false;
  };

  std::map<tagint, MolInfo> mol_atoms;

  for (int i = 0; i < ntotal; i++) {
    tagint mol_id = molecule[i];
    if (mol_id <= 0) continue;
    bool local = (i < nlocal);

    if (is_cg_type(type[i])) {
      auto &entry = mol_atoms[mol_id].cg_by_tag[tag[i]];
      if (entry.gtag == 0 || (local && !entry.is_local)) {
        entry = {tag[i], i, local};
      }
      if (local) mol_atoms[mol_id].has_local_cg = true;
    } else {
      auto &entry = mol_atoms[mol_id].at_by_tag[tag[i]];
      if (entry.gtag == 0 || (local && !entry.is_local)) {
        entry = {tag[i], i, local};
      }
    }
  }

  for (auto &[mol_id, info] : mol_atoms) {
    if (!info.has_local_cg) continue;
    if (info.cg_by_tag.empty()) continue;

    // Flatten deduplicated maps into sorted vectors
    std::vector<AtomRef> cg_atoms, at_atoms;
    for (auto &[t, ref] : info.cg_by_tag) cg_atoms.push_back(ref);
    for (auto &[t, ref] : info.at_by_tag) at_atoms.push_back(ref);

    int n_cg = static_cast<int>(cg_atoms.size());
    int n_at = static_cast<int>(at_atoms.size());

    if (n_at == 0) continue;

    if (n_at % n_cg != 0)
      error->all(FLERR,
                 "fix backmap: molecule {} has {} AT atoms and {} CG beads, "
                 "AT count must be divisible by CG count",
                 mol_id, n_at, n_cg);

    int apb = n_at / n_cg;

    // Sort by global tag to get consistent ordering
    std::sort(
        cg_atoms.begin(), cg_atoms.end(),
        [](const AtomRef &a, const AtomRef &b) { return a.gtag < b.gtag; });
    std::sort(
        at_atoms.begin(), at_atoms.end(),
        [](const AtomRef &a, const AtomRef &b) { return a.gtag < b.gtag; });

    for (int ci = 0; ci < n_cg; ci++) {
      int cg_idx = cg_atoms[ci].local_idx;
      if (cg_idx >= nlocal) continue;

      BeadMap bm;
      bm.cg_local = cg_idx;
      bm.at_mass_sum = 0.0;

      for (int ai = ci * apb; ai < (ci + 1) * apb; ai++) {
        int at_idx = at_atoms[ai].local_idx;
        bm.at_local.push_back(at_idx);
        bm.at_mass_sum += atom->mass[type[at_idx]];
      }

      bead_map.push_back(std::move(bm));
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmap::validate_masses() {
  double *mass_type = atom->mass;
  int *type = atom->type;
  constexpr double tol = 1.0e-4;
  int nwarn = 0;

  for (auto &bm : bead_map) {
    int cg = bm.cg_local;
    if (cg < 0) continue;
    double cg_mass = mass_type[type[cg]];
    double at_mass_sum = bm.at_mass_sum;
    if (std::fabs(cg_mass - at_mass_sum) > tol) {
      nwarn++;
      if (nwarn <= 5)
        error->warning(FLERR,
                       "fix backmap: CG bead (tag {}, type {}, mass {:.4f}) != "
                       "AT mass sum ({:.4f})",
                       atom->tag[cg], type[cg], cg_mass, at_mass_sum);
    }
  }
  if (nwarn > 5)
    error->warning(FLERR,
                   "fix backmap: {} additional mass mismatches not shown",
                   nwarn - 5);
}
