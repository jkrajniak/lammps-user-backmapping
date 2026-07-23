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
      lambda_global(0.0),
      maxatom(0),
      atom2cg(nullptr),
      com_buf(nullptr),
      cg_denom(nullptr),
      cg_fwd(nullptr) {
  if (narg < 7) utils::missing_cmd_args(FLERR, "fix backmap", error);

  restart_peratom = 1;
  restart_global = 0;
  peratom_flag = 1;
  size_peratom_cols = 0;
  peratom_freq = 1;
  scalar_flag = 1;
  extscalar = 0;
  comm_forward = 5;  // lambda + CG (fx, fy, fz, at_mass_sum)
  comm_reverse = 4;  // COM (mass, mdx, mdy, mdz)
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
  memory->create(atom2cg, maxatom, "backmap:atom2cg");
  memory->create(com_buf, maxatom * 4, "backmap:com_buf");
  memory->create(cg_denom, maxatom, "backmap:cg_denom");
  memory->create(cg_fwd, maxatom * 4, "backmap:cg_fwd");
  for (int i = 0; i < maxatom; i++) {
    atom2cg[i] = -1;
    cg_denom[i] = 0.0;
  }
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

  // lambda_global follows the uniform ramp rule regardless of nonuniform
  // per-atom staggering; it is the single source of truth for weighting.
  lambda_global = lambda0;
}

/* ---------------------------------------------------------------------- */

FixBackmap::~FixBackmap() {
  atom->delete_callback(id, Atom::GROW);
  atom->delete_callback(id, Atom::RESTART);
  memory->destroy(lambda);
  memory->destroy(atom2cg);
  memory->destroy(com_buf);
  memory->destroy(cg_denom);
  memory->destroy(cg_fwd);
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
  // CG–AT kinematic coupling always runs.  fix_modify active yes/no only
  // controls the lambda ramp in end_of_step (matches E++ DynamicResolution
  // vs VelocityVerletHybrid and AdResS virtual-site semantics).

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
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;
  double *prd = domain->prd;

  // --- COM update (MPI-correct): each local AT atom contributes to its CG
  // bead's accumulator; reverse_comm sums ghost-CG contributions to CG owners.
  zero_com_buf();

  for (int i = 0; i < nlocal; i++) {
    if (is_cg_type(type[i])) continue;
    int cg = atom2cg[i];
    if (cg < 0 || cg >= ntotal) continue;
    double m_i = mass_type[type[i]];
    double dx = x[i][0] - x[cg][0];
    double dy = x[i][1] - x[cg][1];
    double dz = x[i][2] - x[cg][2];
    if (domain->xperiodic) dx -= prd[0] * round(dx / prd[0]);
    if (domain->yperiodic) dy -= prd[1] * round(dy / prd[1]);
    if (domain->zperiodic) dz -= prd[2] * round(dz / prd[2]);
    int k = cg * 4;
    com_buf[k + 0] += m_i;
    com_buf[k + 1] += m_i * dx;
    com_buf[k + 2] += m_i * dy;
    com_buf[k + 3] += m_i * dz;
  }

  // Sum ghost-CG contributions back to local CG owners.
  comm->reverse_comm(this, 4);

  // Apply COM shift on local CG beads.
  for (int cg = 0; cg < nlocal; cg++) {
    if (!is_cg_type(type[cg])) continue;
    int k = cg * 4;
    double m_total = com_buf[k + 0];
    if (m_total <= 0.0) continue;
    x[cg][0] += com_buf[k + 1] / m_total;
    x[cg][1] += com_buf[k + 2] / m_total;
    x[cg][2] += com_buf[k + 3] / m_total;
  }

  // Zero CG velocities and forces so they don't get integrated by NVE
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

void FixBackmap::setup_pre_force(int vflag) {
  // Verlet::setup() calls modify->setup_pre_force() BEFORE the initial
  // bond/angle/pair force evaluation, but modify->setup() -- where
  // build_bead_map() historically only ran -- executes AFTER that
  // evaluation. Without this override, atom2cg is still freshly
  // constructed (all -1) for the very first force evaluation of every
  // run, so every same-bead AT bond/angle is misclassified as inter-bead
  // by BackmapLambda::same_bead() and gets zero-weighted at
  // lambda_global=0 (compute_weight3 falls through to the
  // different-bead branch, w=lambda_global). That corrupts the very
  // first velocity half-kick every AT atom receives. Delegate to
  // pre_force() so the bead map is ready before any force computation,
  // matching every subsequent step.
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixBackmap::post_force(int /*vflag*/) {
  // CG force distribution always runs; lambda weighting is applied by the
  // backmap interaction styles before forces reach this hook.

  // Forward-comm lambda + CG (fx, fy, fz, at_mass_sum) so each local AT atom
  // can read its (local or ghost) CG bead's complete force.
  comm->forward_comm(this);

  double **f = atom->f;
  double *mass = atom->mass;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nghost = atom->nghost;
  int ntotal = nlocal + nghost;

  // Distribute CG force to each LOCAL AT atom. All additions target local
  // atoms, so no reverse_comm is needed.
  for (int i = 0; i < nlocal; i++) {
    if (is_cg_type(type[i])) continue;
    int cg = atom2cg[i];
    if (cg < 0 || cg >= ntotal) continue;

    double fx_cg, fy_cg, fz_cg, m_bead;
    if (cg < nlocal) {
      fx_cg = f[cg][0];
      fy_cg = f[cg][1];
      fz_cg = f[cg][2];
      m_bead = cg_denom[cg];
    } else {
      const double *p = &cg_fwd[cg * 4];
      fx_cg = p[0];
      fy_cg = p[1];
      fz_cg = p[2];
      m_bead = p[3];
    }
    if (m_bead <= 0.0) continue;
    double ratio = mass[type[i]] / m_bead;
    f[i][0] += ratio * fx_cg;
    f[i][1] += ratio * fy_cg;
    f[i][2] += ratio * fz_cg;
  }

  // Zero CG forces on local CG beads after distribution.
  for (int i = 0; i < nlocal; i++) {
    if (is_cg_type(type[i])) f[i][0] = f[i][1] = f[i][2] = 0.0;
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
    lambda_global += alpha;
    if (lambda_global > 1.0) lambda_global = 1.0;
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

double FixBackmap::compute_scalar() {
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
  double sum = 0.0;
  int n = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      sum += lambda[i];
      n++;
    }
  }
  double all_sum = 0.0;
  int all_n = 0;
  MPI_Allreduce(&sum, &all_sum, 1, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(&n, &all_n, 1, MPI_INT, MPI_SUM, world);
  if (all_n == 0) return 0.0;
  return all_sum / static_cast<double>(all_n);
}

/* ---------------------------------------------------------------------- */

void *FixBackmap::extract(const char *str, int &dim) {
  if (strcmp(str, "lambda") == 0) {
    dim = 1;
    return static_cast<void *>(lambda);
  }
  if (strcmp(str, "atom2cg") == 0) {
    dim = 1;
    return static_cast<void *>(atom2cg);
  }
  if (strcmp(str, "lambda_global") == 0) {
    dim = 0;
    return static_cast<void *>(&lambda_global);
  }
  return nullptr;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::grow_arrays(int nmax) {
  int old_maxatom = maxatom;
  maxatom = nmax;
  memory->grow(lambda, maxatom, "backmap:lambda");
  memory->grow(atom2cg, maxatom, "backmap:atom2cg");
  memory->grow(com_buf, maxatom * 4, "backmap:com_buf");
  memory->grow(cg_denom, maxatom, "backmap:cg_denom");
  memory->grow(cg_fwd, maxatom * 4, "backmap:cg_fwd");
  for (int i = old_maxatom; i < maxatom; i++) {
    lambda[i] = 0.0;
    atom2cg[i] = -1;
    cg_denom[i] = 0.0;
  }
  vector_atom = lambda;
}

/* ---------------------------------------------------------------------- */

void FixBackmap::copy_arrays(int i, int j, int /*delflag*/) {
  lambda[j] = lambda[i];
  atom2cg[j] = atom2cg[i];
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
  double **f = atom->f;
  int *type = atom->type;
  int m = 0;
  for (int k = 0; k < n; k++) {
    int i = list[k];
    buf[m++] = lambda[i];
    if (is_cg_type(type[i])) {
      buf[m++] = f[i][0];
      buf[m++] = f[i][1];
      buf[m++] = f[i][2];
      buf[m++] = cg_denom[i];
    } else {
      buf[m++] = 0.0;
      buf[m++] = 0.0;
      buf[m++] = 0.0;
      buf[m++] = 0.0;
    }
  }
  return m;
}

void FixBackmap::unpack_forward_comm(int n, int first, double *buf) {
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    lambda[i] = buf[m++];
    int k = i * 4;
    cg_fwd[k + 0] = buf[m++];
    cg_fwd[k + 1] = buf[m++];
    cg_fwd[k + 2] = buf[m++];
    cg_fwd[k + 3] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int FixBackmap::pack_reverse_comm(int n, int first, double *buf) {
  int m = 0;
  int last = first + n;
  for (int i = first; i < last; i++) {
    const double *p = &com_buf[i * 4];
    buf[m++] = p[0];
    buf[m++] = p[1];
    buf[m++] = p[2];
    buf[m++] = p[3];
  }
  return m;
}

void FixBackmap::unpack_reverse_comm(int n, int *list, double *buf) {
  int m = 0;
  for (int k = 0; k < n; k++) {
    int i = list[k];
    double *p = &com_buf[i * 4];
    p[0] += buf[m++];
    p[1] += buf[m++];
    p[2] += buf[m++];
    p[3] += buf[m++];
  }
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
  return static_cast<double>(maxatom) *
         (2 * sizeof(double) + sizeof(int) + 4 * sizeof(double) +
          4 * sizeof(double));
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

  // Reset the per-atom CG partner map for all local+ghost atoms.
  for (int i = 0; i < ntotal; i++) {
    atom2cg[i] = -1;
    cg_denom[i] = 0.0;
  }

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
    if (info.cg_by_tag.empty()) continue;
    if (info.at_by_tag.empty()) continue;

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

      // Populate atom2cg for every AT atom of this bead (local or ghost).
      for (int ai = ci * apb; ai < (ci + 1) * apb; ai++) {
        int at_idx = at_atoms[ai].local_idx;
        atom2cg[at_idx] = cg_idx;
      }

      // bead_map entry kept for setup logging + validate_masses (local CG
      // only).
      if (cg_idx >= nlocal) continue;

      BeadMap bm;
      bm.cg_local = cg_idx;
      bm.at_mass_sum = 0.0;
      for (int ai = ci * apb; ai < (ci + 1) * apb; ai++) {
        int at_idx = at_atoms[ai].local_idx;
        bm.at_local.push_back(at_idx);
        bm.at_mass_sum += atom->mass[type[at_idx]];
      }
      cg_denom[cg_idx] = bm.at_mass_sum;
      bead_map.push_back(std::move(bm));
    }
  }

  // Warn if any local AT atom has no resolvable CG partner on this rank,
  // which indicates comm_modify cutoff is too small for the bead extent.
  int nmissing = 0;
  for (int i = 0; i < nlocal; i++) {
    if (!is_cg_type(type[i]) && atom2cg[i] < 0) nmissing++;
  }
  if (nmissing > 0)
    error->warning(FLERR,
                   "fix backmap: {} local AT atom(s) have no CG partner in "
                   "local+ghost range; increase comm_modify cutoff",
                   nmissing);
}

/* ---------------------------------------------------------------------- */

void FixBackmap::zero_com_buf() {
  int nghost = atom->nghost;
  int ntotal = atom->nlocal + nghost;
  for (int i = 0; i < ntotal * 4; i++) com_buf[i] = 0.0;
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
