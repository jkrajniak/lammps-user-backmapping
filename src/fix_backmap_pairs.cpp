/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_backmap_pairs.h"

#include <cmath>
#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "text_file_reader.h"
#include "tokenizer.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixBackmapPairs::FixBackmapPairs(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg),
      cut(0.0),
      cutsq(0.0),
      has_coulomb(false),
      energy_local{0.0, 0.0},
      fix_backmap(nullptr) {
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix backmap/pairs", error);

  // The 1-4 terms are part of the potential energy and the virial, as for
  // the pair style they replace (special_bonds excludes them there).
  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  global_freq = 1;
  extscalar = 1;
  extvector = 1;
  energy_global_flag = 1;
  virial_global_flag = virial_peratom_flag = 1;
  thermo_energy = thermo_virial = 1;

  int iarg = 3;
  if (strcmp(arg[iarg], "at") != 0 && strcmp(arg[iarg], "cg") != 0)
    error->all(FLERR, "fix backmap/pairs keyword must be at or cg");
  int default_is_cg = (strcmp(arg[iarg], "cg") == 0) ? 1 : 0;
  iarg++;

  if (strcmp(arg[iarg], "file") != 0)
    error->all(FLERR, "Expected 'file' after fix backmap/pairs keyword");
  iarg++;
  if (iarg >= narg) error->all(FLERR, "Missing filename for fix backmap/pairs");
  read_file(arg[iarg]);
  iarg++;

  // 1-4 pairs are bonded terms: no cutoff unless one is asked for.
  cut = 0.0;
  if (iarg < narg && strcmp(arg[iarg], "cut") == 0) {
    if (iarg + 1 >= narg)
      error->all(FLERR, "Missing cutoff for fix backmap/pairs cut");
    cut = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    iarg += 2;
  }
  cutsq = (cut > 0.0) ? cut * cut : -1.0;

  for (auto &entry : pairs) {
    entry.is_cg = default_is_cg;
    if (entry.qq_scale != 0.0) has_coulomb = true;
  }
}

FixBackmapPairs::~FixBackmapPairs() {}

/* ---------------------------------------------------------------------- */

int FixBackmapPairs::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmapPairs::init() {
  fix_backmap = BackmapLambda::find_fix_backmap(lmp, "fix backmap/pairs");
  if (has_coulomb && !atom->q_flag)
    error->all(
        FLERR,
        "fix backmap/pairs: the pairs file has 1-4 Coulomb scales but the "
        "atom style has no charges");
}

/* ---------------------------------------------------------------------- */

void FixBackmapPairs::setup(int vflag) { post_force(vflag); }

void FixBackmapPairs::min_setup(int vflag) { post_force(vflag); }

void FixBackmapPairs::min_post_force(int vflag) { post_force(vflag); }

/* ----------------------------------------------------------------------
   pairs file: first line N, then N lines
     id1 id2 sigma epsilon [qq_scale]
   qq_scale (optional, default 0) scales the 1-4 Coulomb term q_i q_j / r.
------------------------------------------------------------------------- */

void FixBackmapPairs::read_file(const char *filename) {
  std::string filecontent = utils::get_potential_file_path(filename);
  if (filecontent.empty())
    error->all(FLERR, "Cannot open backmap pairs file {}", filename);

  auto reader = TextFileReader(filecontent, "backmap pairs");
  reader.ignore_comments = true;

  char *line = reader.next_line();
  if (!line) error->all(FLERR, "Empty backmap pairs file {}", filename);
  ValueTokenizer header(line);
  int npairs = header.next_int();
  if (npairs <= 0)
    error->all(FLERR, "Invalid pair count in backmap pairs file {}", filename);

  pairs.clear();
  pairs.reserve(npairs);

  for (int n = 0; n < npairs; n++) {
    line = reader.next_line();
    if (!line)
      error->all(FLERR, "Premature end of backmap pairs file {}", filename);
    ValueTokenizer values(line);
    int ncols = static_cast<int>(values.count());
    if (ncols != 4 && ncols != 5)
      error->all(FLERR,
                 "Backmap pairs file {}: expected 4 or 5 columns, got {}",
                 filename, ncols);
    PairEntry entry;
    entry.id1 = values.next_tagint();
    entry.id2 = values.next_tagint();
    entry.sigma = values.next_double();
    entry.epsilon = values.next_double();
    entry.qq_scale = (ncols == 5) ? values.next_double() : 0.0;
    entry.is_cg = 0;
    pairs.push_back(entry);
  }
}

/* ----------------------------------------------------------------------
   Each rank evaluates every pair with at least one local atom and adds the
   force to its local atoms only, so no reverse communication is needed
   (post_force runs after it). Energy and virial are split evenly between the
   two atoms, so a pair spanning two ranks is counted once in total.
------------------------------------------------------------------------- */

void FixBackmapPairs::post_force(int vflag) {
  if (!fix_backmap)
    fix_backmap = BackmapLambda::find_fix_backmap(lmp, "fix backmap/pairs");

  int *atom2cg = BackmapLambda::extract_atom2cg(fix_backmap);
  double *lam_global_ptr = BackmapLambda::extract_lambda_global(fix_backmap);
  if (!atom2cg || !lam_global_ptr)
    error->all(FLERR,
               "fix backmap/pairs: cannot extract atom2cg/lambda_global "
               "from fix backmap");
  double lambda_global = BackmapLambda::clamp_lambda(*lam_global_ptr);

  energy_local[0] = energy_local[1] = 0.0;
  // Global energy is accumulated unconditionally (compute_scalar/vector);
  // ev_init/ev_tally handle the virial only (no per-atom energy).
  ev_init(0, vflag);

  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  int nlocal = atom->nlocal;
  double qqrd2e = force->qqrd2e;

  for (const auto &par : pairs) {
    // atom->map returns the owned (local) index when this rank owns the atom.
    int i = atom->map(par.id1);
    int j = atom->map(par.id2);
    bool own_i = i >= 0 && i < nlocal;
    bool own_j = j >= 0 && j < nlocal;
    if (!own_i && !own_j) continue;
    if (i < 0 || j < 0)
      error->one(
          FLERR,
          "fix backmap/pairs: partner of 1-4 pair {} {} missing on this proc; "
          "increase comm_modify cutoff",
          par.id1, par.id2);

    // Geometry from the closest image of the partner; forces go to the owned
    // atoms themselves. (Using the image index for the force would drop the
    // force on an owned partner whose closest image is a periodic ghost.)
    double dx, dy, dz;
    if (own_i) {
      int jj = domain->closest_image(i, j);
      dx = x[i][0] - x[jj][0];
      dy = x[i][1] - x[jj][1];
      dz = x[i][2] - x[jj][2];
    } else {
      int ii = domain->closest_image(j, i);
      dx = x[ii][0] - x[j][0];
      dy = x[ii][1] - x[j][1];
      dz = x[ii][2] - x[j][2];
    }
    double rsq = dx * dx + dy * dy + dz * dz;
    if ((cutsq > 0.0 && rsq >= cutsq) || rsq <= 0.0) continue;

    bool same_bead = BackmapLambda::same_bead(atom2cg, i, j);
    double w =
        BackmapLambda::compute_weight3(same_bead, par.is_cg, lambda_global);
    if (BackmapLambda::is_almost_zero(w)) continue;

    double r2inv = 1.0 / rsq;
    double r6inv = r2inv * r2inv * r2inv;
    double sig6 = std::pow(par.sigma, 6.0);
    double lj6 = par.epsilon * sig6 * r6inv;
    double lj12 = lj6 * sig6 * r6inv;
    double e_lj = 4.0 * (lj12 - lj6);
    double f_lj = 24.0 * (2.0 * lj12 - lj6) * r2inv;

    double e_coul = 0.0, f_coul = 0.0;
    if (par.qq_scale != 0.0) {
      e_coul = qqrd2e * par.qq_scale * q[i] * q[j] / sqrt(rsq);
      f_coul = e_coul * r2inv;
    }
    double fpair = w * (f_lj + f_coul);

    int list[2];
    int nlist = 0;
    if (own_i) {
      f[i][0] += dx * fpair;
      f[i][1] += dy * fpair;
      f[i][2] += dz * fpair;
      list[nlist++] = i;
    }
    if (own_j) {
      f[j][0] -= dx * fpair;
      f[j][1] -= dy * fpair;
      f[j][2] -= dz * fpair;
      list[nlist++] = j;
    }

    double share = 0.5 * nlist;
    energy_local[0] += share * w * e_lj;
    energy_local[1] += share * w * e_coul;

    if (evflag) {
      double v[6] = {dx * dx * fpair, dy * dy * fpair, dz * dz * fpair,
                     dx * dy * fpair, dx * dz * fpair, dy * dz * fpair};
      ev_tally(nlist, list, 2.0, w * (e_lj + e_coul), v);
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixBackmapPairs::compute_scalar() {
  double all[2];
  MPI_Allreduce(energy_local, all, 2, MPI_DOUBLE, MPI_SUM, world);
  return all[0] + all[1];
}

/* ---------------------------------------------------------------------- */

double FixBackmapPairs::compute_vector(int n) {
  double all[2];
  MPI_Allreduce(energy_local, all, 2, MPI_DOUBLE, MPI_SUM, world);
  return all[n];
}

/* ---------------------------------------------------------------------- */

double FixBackmapPairs::memory_usage() {
  return static_cast<double>(pairs.size() * sizeof(PairEntry));
}
