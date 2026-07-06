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
    : Fix(lmp, narg, arg), cut(0.0), cutsq(0.0), fix_backmap(nullptr) {
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix backmap/pairs", error);

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

  cut = force->pair ? force->pair->cutforce : 12.0;
  if (iarg < narg && strcmp(arg[iarg], "cut") == 0) {
    if (iarg + 1 >= narg)
      error->all(FLERR, "Missing cutoff for fix backmap/pairs cut");
    cut = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
    iarg += 2;
  }
  cutsq = cut * cut;

  for (auto &entry : pairs) entry.is_cg = default_is_cg;
}

FixBackmapPairs::~FixBackmapPairs() {}

/* ---------------------------------------------------------------------- */

int FixBackmapPairs::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmapPairs::init() {
  fix_backmap = BackmapLambda::find_fix_backmap(lmp, "fix backmap/pairs");
}

/* ---------------------------------------------------------------------- */

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
    PairEntry entry;
    entry.id1 = values.next_tagint();
    entry.id2 = values.next_tagint();
    entry.sigma = values.next_double();
    entry.epsilon = values.next_double();
    entry.is_cg = 0;
    pairs.push_back(entry);
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmapPairs::post_force(int /*vflag*/) {
  if (!fix_backmap)
    fix_backmap = BackmapLambda::find_fix_backmap(lmp, "fix backmap/pairs");

  double *lambda = BackmapLambda::extract_lambda(fix_backmap);
  if (!lambda)
    error->all(FLERR,
               "fix backmap/pairs: cannot extract lambda from fix backmap");

  double **x = atom->x;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  for (const auto &par : pairs) {
    int i = atom->map(par.id1);
    int j = atom->map(par.id2);
    if (i < 0 || j < 0) continue;
    if (i >= nlocal && j >= nlocal) continue;

    if (newton_pair) {
      if ((i >= nlocal) && ((par.id1 + par.id2) & 1) == 0) continue;
      if ((j >= nlocal) && ((par.id1 + par.id2) & 1) == 1) continue;
    }

    double dx = x[i][0] - x[j][0];
    double dy = x[i][1] - x[j][1];
    double dz = x[i][2] - x[j][2];
    domain->minimum_image(FLERR, dx, dy, dz);
    double rsq = dx * dx + dy * dy + dz * dz;
    if (rsq >= cutsq || rsq <= 0.0) continue;

    double li = BackmapLambda::clamp_lambda(lambda[i]);
    double lj = BackmapLambda::clamp_lambda(lambda[j]);
    double w = BackmapLambda::compute_weight(li, lj, par.is_cg);
    if (BackmapLambda::is_almost_zero(w)) continue;

    double r2inv = 1.0 / rsq;
    double r6inv = r2inv * r2inv * r2inv;
    double sig6 = std::pow(par.sigma, 6.0);
    double fpair = w * 24.0 * par.epsilon * r6inv *
                   (2.0 * sig6 * sig6 * r6inv - sig6) * r2inv;

    if (newton_pair || i < nlocal) {
      f[i][0] += dx * fpair;
      f[i][1] += dy * fpair;
      f[i][2] += dz * fpair;
    }
    if (newton_pair || j < nlocal) {
      f[j][0] -= dx * fpair;
      f[j][1] -= dy * fpair;
      f[j][2] -= dz * fpair;
    }
  }
}

/* ---------------------------------------------------------------------- */

double FixBackmapPairs::memory_usage() {
  return static_cast<double>(pairs.size() * sizeof(PairEntry));
}
