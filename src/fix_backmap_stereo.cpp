/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_backmap_stereo.h"

#include <mpi.h>

#include <cmath>
#include <fstream>
#include <string>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "math_const.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using MathConst::DEG2RAD;
using MathConst::MY_PI;

/* ---------------------------------------------------------------------- */

FixBackmapStereo::FixBackmapStereo(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg), kspring(0.0), vector_reduced(false) {
  if (narg != 5)
    error->all(FLERR,
               "Illegal fix backmap/stereo command: fix ID group "
               "backmap/stereo K file");

  kspring = utils::numeric(FLERR, arg[3], false, lmp);
  if (kspring < 0.0) error->all(FLERR, "fix backmap/stereo K must be >= 0");
  read_restraints(arg[4]);

  const int n = static_cast<int>(cos_shift.size());
  energy_local.assign(n, 0.0);
  energy_all.assign(n, 0.0);

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = n;
  extscalar = 1;
  extvector = 1;
  global_freq = 1;
  energy_global_flag = 1;
}

/* ----------------------------------------------------------------------
   rank 0 reads "atom1 atom2 atom3 atom4 phi_t" lines and broadcasts them
------------------------------------------------------------------------- */

void FixBackmapStereo::read_restraints(const char *file) {
  std::vector<double> targets;
  std::string errmsg;

  if (comm->me == 0) {
    std::ifstream in(file);
    if (!in)
      errmsg = fmt::format("fix backmap/stereo cannot open file {}", file);
    std::string line;
    int lineno = 0;
    while (errmsg.empty() && std::getline(in, line)) {
      ++lineno;
      auto words = utils::split_words(utils::trim_comment(line));
      if (words.empty()) continue;
      if (words.size() != 5) {
        errmsg = fmt::format(
            "fix backmap/stereo {} line {}: expected 4 atom IDs and a target "
            "angle",
            file, lineno);
        break;
      }
      bool valid = utils::is_double(words[4]);
      for (int k = 0; k < 4; ++k) valid = valid && utils::is_integer(words[k]);
      if (!valid) {
        errmsg = fmt::format("fix backmap/stereo {} line {}: invalid number",
                             file, lineno);
        break;
      }
      for (int k = 0; k < 4; ++k)
        ids.push_back(static_cast<tagint>(std::stoll(words[k])));
      targets.push_back(std::stod(words[4]));
    }
    if (errmsg.empty() && targets.empty())
      errmsg = fmt::format("fix backmap/stereo {} holds no restraints", file);
  }

  int len = static_cast<int>(errmsg.size());
  MPI_Bcast(&len, 1, MPI_INT, 0, world);
  if (len > 0) {
    errmsg.resize(len);
    MPI_Bcast(&errmsg[0], len, MPI_CHAR, 0, world);
    error->all(FLERR, errmsg);
  }

  int n = static_cast<int>(targets.size());
  MPI_Bcast(&n, 1, MPI_INT, 0, world);
  ids.resize(4 * static_cast<size_t>(n));
  targets.resize(n);
  MPI_Bcast(ids.data(), 4 * n, MPI_LMP_TAGINT, 0, world);
  MPI_Bcast(targets.data(), n, MPI_DOUBLE, 0, world);

  // E = K [1 + cos(phi - d)] with d = phi_t - 180 is K [1 - cos(phi - phi_t)]
  cos_shift.resize(n);
  sin_shift.resize(n);
  for (int m = 0; m < n; ++m) {
    const double d = DEG2RAD * targets[m] - MY_PI;
    cos_shift[m] = cos(d);
    sin_shift[m] = sin(d);
  }
}

/* ---------------------------------------------------------------------- */

int FixBackmapStereo::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmapStereo::init() {
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,
               "fix backmap/stereo requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

void FixBackmapStereo::setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapStereo::min_setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapStereo::min_post_force(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapStereo::post_force(int /*vflag*/) {
  vector_reduced = false;
  const int n = static_cast<int>(cos_shift.size());
  for (int m = 0; m < n; ++m) energy_local[m] = restrain(m);
}

/* ----------------------------------------------------------------------
   one restraint; the force expressions follow fix restrain dihedral
   (adopted from dihedral_charmm) with multiplicity 1
------------------------------------------------------------------------- */

double FixBackmapStereo::restrain(int m) {
  double **x = atom->x;
  double **f = atom->f;
  const int nlocal = atom->nlocal;
  const tagint *id = &ids[4 * static_cast<size_t>(m)];

  const int i1 = atom->map(id[0]);
  const int i2 = atom->map(id[1]);
  const int i3 = atom->map(id[2]);
  const int i4 = atom->map(id[3]);

  // fix forces are added after ghost forces were communicated back to their
  // owners, so every processor owning one of the atoms computes the restraint
  // and applies forces only to its own atoms

  if ((i1 < 0 || i1 >= nlocal) && (i2 < 0 || i2 >= nlocal) &&
      (i3 < 0 || i3 >= nlocal) && (i4 < 0 || i4 >= nlocal))
    return 0.0;
  if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 0)
    error->one(
        FLERR,
        "fix backmap/stereo atoms {} {} {} {} missing on proc {} at step {}",
        id[0], id[1], id[2], id[3], comm->me, update->ntimestep);

  double vb1x = x[i1][0] - x[i2][0];
  double vb1y = x[i1][1] - x[i2][1];
  double vb1z = x[i1][2] - x[i2][2];
  domain->minimum_image(FLERR, vb1x, vb1y, vb1z);

  double vb2x = x[i3][0] - x[i2][0];
  double vb2y = x[i3][1] - x[i2][1];
  double vb2z = x[i3][2] - x[i2][2];
  domain->minimum_image(FLERR, vb2x, vb2y, vb2z);

  double vb2xm = -vb2x;
  double vb2ym = -vb2y;
  double vb2zm = -vb2z;
  domain->minimum_image(FLERR, vb2xm, vb2ym, vb2zm);

  double vb3x = x[i4][0] - x[i3][0];
  double vb3y = x[i4][1] - x[i3][1];
  double vb3z = x[i4][2] - x[i3][2];
  domain->minimum_image(FLERR, vb3x, vb3y, vb3z);

  const double ax = vb1y * vb2zm - vb1z * vb2ym;
  const double ay = vb1z * vb2xm - vb1x * vb2zm;
  const double az = vb1x * vb2ym - vb1y * vb2xm;
  const double bx = vb3y * vb2zm - vb3z * vb2ym;
  const double by = vb3z * vb2xm - vb3x * vb2zm;
  const double bz = vb3x * vb2ym - vb3y * vb2xm;

  const double rasq = ax * ax + ay * ay + az * az;
  const double rbsq = bx * bx + by * by + bz * bz;
  const double rgsq = vb2xm * vb2xm + vb2ym * vb2ym + vb2zm * vb2zm;
  const double rg = sqrt(rgsq);

  double rginv = 0.0, ra2inv = 0.0, rb2inv = 0.0;
  if (rg > 0) rginv = 1.0 / rg;
  if (rasq > 0) ra2inv = 1.0 / rasq;
  if (rbsq > 0) rb2inv = 1.0 / rbsq;
  const double rabinv = sqrt(ra2inv * rb2inv);

  double c = (ax * bx + ay * by + az * bz) * rabinv;
  const double s = rg * rabinv * (ax * vb3x + ay * vb3y + az * vb3z);
  if (c > 1.0) c = 1.0;
  if (c < -1.0) c = -1.0;

  const double p = 1.0 + c * cos_shift[m] + s * sin_shift[m];
  const double df1 = -(s * cos_shift[m] - c * sin_shift[m]);

  const double fg = vb1x * vb2xm + vb1y * vb2ym + vb1z * vb2zm;
  const double hg = vb3x * vb2xm + vb3y * vb2ym + vb3z * vb2zm;
  const double fga = fg * ra2inv * rginv;
  const double hgb = hg * rb2inv * rginv;
  const double gaa = -ra2inv * rg;
  const double gbb = rb2inv * rg;

  const double dtfx = gaa * ax, dtfy = gaa * ay, dtfz = gaa * az;
  const double dtgx = fga * ax - hgb * bx;
  const double dtgy = fga * ay - hgb * by;
  const double dtgz = fga * az - hgb * bz;
  const double dthx = gbb * bx, dthy = gbb * by, dthz = gbb * bz;

  const double df = -kspring * df1;

  const double sx2 = df * dtgx, sy2 = df * dtgy, sz2 = df * dtgz;
  const double f1[3] = {df * dtfx, df * dtfy, df * dtfz};
  const double f2[3] = {sx2 - f1[0], sy2 - f1[1], sz2 - f1[2]};
  const double f4[3] = {df * dthx, df * dthy, df * dthz};
  const double f3[3] = {-sx2 - f4[0], -sy2 - f4[1], -sz2 - f4[2]};

  const int idx[4] = {i1, i2, i3, i4};
  const double *fi[4] = {f1, f2, f3, f4};
  int nown = 0;
  for (int k = 0; k < 4; ++k) {
    if (idx[k] >= nlocal) continue;
    ++nown;
    f[idx[k]][0] += fi[k][0];
    f[idx[k]][1] += fi[k][1];
    f[idx[k]][2] += fi[k][2];
  }

  // each processor tallies the energy share of its own atoms
  return 0.25 * nown * kspring * p;
}

/* ---------------------------------------------------------------------- */

double FixBackmapStereo::compute_scalar() {
  compute_vector(0);
  double total = 0.0;
  for (double e : energy_all) total += e;
  return total;
}

/* ---------------------------------------------------------------------- */

double FixBackmapStereo::compute_vector(int n) {
  if (!vector_reduced) {
    MPI_Allreduce(energy_local.data(), energy_all.data(),
                  static_cast<int>(energy_local.size()), MPI_DOUBLE, MPI_SUM,
                  world);
    vector_reduced = true;
  }
  return energy_all[n];
}
