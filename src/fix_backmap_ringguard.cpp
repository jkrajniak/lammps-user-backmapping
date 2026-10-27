/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_backmap_ringguard.h"

#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "memory.h"
#include "update.h"
#include "utils.h"

using namespace LAMMPS_NS;
using namespace FixConst;

static constexpr int MIN_RING = 3;
static constexpr int MAX_RING = 12;
static constexpr int MAX_EXCLUDED = 64;

/* ---------------------------------------------------------------------- */

FixBackmapRingGuard::FixBackmapRingGuard(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg),
      kwall(0.0),
      rc(0.0),
      fbuf(nullptr),
      nmax(0),
      energy_local(0.0),
      contacts_local(0.0),
      sums{0.0, 0.0},
      reduced(false),
      cell_lo{0.0, 0.0, 0.0},
      cell_n{0, 0, 0} {
  if (narg != 6)
    error->all(FLERR,
               "Illegal fix backmap/ringguard command: fix ID group "
               "backmap/ringguard K rc file");

  kwall = utils::numeric(FLERR, arg[3], false, lmp);
  rc = utils::numeric(FLERR, arg[4], false, lmp);
  if (kwall < 0.0) error->all(FLERR, "fix backmap/ringguard K must be >= 0");
  if (rc <= 0.0) error->all(FLERR, "fix backmap/ringguard rc must be > 0");
  read_rings(arg[5]);

  scalar_flag = 1;
  vector_flag = 1;
  size_vector = 2;
  extscalar = 1;
  extvector = 1;
  global_freq = 1;
  energy_global_flag = 1;
  comm_reverse = 3;
}

/* ---------------------------------------------------------------------- */

FixBackmapRingGuard::~FixBackmapRingGuard() { memory->destroy(fbuf); }

/* ----------------------------------------------------------------------
   rank 0 reads one ring per line (atom IDs, then optionally | and the IDs of
   atoms the ring does not repel) and broadcasts them
------------------------------------------------------------------------- */

void FixBackmapRingGuard::read_rings(const char *file) {
  std::string errmsg;
  ring_start.assign(1, 0);
  excl_start.assign(1, 0);

  if (comm->me == 0) {
    std::ifstream in(file);
    if (!in)
      errmsg = fmt::format("fix backmap/ringguard cannot open file {}", file);
    std::string line;
    int lineno = 0;
    while (errmsg.empty() && std::getline(in, line)) {
      ++lineno;
      auto words = utils::split_words(utils::trim_comment(line));
      if (words.empty()) continue;
      const auto bar = std::find(words.begin(), words.end(), "|");
      const int n = static_cast<int>(bar - words.begin());
      const int nexcl =
          bar == words.end() ? 0 : static_cast<int>(words.end() - bar) - 1;
      if (n < MIN_RING || n > MAX_RING) {
        errmsg = fmt::format(
            "fix backmap/ringguard {} line {}: a ring needs {} to {} atom IDs",
            file, lineno, MIN_RING, MAX_RING);
        break;
      }
      if (nexcl > MAX_EXCLUDED) {
        errmsg = fmt::format(
            "fix backmap/ringguard {} line {}: more than {} excluded atom IDs",
            file, lineno, MAX_EXCLUDED);
        break;
      }
      for (auto it = words.begin(); it != words.end(); ++it) {
        if (it == bar) continue;
        if (!utils::is_integer(*it)) {
          errmsg =
              fmt::format("fix backmap/ringguard {} line {}: invalid atom ID",
                          file, lineno);
          break;
        }
        const auto id = static_cast<tagint>(std::stoll(*it));
        if (it < bar)
          ring_atoms.push_back(id);
        else
          excl_atoms.push_back(id);
      }
      ring_start.push_back(static_cast<int>(ring_atoms.size()));
      excl_start.push_back(static_cast<int>(excl_atoms.size()));
    }
    if (errmsg.empty() && ring_start.size() < 2)
      errmsg = fmt::format("fix backmap/ringguard {} holds no rings", file);
  }

  int len = static_cast<int>(errmsg.size());
  MPI_Bcast(&len, 1, MPI_INT, 0, world);
  if (len > 0) {
    errmsg.resize(len);
    MPI_Bcast(&errmsg[0], len, MPI_CHAR, 0, world);
    error->all(FLERR, errmsg);
  }

  int nrings = static_cast<int>(ring_start.size()) - 1;
  int natoms = static_cast<int>(ring_atoms.size());
  int nexcl = static_cast<int>(excl_atoms.size());
  MPI_Bcast(&nrings, 1, MPI_INT, 0, world);
  MPI_Bcast(&natoms, 1, MPI_INT, 0, world);
  MPI_Bcast(&nexcl, 1, MPI_INT, 0, world);
  ring_start.resize(nrings + 1);
  ring_atoms.resize(natoms);
  excl_start.resize(nrings + 1);
  excl_atoms.resize(nexcl);
  MPI_Bcast(ring_start.data(), nrings + 1, MPI_INT, 0, world);
  MPI_Bcast(ring_atoms.data(), natoms, MPI_LMP_TAGINT, 0, world);
  MPI_Bcast(excl_start.data(), nrings + 1, MPI_INT, 0, world);
  if (nexcl > 0) MPI_Bcast(excl_atoms.data(), nexcl, MPI_LMP_TAGINT, 0, world);
}

/* ---------------------------------------------------------------------- */

int FixBackmapRingGuard::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::init() {
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR,
               "fix backmap/ringguard requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::min_setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::min_post_force(int vflag) { post_force(vflag); }

/* ----------------------------------------------------------------------
   linked-list cell grid with cell size rc over owned and ghost group atoms
------------------------------------------------------------------------- */

void FixBackmapRingGuard::build_cells() {
  double **x = atom->x;
  int *mask = atom->mask;
  const int nall = atom->nlocal + atom->nghost;

  double lo[3] = {1.0e300, 1.0e300, 1.0e300};
  double hi[3] = {-1.0e300, -1.0e300, -1.0e300};
  for (int i = 0; i < nall; ++i) {
    if (!(mask[i] & groupbit)) continue;
    for (int k = 0; k < 3; ++k) {
      lo[k] = std::min(lo[k], x[i][k]);
      hi[k] = std::max(hi[k], x[i][k]);
    }
  }
  long ncell = 1;
  for (int k = 0; k < 3; ++k) {
    if (hi[k] < lo[k]) hi[k] = lo[k] = 0.0;
    cell_lo[k] = lo[k];
    cell_n[k] = static_cast<int>((hi[k] - lo[k]) / rc) + 1;
    ncell *= cell_n[k];
  }
  cell_head.assign(ncell, -1);
  cell_next.assign(nall, -1);
  for (int i = 0; i < nall; ++i) {
    if (!(mask[i] & groupbit)) continue;
    int c[3];
    for (int k = 0; k < 3; ++k)
      c[k] = std::min(static_cast<int>((x[i][k] - cell_lo[k]) / rc),
                      cell_n[k] - 1);
    const long id =
        (static_cast<long>(c[2]) * cell_n[1] + c[1]) * cell_n[0] + c[0];
    cell_next[i] = cell_head[id];
    cell_head[id] = i;
  }
}

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::post_force(int /*vflag*/) {
  reduced = false;
  energy_local = contacts_local = 0.0;

  if (atom->nmax > nmax) {
    nmax = atom->nmax;
    memory->destroy(fbuf);
    memory->create(fbuf, nmax, 3, "backmap/ringguard:fbuf");
  }
  const int nlocal = atom->nlocal;
  const int nall = nlocal + atom->nghost;
  for (int i = 0; i < nall; ++i) fbuf[i][0] = fbuf[i][1] = fbuf[i][2] = 0.0;

  build_cells();

  double **x = atom->x;
  tagint *tag = atom->tag;
  const double rcsq = rc * rc;
  const int nrings = static_cast<int>(ring_start.size()) - 1;
  int idx[MAX_RING];

  for (int m = 0; m < nrings; ++m) {
    const int first = ring_start[m];
    const int n = ring_start[m + 1] - first;
    const tagint *ring = &ring_atoms[first];
    const tagint *excl = excl_atoms.data() + excl_start[m];
    const int nexcl = excl_start[m + 1] - excl_start[m];

    // the owner of the first ring atom evaluates the ring
    const int anchor = atom->map(ring[0]);
    if (anchor < 0 || anchor >= nlocal) continue;

    double cen[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < n; ++k) {
      const int j = atom->map(ring[k]);
      if (j < 0)
        error->one(
            FLERR,
            "fix backmap/ringguard ring atom {} missing on proc {} at step {}",
            ring[k], comm->me, update->ntimestep);
      idx[k] = domain->closest_image(anchor, j);
      cen[0] += x[idx[k]][0];
      cen[1] += x[idx[k]][1];
      cen[2] += x[idx[k]][2];
    }
    cen[0] /= n;
    cen[1] /= n;
    cen[2] /= n;

    int c[3];
    for (int k = 0; k < 3; ++k)
      c[k] = static_cast<int>(std::floor((cen[k] - cell_lo[k]) / rc));

    for (int cz = c[2] - 1; cz <= c[2] + 1; ++cz) {
      if (cz < 0 || cz >= cell_n[2]) continue;
      for (int cy = c[1] - 1; cy <= c[1] + 1; ++cy) {
        if (cy < 0 || cy >= cell_n[1]) continue;
        for (int cx = c[0] - 1; cx <= c[0] + 1; ++cx) {
          if (cx < 0 || cx >= cell_n[0]) continue;
          const long id =
              (static_cast<long>(cz) * cell_n[1] + cy) * cell_n[0] + cx;
          for (int j = cell_head[id]; j >= 0; j = cell_next[j]) {
            const double dx = x[j][0] - cen[0];
            const double dy = x[j][1] - cen[1];
            const double dz = x[j][2] - cen[2];
            const double rsq = dx * dx + dy * dy + dz * dz;
            if (rsq >= rcsq || rsq == 0.0) continue;
            bool member = false;
            for (int k = 0; k < n && !member; ++k) member = (tag[j] == ring[k]);
            for (int k = 0; k < nexcl && !member; ++k)
              member = (tag[j] == excl[k]);
            if (member) continue;

            const double r = sqrt(rsq);
            const double gap = rc - r;
            energy_local += kwall * gap * gap;
            contacts_local += 1.0;
            // force on the atom along the centroid-atom axis, reaction shared
            // equally by the ring atoms
            const double fpair = 2.0 * kwall * gap / r;
            fbuf[j][0] += fpair * dx;
            fbuf[j][1] += fpair * dy;
            fbuf[j][2] += fpair * dz;
            for (int k = 0; k < n; ++k) {
              fbuf[idx[k]][0] -= fpair * dx / n;
              fbuf[idx[k]][1] -= fpair * dy / n;
              fbuf[idx[k]][2] -= fpair * dz / n;
            }
          }
        }
      }
    }
  }

  // forces added after the reverse communication of the force computation:
  // sum the ghost contributions onto their owners here
  comm->reverse_comm(this);

  double **f = atom->f;
  for (int i = 0; i < nlocal; ++i) {
    f[i][0] += fbuf[i][0];
    f[i][1] += fbuf[i][1];
    f[i][2] += fbuf[i][2];
  }
}

/* ---------------------------------------------------------------------- */

int FixBackmapRingGuard::pack_reverse_comm(int n, int first, double *buf) {
  int m = 0;
  for (int i = first; i < first + n; ++i) {
    buf[m++] = fbuf[i][0];
    buf[m++] = fbuf[i][1];
    buf[m++] = fbuf[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::unpack_reverse_comm(int n, int *list, double *buf) {
  int m = 0;
  for (int i = 0; i < n; ++i) {
    const int j = list[i];
    fbuf[j][0] += buf[m++];
    fbuf[j][1] += buf[m++];
    fbuf[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

double FixBackmapRingGuard::compute_scalar() { return compute_vector(0); }

/* ---------------------------------------------------------------------- */

double FixBackmapRingGuard::compute_vector(int n) {
  if (!reduced) {
    const double local[2] = {energy_local, contacts_local};
    MPI_Allreduce(local, sums, 2, MPI_DOUBLE, MPI_SUM, world);
    reduced = true;
  }
  return sums[n];
}

/* ---------------------------------------------------------------------- */

double FixBackmapRingGuard::memory_usage() {
  return 3.0 * nmax * sizeof(double) +
         (cell_head.capacity() + cell_next.capacity()) * sizeof(int) +
         (ring_atoms.capacity() + excl_atoms.capacity()) * sizeof(tagint);
}
