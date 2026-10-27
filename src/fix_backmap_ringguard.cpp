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
#include "force.h"
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
      rcb(0.0),
      fbuf(nullptr),
      nmax(0),
      energy_local(0.0),
      contacts_local(0.0),
      sums{0.0, 0.0},
      reduced(false),
      cell_lo{0.0, 0.0, 0.0},
      cell_w{1.0, 1.0, 1.0},
      cell_n{1, 1, 1} {
  if (narg != 7)
    error->all(FLERR,
               "Illegal fix backmap/ringguard command: fix ID group "
               "backmap/ringguard K rc rcb file");

  kwall = utils::numeric(FLERR, arg[3], false, lmp);
  rc = utils::numeric(FLERR, arg[4], false, lmp);
  rcb = utils::numeric(FLERR, arg[5], false, lmp);
  if (kwall < 0.0) error->all(FLERR, "fix backmap/ringguard K must be >= 0");
  if (rc <= 0.0) error->all(FLERR, "fix backmap/ringguard rc must be > 0");
  if (rcb < 0.0 || rcb > rc)
    error->all(FLERR, "fix backmap/ringguard rcb must be in [0, rc]");
  read_rings(arg[6]);

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
   the atoms bonded to the ring, which get the radius rcb) and broadcasts them
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
  if (domain->triclinic)
    error->all(FLERR, "fix backmap/ringguard needs an orthogonal box");
}

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::min_setup(int vflag) { post_force(vflag); }

/* ---------------------------------------------------------------------- */

void FixBackmapRingGuard::min_post_force(int vflag) { post_force(vflag); }

/* ----------------------------------------------------------------------
   centroids of the rings whose atoms are all present on this proc (owned or
   ghost), binned on a periodic grid of cell size >= rc over the global box
------------------------------------------------------------------------- */

void FixBackmapRingGuard::build_rings() {
  double **x = atom->x;
  const int nrings = static_cast<int>(ring_start.size()) - 1;
  ring_idx.resize(ring_atoms.size());
  ring_cen.resize(3 * nrings);

  for (int k = 0; k < 3; ++k) {
    cell_lo[k] = domain->boxlo[k];
    cell_n[k] = std::max(1, static_cast<int>(domain->prd[k] / rc));
    cell_w[k] = domain->prd[k] / cell_n[k];
  }
  cell_head.assign(static_cast<long>(cell_n[0]) * cell_n[1] * cell_n[2], -1);
  cell_next.assign(nrings, -1);

  for (int m = 0; m < nrings; ++m) {
    const int first = ring_start[m];
    const int n = ring_start[m + 1] - first;
    const int anchor = ring_anchor(ring_atoms[first]);
    if (anchor < 0) continue;
    double *cen = &ring_cen[3 * m];
    cen[0] = cen[1] = cen[2] = 0.0;
    bool present = true;
    for (int k = 0; k < n && present; ++k) {
      const int j = atom->map(ring_atoms[first + k]);
      if (j < 0) {
        present = false;
        break;
      }
      const int img = domain->closest_image(anchor, j);
      // a ghost anchor at the edge of the ghost shell may only see images of
      // the other ring atoms across the box: then the ring is not here
      double del[3], mi[3];
      for (int d = 0; d < 3; ++d) mi[d] = del[d] = x[img][d] - x[anchor][d];
      domain->minimum_image(FLERR, mi);
      for (int d = 0; d < 3; ++d)
        if (fabs(mi[d] - del[d]) > 0.25 * domain->prd[d]) present = false;
      ring_idx[first + k] = img;
      for (int d = 0; d < 3; ++d) cen[d] += x[img][d];
    }
    if (!present) continue;
    int c[3];
    for (int d = 0; d < 3; ++d) {
      cen[d] /= n;
      c[d] = cell_of(cen[d], d);
    }
    const long id =
        (static_cast<long>(c[2]) * cell_n[1] + c[1]) * cell_n[0] + c[0];
    cell_next[m] = cell_head[id];
    cell_head[id] = m;
  }
}

/* ----------------------------------------------------------------------
   local image of atom tag t to build its ring around: the owned atom if
   there is one, otherwise the ghost image nearest the centre of the
   subdomain, -1 if the atom is not on this proc
------------------------------------------------------------------------- */

int FixBackmapRingGuard::ring_anchor(tagint t) const {
  int i = atom->map(t);
  if (i < 0) return -1;
  double **x = atom->x;
  int *sametag = atom->sametag;
  const int nlocal = atom->nlocal;
  double mid[3];
  for (int d = 0; d < 3; ++d)
    mid[d] = 0.5 * (domain->sublo[d] + domain->subhi[d]);
  int best = -1;
  double best_r2 = 0.0;
  for (; i >= 0; i = sametag[i]) {
    if (i < nlocal) return i;
    double r2 = 0.0;
    for (int d = 0; d < 3; ++d) r2 += (x[i][d] - mid[d]) * (x[i][d] - mid[d]);
    if (best < 0 || r2 < best_r2) {
      best = i;
      best_r2 = r2;
    }
  }
  return best;
}

/* ----------------------------------------------------------------------
   cell index of coordinate v along dimension d: wrapped when periodic,
   clamped otherwise
------------------------------------------------------------------------- */

int FixBackmapRingGuard::cell_of(double v, int d) const {
  int c = static_cast<int>(std::floor((v - cell_lo[d]) / cell_w[d]));
  if (domain->periodicity[d]) {
    c %= cell_n[d];
    if (c < 0) c += cell_n[d];
  } else {
    c = std::min(std::max(c, 0), cell_n[d] - 1);
  }
  return c;
}

/* ----------------------------------------------------------------------
   rings whose centroid may lie within the box [lo, hi] (absolute coords)
------------------------------------------------------------------------- */

void FixBackmapRingGuard::rings_near(const double *lo, const double *hi,
                                     std::vector<int> &out) const {
  out.clear();
  int c0[3], nc[3];
  for (int d = 0; d < 3; ++d) {
    const int a =
        static_cast<int>(std::floor((lo[d] - cell_lo[d]) / cell_w[d]));
    const int b =
        static_cast<int>(std::floor((hi[d] - cell_lo[d]) / cell_w[d]));
    c0[d] = a;
    nc[d] = std::min(b - a + 1, cell_n[d]);
  }
  for (int iz = 0; iz < nc[2]; ++iz) {
    const int cz = cell_of(cell_lo[2] + (c0[2] + iz + 0.5) * cell_w[2], 2);
    for (int iy = 0; iy < nc[1]; ++iy) {
      const int cy = cell_of(cell_lo[1] + (c0[1] + iy + 0.5) * cell_w[1], 1);
      for (int ix = 0; ix < nc[0]; ++ix) {
        const int cx = cell_of(cell_lo[0] + (c0[0] + ix + 0.5) * cell_w[0], 0);
        const long id =
            (static_cast<long>(cz) * cell_n[1] + cy) * cell_n[0] + cx;
        for (int m = cell_head[id]; m >= 0; m = cell_next[m]) out.push_back(m);
      }
    }
  }
  // a clamped non-periodic range can visit a cell twice
  std::sort(out.begin(), out.end());
  out.erase(std::unique(out.begin(), out.end()), out.end());
}

/* ----------------------------------------------------------------------
   0 = not related to ring m, 1 = ring atom, 2 = atom bonded to the ring
------------------------------------------------------------------------- */

int FixBackmapRingGuard::relation(int m, tagint t) const {
  for (int k = ring_start[m]; k < ring_start[m + 1]; ++k)
    if (ring_atoms[k] == t) return 1;
  for (int k = excl_start[m]; k < excl_start[m + 1]; ++k)
    if (excl_atoms[k] == t) return 2;
  return 0;
}

/* ----------------------------------------------------------------------
   add force f to the ring atoms of ring m, shared equally
------------------------------------------------------------------------- */

void FixBackmapRingGuard::push_ring(int m, const double *f) {
  const int first = ring_start[m];
  const int n = ring_start[m + 1] - first;
  for (int k = 0; k < n; ++k)
    for (int d = 0; d < 3; ++d) fbuf[ring_idx[first + k]][d] += f[d] / n;
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

  build_rings();

  double **x = atom->x;
  tagint *tag = atom->tag;
  int *mask = atom->mask;
  std::vector<int> near;

  // atom walls, evaluated by the owner of the atom
  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit)) continue;
    double lo[3], hi[3];
    for (int d = 0; d < 3; ++d) {
      lo[d] = x[i][d] - rc;
      hi[d] = x[i][d] + rc;
    }
    rings_near(lo, hi, near);
    for (int m : near) {
      const int rel = relation(m, tag[i]);
      if (rel == 1) continue;
      const double rw = rel == 2 ? rcb : rc;
      double del[3];
      for (int d = 0; d < 3; ++d) del[d] = x[i][d] - ring_cen[3 * m + d];
      domain->minimum_image(FLERR, del);
      const double r =
          sqrt(del[0] * del[0] + del[1] * del[1] + del[2] * del[2]);
      if (r >= rw || r == 0.0) continue;
      const double gap = rw - r;
      energy_local += kwall * gap * gap;
      contacts_local += 1.0;
      const double fpair = 2.0 * kwall * gap / r;
      double f[3];
      for (int d = 0; d < 3; ++d) {
        f[d] = fpair * del[d];
        fbuf[i][d] += f[d];
        f[d] = -f[d];
      }
      push_ring(m, f);
    }
  }

  // bond walls: the distance s of the bond segment from the centroid, beyond
  // the atom walls, E = w(s) - w(d_near) with w(x) = K (rw - x)^2 for x < rw
  // and d_near the distance of the nearer bond atom. Zero when the nearest
  // point of the segment is an end, so the energy is continuous. Evaluated
  // once per bond, by the owner of the atom that stores it.
  if (atom->molecular == Atom::MOLECULAR && atom->num_bond) {
    int *num_bond = atom->num_bond;
    int **bond_type = atom->bond_type;
    tagint **bond_atom = atom->bond_atom;
    const int newton_bond = force->newton_bond;
    for (int i = 0; i < nlocal; ++i) {
      if (!(mask[i] & groupbit)) continue;
      for (int b = 0; b < num_bond[i]; ++b) {
        if (bond_type[i][b] <= 0) continue;
        const tagint tj = bond_atom[i][b];
        if (!newton_bond && tag[i] > tj) continue;
        int j = atom->map(tj);
        if (j < 0)
          error->one(FLERR,
                     "fix backmap/ringguard bond atom {} missing on proc {} at "
                     "step {}",
                     tj, comm->me, update->ntimestep);
        j = domain->closest_image(i, j);
        if (!(mask[j] & groupbit)) continue;
        double u[3], lo[3], hi[3];
        for (int d = 0; d < 3; ++d) u[d] = x[j][d] - x[i][d];
        domain->minimum_image(FLERR, u);
        for (int d = 0; d < 3; ++d) {
          lo[d] = x[i][d] + std::min(0.0, u[d]) - rc;
          hi[d] = x[i][d] + std::max(0.0, u[d]) + rc;
        }
        const double uu = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
        if (uu == 0.0) continue;
        rings_near(lo, hi, near);
        for (int m : near) {
          const int ri = relation(m, tag[i]);
          const int rj = relation(m, tj);
          if (ri == 1 || rj == 1) continue;
          const double rw = (ri == 2 || rj == 2) ? rcb : rc;
          // centroid relative to atom i
          double c[3];
          for (int d = 0; d < 3; ++d) c[d] = ring_cen[3 * m + d] - x[i][d];
          domain->minimum_image(FLERR, c);
          const double t = (c[0] * u[0] + c[1] * u[1] + c[2] * u[2]) / uu;
          if (t <= 0.0 || t >= 1.0) continue;
          double p[3], dn[3];
          for (int d = 0; d < 3; ++d) {
            p[d] = t * u[d] - c[d];                 // centroid -> nearest point
            dn[d] = (t < 0.5 ? 0.0 : u[d]) - c[d];  // centroid -> nearer atom
          }
          const double s = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
          if (s >= rw || s == 0.0) continue;
          const double rn = sqrt(dn[0] * dn[0] + dn[1] * dn[1] + dn[2] * dn[2]);
          const double gap = rw - s;
          const double gapn = rn < rw ? rw - rn : 0.0;
          energy_local += kwall * (gap * gap - gapn * gapn);
          contacts_local += 1.0;
          // w(s): pushes the nearest point out, shared by the bond atoms
          const double fs = 2.0 * kwall * gap / s;
          // -w(d_near): pulls the nearer atom back in by its atom wall force
          const double fn = rn > 0.0 ? 2.0 * kwall * gapn / rn : 0.0;
          const int near_atom = t < 0.5 ? i : j;
          double fring[3];
          for (int d = 0; d < 3; ++d) {
            const double f = fs * p[d];
            fbuf[i][d] += (1.0 - t) * f;
            fbuf[j][d] += t * f;
            fbuf[near_atom][d] -= fn * dn[d];
            fring[d] = -f + fn * dn[d];
          }
          push_ring(m, fring);
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
         (cell_head.capacity() + cell_next.capacity() + ring_idx.capacity()) *
             sizeof(int) +
         ring_cen.capacity() * sizeof(double) +
         (ring_atoms.capacity() + excl_atoms.capacity()) * sizeof(tagint);
}
