/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* dihedral_style backmap/opls -- lambda-weighted OPLS dihedral.

   E = w [K1/2 (1 + cos phi) + K2/2 (1 - cos 2phi)
          + K3/2 (1 + cos 3phi) + K4/2 (1 - cos 4phi)]

   Same form and phi convention (trans = 180 deg) as LAMMPS dihedral_style
   opls. The series is exact in cos(phi), so coeff() converts K1..K4 once to
   the backmap/ryckaert coefficients E = w sum_n Cn cos^n(phi):
     C0 = K1/2 + K2 + K3/2   C1 = K1/2 - 3 K3/2   C2 = -K2 + 4 K4
     C3 = 2 K3               C4 = -4 K4           C5 = 0
   and the RB kernel computes energy and forces. w as for every backmap
   style: 1 intra-bead, lambda for at terms across beads, 1 - lambda for cg.

   Syntax:
     dihedral_style backmap/opls
     dihedral_coeff N at/cg K1 K2 K3 K4
------------------------------------------------------------------------- */

#include "dihedral_backmap_opls.h"

#include <cstring>

#include "atom.h"
#include "backmap_lambda.h"
#include "comm.h"
#include "error.h"
#include "memory.h"
#include "utils.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

DihedralBackmapOpls::DihedralBackmapOpls(LAMMPS *lmp)
    : DihedralBackmapRyckaert(lmp),
      k1(nullptr),
      k2(nullptr),
      k3(nullptr),
      k4(nullptr) {}

/* ---------------------------------------------------------------------- */

DihedralBackmapOpls::~DihedralBackmapOpls() {
  if (allocated && !copymode) {
    memory->destroy(k1);
    memory->destroy(k2);
    memory->destroy(k3);
    memory->destroy(k4);
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapOpls::allocate() {
  DihedralBackmapRyckaert::allocate();
  int n = atom->ndihedraltypes + 1;
  memory->create(k1, n, "dihedral:k1");
  memory->create(k2, n, "dihedral:k2");
  memory->create(k3, n, "dihedral:k3");
  memory->create(k4, n, "dihedral:k4");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapOpls::set_rb(int i) {
  c0[i] = 0.5 * k1[i] + k2[i] + 0.5 * k3[i];
  c1[i] = 0.5 * k1[i] - 1.5 * k3[i];
  c2[i] = -k2[i] + 4.0 * k4[i];
  c3[i] = 2.0 * k3[i];
  c4[i] = -4.0 * k4[i];
  c5[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapOpls::coeff(int narg, char **arg) {
  if (narg != 6)
    error->all(FLERR, "Incorrect args for dihedral_coeff backmap/opls");
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->ndihedraltypes, ilo, ihi, error);

  int cg_flag;
  if (strcmp(arg[1], "at") == 0)
    cg_flag = 0;
  else if (strcmp(arg[1], "cg") == 0)
    cg_flag = 1;
  else
    error->all(FLERR,
               "dihedral_coeff backmap/opls: 2nd arg must be 'at' or 'cg'");

  double k1_one = utils::numeric(FLERR, arg[2], false, lmp);
  double k2_one = utils::numeric(FLERR, arg[3], false, lmp);
  double k3_one = utils::numeric(FLERR, arg[4], false, lmp);
  double k4_one = utils::numeric(FLERR, arg[5], false, lmp);

  for (int i = ilo; i <= ihi; i++) {
    k1[i] = k1_one;
    k2[i] = k2_one;
    k3[i] = k3_one;
    k4[i] = k4_one;
    set_rb(i);
    is_cg[i] = cg_flag;
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapOpls::init_style() {
  fix_backmap =
      BackmapLambda::find_fix_backmap(lmp, "dihedral_style backmap/opls");
}

/* ---------------------------------------------------------------------- */

void DihedralBackmapOpls::write_restart(FILE *fp) {
  fwrite(k1 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(k2 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(k3 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(k4 + 1, sizeof(double), atom->ndihedraltypes, fp);
  fwrite(is_cg + 1, sizeof(int), atom->ndihedraltypes, fp);
}

void DihedralBackmapOpls::read_restart(FILE *fp) {
  allocate();
  if (comm->me == 0) {
    utils::sfread(FLERR, k1 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, k2 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, k3 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, k4 + 1, sizeof(double), atom->ndihedraltypes, fp,
                  nullptr, error);
    utils::sfread(FLERR, is_cg + 1, sizeof(int), atom->ndihedraltypes, fp,
                  nullptr, error);
  }
  MPI_Bcast(k1 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(k2 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(k3 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(k4 + 1, atom->ndihedraltypes, MPI_DOUBLE, 0, world);
  MPI_Bcast(is_cg + 1, atom->ndihedraltypes, MPI_INT, 0, world);

  for (int i = 1; i <= atom->ndihedraltypes; i++) {
    set_rb(i);
    setflag[i] = 1;
  }
}

/* ---------------------------------------------------------------------- */

/* proc 0 writes to data file, one line per type in coeff() argument order */

void DihedralBackmapOpls::write_data(FILE *fp) {
  for (int i = 1; i <= atom->ndihedraltypes; i++)
    fprintf(fp, "%d %s %.15g %.15g %.15g %.15g\n", i, is_cg[i] ? "cg" : "at",
            k1[i], k2[i], k3[i], k4[i]);
}
