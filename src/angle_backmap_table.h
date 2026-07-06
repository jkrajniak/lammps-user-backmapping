/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ANGLE_CLASS
// clang-format off
AngleStyle(backmap/table,AngleBackmapTable);
// clang-format on
#else

#ifndef LMP_ANGLE_BACKMAP_TABLE_H
#define LMP_ANGLE_BACKMAP_TABLE_H

#include "angle.h"

namespace LAMMPS_NS {

class AngleBackmapTable : public Angle {
 public:
  AngleBackmapTable(class LAMMPS *);
  ~AngleBackmapTable() override;
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double equilibrium_angle(int) override;
  void write_restart(FILE *) override;
  void read_restart(FILE *) override;
  double single(int, int, int, int) override;

 protected:
  int tabstyle, tablength;

  struct Table {
    int ninput;
    double lo, hi;
    double *rfile, *efile, *ffile;
    double *e2file, *f2file;
    double delta, invdelta;
    double *r, *e, *f;
    double *de, *df;
  };

  int ntables;
  Table *tables;
  int *tabindex;
  int *is_cg;

  class Fix *fix_backmap;

  void null_table(Table *);
  void free_table(Table *);
  void read_table(Table *, const char *, const char *);
  void bcast_table(Table *);
  void spline_table(Table *);
  void compute_table(Table *);
  double uf_lookup(int, double, double &);

  virtual void allocate();
};

}  // namespace LAMMPS_NS

#endif
#endif
