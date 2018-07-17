/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef BOND_CLASS

BondStyle(gromos,BondGromos)

#else

#ifndef LMP_BOND_GROMOS_H
#define LMP_BOND_GROMOS_H

#include <cstdio>
#include "bond.h"

namespace LAMMPS_NS {

class BondGromos : public Bond {
 public:
  BondGromos(class LAMMPS *);
  virtual ~BondGromos();
  virtual void compute(int, int);
  void coeff(int, char **);
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_data(FILE *, char **);
  double single(int, double, int, int, double &);
  virtual void *extract(char *, int &);

 protected:
  double *k,*r0;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Incorrect args for bond coefficients

Self-explanatory.  Check the input script or data file.

*/
