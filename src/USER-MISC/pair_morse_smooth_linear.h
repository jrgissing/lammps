/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
  -------------------------------------------------------------------------
  Contributed by Stefan Paquay @ Eindhoven University of Technology
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS

PairStyle(morse/smooth/linear,PairMorseSmoothLinear)

#else

#ifndef LMP_PAIR_MORSE_SMOOTH_LINEAR_H
#define LMP_PAIR_MORSE_SMOOTH_LINEAR_H

#include "pair.h"

namespace LAMMPS_NS {

class PairMorseSmoothLinear : public Pair {
 public:
  PairMorseSmoothLinear(class LAMMPS *);
  virtual ~PairMorseSmoothLinear();

  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  void write_data(FILE *, char **);
  void write_data_all(FILE *, char **);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(const char *, int &);

 protected:
  double cut_global;
  double **cut;
  double **d0,**alpha,**r0;
  double **morse1;
  double **der_at_cutoff;
  double **offset;

  void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Incorrect args for pair coefficients

Self-explanatory.  Check the input script or data file.

E: All pair coeffs are not set

All pair coefficients must be set in the data file or by the
pair_coeff command before running a simulation.

*/
