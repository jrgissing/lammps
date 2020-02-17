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

/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (jewett@caltech.edu)
   adapted from fix_restrain.cpp (by Craig Tenney and Andres Jaramillo-Botero),
   and fix_ave_atom.cpp
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(twist,FixTwist)

#else

#ifndef LMP_FIX_TWIST_H
#define LMP_FIX_TWIST_H

#include "fix.h"

namespace LAMMPS_NS {

class FixTwist : public Fix {
 public:
  FixTwist(class LAMMPS *, int, char **);
  ~FixTwist();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

  void setup_pre_force(int);
  void setup_min_pre_force(int);
  // The next two functions need to be defined so that LAMMPS does not ignore
  // the setup_pre_force() and setup_min_pre_force() functions.  (IE. so that
  // "n_pre_force" and "n_min_pre_force" are both > 0.  See modify.cpp)
  void pre_force(int);
  void min_pre_force(int);

  double compute_scalar();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);

  void calc_phi_prev(); //initialize the phi_prev[] array

 private:
  int nlevels_respa;
  tagint ntwist,maxtwist;
  tagint **ids;
  double *k;
  double *param_start,*param_stop;
  double energy,energy_all;

  // The next array stores the value of the "phi" angle between subequent
  // invocations of compute_force().  It is indexed by the "central" atom (2nd
  // atom) in the list of dihedral atoms.  As this atom migrates to different
  // processors, the entries in the phi_prev[] array must migrate as well.
  double *phi_prev;

  double Phi(double const *x1, //array holding x,y,z coords atom 1
	     double const *x2, // :       :      :      :        2
	     double const *x3, // :       :      :      :        3
	     double const *x4, // :       :      :      :        4
	     Domain *domain, //<-periodic boundary information
	     double *vb12,  // will store x2-x1
	     double *vb23,  // will store x3-x2
	     double *vb34,  // will store x4-x3
	     double *n123,  // will store normal to plane x1,x2,x3
	     double *n234); // will store normal to plane x2,x3,x4

  void compute_force(tagint);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Fix twist requires an atom map, see atom_modify

Self-explanatory.

E: Twist atoms %d %d missing on proc %d at step %ld

The 2 atoms in a twist bond specified by the fix twist
command are not all accessible to a processor.  This probably means an
atom has moved too far.

E: Twist atoms %d %d %d missing on proc %d at step %ld

The 3 atoms in a twist angle specified by the fix twist
command are not all accessible to a processor.  This probably means an
atom has moved too far.

E: Twist atoms %d %d %d %d missing on proc %d at step %ld

The 4 atoms in a twist dihedral specified by the fix twist
command are not all accessible to a processor.  This probably means an
atom has moved too far.

W: Twist problem: %d %ld %d %d %d %d

Conformation of the 4 listed dihedral atoms is extreme; you may want
to check your simulation geometry.

*/
