/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(guess_bonds,FixGuessBonds);
// clang-format on
#else

#ifndef LMP_FIX_GUESS_BONDS_H
#define LMP_FIX_GUESS_BONDS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGuessBonds : public Fix {
 public:
  FixGuessBonds(class LAMMPS *, int, char **);
  //~FixGuessBonds() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void min_setup(int) override;
  void post_force(int) override;
  void post_force_respa(int, int, int) override;
  void min_post_force(int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  double prefactor;
  std::vector<double> radii;
  std::vector<std::vector<double>> cutsq;
  int maxlocal;        // size of atom selection and variable arrays
  int nchoose;         // # of selected atoms
  double *buf;         // memory for atom quantities
  int *choose;         // local indices of selected atoms
  int *clist;          // compressed list of indices of selected atoms
  int *chooseghost;    // extended choose array for comm
  double **bufcopy;    // buffer for communicating bond/atom info
  int maxbufcopy;
  int nlevels_respa;
};

}    // namespace LAMMPS_NS

#endif
#endif
