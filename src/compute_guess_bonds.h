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

#ifdef COMPUTE_CLASS
// clang-format off
ComputeStyle(guess_bonds,ComputeGuessBonds);
// clang-format on
#else

#ifndef LMP_COMPUTE_GUESS_BONDS_H
#define LMP_COMPUTE_GUESS_BONDS_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeGuessBonds : public Compute {
 public:
  ComputeGuessBonds(class LAMMPS *, int, char **);
  ~ComputeGuessBonds() override;
  void init() override;
  void init_list(int, class NeighList *) override;
  void compute_peratom() override; // I guess this compute should be guess/bonds/atom
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

 private:
  class NeighList *list;
  double prefactor;
  std::vector<double> radii;
  std::vector<std::vector<double>> cutsq;
  int maxlocal;        // size of atom selection and variable arrays
  int ncol;            // number of columns is atom->bond_per_atom + 1
  double **carray;     // first column is num_bonds, rest list tags of bonded atoms
  int nchoose;         // # of selected atoms
  double *buf;         // memory for atom quantities
  int *choose;         // local indices of selected atoms
  int *clist;          // compressed list of indices of selected atoms
  int *chooseghost;    // extended choose array for comm
  double **bufcopy;    // buffer for communicating bond/atom info
  int maxbufcopy;
};

}    // namespace LAMMPS_NS

#endif
#endif
