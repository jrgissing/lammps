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
  void end_of_step() override;

 private:
  double prefactor;
  std::vector<double> radii;
  std::vector<std::vector<double>> cutsq;
  class ComputeGuessBonds *cgb;
};

}    // namespace LAMMPS_NS

#endif
#endif
