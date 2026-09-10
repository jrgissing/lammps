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
FixStyle(assign_bonds,FixAssignBonds);
// clang-format on
#else

#ifndef LMP_FIX_ASSIGN_BONDS_H
#define LMP_FIX_ASSIGN_BONDS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAssignBonds : public Fix {
 public:
  FixAssignBonds(class LAMMPS *, int, char **);
  //~FixAssignBonds() override;
  void post_constructor() override;
  int setmask() override;
  void end_of_step() override;

 private:
  std::string groupid;
  std::string guess_bonds_id;
  int nevery_history;      // store a history frame once every Nevery steps
  int nrepeat_history;     // # of history frames to store
  int nfreq_history;       // enable output of stored history on these steps
  int bonded_fraction;     // cutoff for fraction of sampled timesteps that bond exists
  class FixStoreState *fss;
};

}    // namespace LAMMPS_NS

#endif
#endif
