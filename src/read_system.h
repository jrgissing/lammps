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

#ifdef COMMAND_CLASS
// clang-format off
CommandStyle(read_system,ReadSystem);
// clang-format on
#else

#ifndef LMP_READ_SYSTEM_H
#define LMP_READ_SYSTEM_H

#include "command.h"

namespace LAMMPS_NS {
class Fix;
class ReadSystem : public Command {
 public:
  ReadSystem(class LAMMPS *);
  ~ReadSystem() override;
  void command(int, char **) override;

 private:
  int me, compressed;
  class ReadData *read_data;
  class Molecule *molecule;
};

}    // namespace LAMMPS_NS

#endif
#endif
