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

/* ----------------------------------------------------------------------
   Contributing Author: Jacob Gissinger (jgissing@stevens.edu)
------------------------------------------------------------------------- */

#ifndef LMP_REACTION_PARSER_H
#define LMP_REACTION_PARSER_H

#include "pointers.h"    // IWYU pragma: export
#include "reaction.h"

#include <array>
#include <string>
#include <vector>

namespace LAMMPS_NS {

class ReactionParser : protected Pointers {
public:
  ReactionParser(LAMMPS *lmp) : Pointers(lmp) {}

  void parse_reaction(char **, int &, int, Reaction &);

  void validate_variable_keyword(const char *, int);
};

}    // namespace LAMMPS_NS

#endif
