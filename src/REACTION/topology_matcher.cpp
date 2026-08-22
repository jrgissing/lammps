/* ----------------------------------------------------------------------
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

#include "topology_matcher.h"

#include "atom.h"
#include "molecule.h"

using namespace LAMMPS_NS;

TopologyMatcher::TopologyMatcher(LAMMPS *lmp) : Pointers(lmp) {
  status = Status::PROCEED;
}

/* ----------------------------------------------------------------------
  Check that newly assigned atoms have correct bonds
  Necessary for certain ringed structures
------------------------------------------------------------------------- */

int TopologyMatcher::ring_check(Reaction &rxn, std::vector<tagint> &glove)
{
  // full special lists - may need correction for unusual special bond settings
  int **nxspecial = atom->nspecial;
  tagint **xspecial = atom->special;

  // ring_check can be made more efficient by re-introducing 'frozen' atoms
  // 'frozen' atoms have been assigned and also are no longer pioneers

  // double check the number of neighbors match for all non-edge atoms
  // otherwise, atoms at 'end' of symmetric ring can behave like edge atoms
  for (int i = 0; i < rxn.reactant->natoms; i++)
    if (rxn.atoms[i].edge == 0 &&
        rxn.reactant->nspecial[i][0] != nxspecial[atom->map(glove[i])][0])
      return 0;

  for (int i = 0; i < rxn.reactant->natoms; i++) {
    for (int j = 0; j < rxn.reactant->nspecial[i][0]; j++) {
      int ring_fail = 1;
      int ispecial = rxn.reactant->special[i][j];
      for (int k = 0; k < nxspecial[atom->map(glove[i])][0]; k++) {
        if (xspecial[atom->map(glove[i])][k] == glove[ispecial-1]) {
          ring_fail = 0;
          break;
        }
      }
      if (ring_fail == 1) return 0;
    }
  }
  return 1;
}
