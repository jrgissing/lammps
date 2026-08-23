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
#include "error.h"
#include "molecule.h"

using namespace LAMMPS_NS;

TopologyMatcher::TopologyMatcher(LAMMPS *lmp) : Pointers(lmp) {

  rxn_constraints = new ReactionConstraints(lmp);

  status = Status::PROCEED;
}

TopologyMatcher::~TopologyMatcher()
{
  delete rxn_constraints;
}

/* ----------------------------------------------------------------------
  We are ready to make a guess. If there are multiple possible choices
  for this guess, keep track of these.
------------------------------------------------------------------------- */

void TopologyMatcher::inner_crosscheck_loop(Superimpose &super, Reaction &rxn)
{
  // full special lists - may need correction for unusual special bond settings
  int **nxspecial = atom->nspecial;
  tagint **xspecial = atom->special;

  Superimpose::StatePoint &sp = super.sp;
  int &avail_guesses = super.avail_guesses;
  std::vector<int> &guess_branch = super.guess_branch;

  int *type = atom->type;
  // arbitrarily limited to 5 identical first neighbors
  tagint tag_choices[5];
  int nfirst_neighs = rxn.reactant->nspecial[sp.pion][0];

  int num_choices = 0;
  for (int i = 0; i < nfirst_neighs; i++) {
    int checktype = type[(int)atom->map(xspecial[atom->map(sp.glove[sp.pion])][i])];
    int reactant_atom = (int) rxn.reactant->special[sp.pion][sp.neigh];
    if (compare_atomtype(checktype, rxn, reactant_atom)) {
      if (num_choices == 5) { // here failed because too many identical first neighbors. but really no limit if situation arises
        status = Status::GUESSFAIL;
        return;
      }
      tag_choices[num_choices++] = xspecial[atom->map(sp.glove[sp.pion])][i];
    }
  }

  // guess branch is for when multiple identical neighbors. then we guess each one in turn
  // guess_branch must work even when avail_guesses = 0. so index accordingly!
  // ...actually, avail_guesses should never be zero here anyway
  if (guess_branch[avail_guesses-1] == 0) guess_branch[avail_guesses-1] = num_choices;

  for (int i=1; i < num_choices; ++i) {
    tagint hold = tag_choices[i];
    int j = i - 1;
    while ((j >= 0) && (tag_choices[j] > hold)) {
      tag_choices[j+1] = tag_choices[j];
      --j;
    }
    tag_choices[j+1] = hold;
  }

  for (int i = guess_branch[avail_guesses-1]-1; i >= 0; i--) {
    int already_assigned = 0;
    for (int j = 0; j < rxn.reactant->natoms; j++) {
      if (sp.glove[j] == tag_choices[i]) {
        already_assigned = 1;
        break;
      }
    }
    if (already_assigned == 1) {
      guess_branch[avail_guesses-1]--;
      if (guess_branch[avail_guesses-1] == 0) {
        status = Status::REJECT;
        return;
      }
    } else {
      sp.glove[rxn.reactant->special[sp.pion][sp.neigh]-1] = tag_choices[i];
      guess_branch[avail_guesses-1]--;
      break;
    }
  }

  //another check for ghost atoms. perhaps remove the one in make_a_guess
  if (atom->map(sp.glove[(int)rxn.reactant->special[sp.pion][sp.neigh]-1]) < 0) {
    error->one(FLERR,"Fix bond/react: Fix bond/react needs ghost atoms from further away");
  }

  if (guess_branch[avail_guesses-1] == 0) avail_guesses--;

  for (int i = 0; i < rxn.reactant->nspecial[rxn.reactant->special[sp.pion][sp.neigh]-1][0]; i++) {
    sp.pioneer_count[rxn.reactant->special[rxn.reactant->special[sp.pion][sp.neigh]-1][i]-1]++;
  }
  sp.glove_counter++;
  if (sp.glove_counter == rxn.reactant->natoms) {
    if (ring_check(rxn, sp.glove) && rxn_constraints->check(rxn, sp.glove)) status = Status::ACCEPT;
    else status = Status::GUESSFAIL;
    return;
  }
  status = Status::CONTINUE;
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

/* ----------------------------------------------------------------------
check if an atom type matches that of a reactant atom
------------------------------------------------------------------------- */

bool TopologyMatcher::compare_atomtype(int checktype, Reaction &rxn, int reactant_atom)
{
  int iatom = reactant_atom - 1; // index of reactant atom
  if (checktype == rxn.reactant->type[iatom] || rxn.atoms[iatom].wildcard) return true;
  else return false;
}
