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

#include "reaction_parser.h"

#include "atom.h"
#include "error.h"
#include "fix_bond_react.h"
#include "group.h"
#include "input.h"
#include "molecule.h"
#include "variable.h"

using namespace LAMMPS_NS;

void ReactionParser::parse_reaction(char **arg, int iarg, int rxn_narg, Reaction &rxn) {

  rxn.name = arg[iarg++];
  if (rxn.name.size()+1 > MAXNAME) error->all(FLERR,"Reaction name (react-ID) is too long (limit: 255 characters)");

  int groupid = group->find(arg[iarg++]);
  if (groupid == -1) error->all(FLERR,"Could not find fix group ID");
  rxn.groupbits = group->bitmask[groupid];

  if (strncmp(arg[iarg],"v_",2) == 0) {
    rxn.v_nevery = input->variable->find(&arg[iarg][2]);
    validate_variable_keyword(&arg[iarg][2], rxn.v_nevery);
  } else {
    rxn.nevery = utils::inumeric(FLERR,arg[iarg],false,lmp);
    if (rxn.nevery <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                     "'Nevery' must be a positive integer");
  }
  iarg++;

  double cutoff;
  if (strncmp(arg[iarg],"v_",2) == 0) {
    rxn.v_rmin = input->variable->find(&arg[iarg][2]);
    validate_variable_keyword(&arg[iarg][2], rxn.v_rmin);
    cutoff = input->variable->compute_equal(rxn.v_rmin);
  } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
    if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command: "
                                 "'Rmin' cannot be negative");
    rxn.rminsq = cutoff*cutoff;
  iarg++;

  if (strncmp(arg[iarg],"v_",2) == 0) {
    rxn.v_rmax = input->variable->find(&arg[iarg][2]);
    validate_variable_keyword(&arg[iarg][2], rxn.v_rmax);
    cutoff = input->variable->compute_equal(rxn.v_rmax);
  } else cutoff = utils::numeric(FLERR,arg[iarg],false,lmp);
    if (cutoff < 0.0) error->all(FLERR,"Illegal fix bond/react command:"
                                 "'Rmax' cannot be negative");
    rxn.rmaxsq = cutoff*cutoff;
  iarg++;

  int mol_idx = atom->find_molecule(arg[iarg++]);
  if (mol_idx == -1) error->all(FLERR,"Pre-reaction molecule template ID for "
                                           "fix bond/react does not exist");
  rxn.reactant = atom->molecules[mol_idx];
  mol_idx = atom->find_molecule(arg[iarg++]);
  if (mol_idx == -1) error->all(FLERR,"Post-reaction molecule template ID for "
                                         "fix bond/react does not exist");
  rxn.product = atom->molecules[mol_idx];

  //read map file
  rxn.mapfilename = arg[iarg];
  iarg++;

  while (iarg < rxn_narg) {
    if (strcmp(arg[iarg],"prob") == 0) {
      if (iarg+3 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'prob' keyword has too few arguments");
      // check if probability is a variable
      if (strncmp(arg[iarg+1],"v_",2) == 0) {
        rxn.v_prob = input->variable->find(&arg[iarg+1][2]);
        validate_variable_keyword(&arg[iarg+1][2], rxn.v_prob);
        rxn.fraction = input->variable->compute_equal(rxn.v_prob);
      } else {
        // otherwise probability should be a number
        rxn.fraction = utils::numeric(FLERR,arg[iarg+1],false,lmp);
      }
      rxn.seed = utils::inumeric(FLERR,arg[iarg+2],false,lmp);
      if (rxn.fraction < 0.0 || rxn.fraction > 1.0)
        error->all(FLERR,"Illegal fix bond/react command: "
                   "probability fraction must between 0 and 1, inclusive");
      if (rxn.seed <= 0) error->all(FLERR,"Illegal fix bond/react command: "
                                     "probability seed must be positive");
      iarg += 3;
    } else if (strcmp(arg[iarg],"stabilize_steps") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'stabilize_steps' has too few arguments");
      rxn.limit_duration = utils::inumeric(FLERR,arg[iarg+1],false,lmp);
      rxn.stabilize_steps_flag = 1;
      iarg += 2;
    } else if (strcmp(arg[iarg],"custom_charges") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'custom_charges' has too few arguments");
      if (strcmp(arg[iarg+1],"no") == 0) rxn.custom_charges_fragid = -1; //default
      else {
        rxn.custom_charges_fragid = rxn.reactant->findfragment(arg[iarg+1]);
        if (rxn.custom_charges_fragid < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                       "'custom_charges' keyword does not exist");
      }
      iarg += 2;
    } else if (strcmp(arg[iarg],"rescale_charges") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'rescale_charges' has too few arguments");
      if (strcmp(arg[iarg+1],"no") == 0) rxn.rescale_charges_flag = 0; //default
      else if (strcmp(arg[iarg+1],"yes") == 0) {
        if (!atom->q_flag) error->all(FLERR,"Illegal fix bond/react command: cannot use "
                                    "'rescale_charges' without atomic charges enabled");
        if (!rxn.product->qflag) error->all(FLERR,"Illegal fix bond/react command: cannot use "
                                    "'rescale_charges' without Charges section in post-reaction template");
        rxn.rescale_charges_flag = 1; // overloaded below to also indicate number of atoms to update
      } else error->one(FLERR,"Bond/react: Illegal option for 'rescale_charges' keyword");
      iarg += 2;
    } else if (strcmp(arg[iarg],"molecule") == 0) {
      if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'molecule' has too few arguments");
      if (strcmp(arg[iarg+1],"off") == 0) rxn.molecule_keyword = Molecule_Keys::OFF; //default
      else if (strcmp(arg[iarg+1],"inter") == 0) rxn.molecule_keyword = Molecule_Keys::INTER;
      else if (strcmp(arg[iarg+1],"intra") == 0) rxn.molecule_keyword = Molecule_Keys::INTRA;
      else error->one(FLERR,"Fix bond/react: Illegal option for 'molecule' keyword");
      iarg += 2;
    } else if (strcmp(arg[iarg],"modify_create") == 0) {
      if (iarg++ > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                    "'modify_create' has too few arguments");
      while (iarg < rxn_narg && strcmp(arg[iarg],"react") != 0) {
        if (strcmp(arg[iarg],"fit") == 0) {
          if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                        "'modify_create' has too few arguments");
          if (strcmp(arg[iarg+1],"all") == 0) rxn.modify_create_fragid = -1; //default
          else {
            rxn.modify_create_fragid = rxn.product->findfragment(arg[iarg+1]);
            if (rxn.modify_create_fragid < 0) error->one(FLERR,"Fix bond/react: Molecule fragment for "
                                                           "'modify_create' keyword does not exist");
          }
          iarg += 2;
        } else if (strcmp(arg[iarg],"overlap") == 0) {
          if (iarg+2 > rxn_narg) error->all(FLERR,"Illegal fix bond/react command: "
                                        "'modify_create' has too few arguments");
          rxn.overlapsq = utils::numeric(FLERR,arg[iarg+1],false,lmp);
          rxn.overlapsq *= rxn.overlapsq;
          iarg += 2;
        } else break;
      }
    } else if (strcmp(arg[iarg],"rate_limit") == 0) {
      error->all(FLERR,"Fix bond/react: 'rate_limit' as an 'individual keyword' has been deprecated. "
                       "Please use the 'rate_limit' common keyword instead, which can be applied to one or more reactions.");
    } else if (strcmp(arg[iarg],"max_rxn") == 0) {
      error->all(FLERR,"Fix bond/react: 'max_rxn' as an 'individual keyword' has been deprecated. "
                       "Please use the 'max_rxn' common keyword instead, which can be applied to one or more reactions.");
    } else error->all(FLERR,"Illegal fix bond/react command: unknown keyword");
  }
};

/* ----------------------------------------------------------------------
add equal-style variable to keyword argument list
------------------------------------------------------------------------- */

void ReactionParser::validate_variable_keyword(const char *myarg, int var_id)
{
  if (var_id < 0)
    error->all(FLERR,"Fix bond/react: Variable name {} does not exist",myarg);
  if (!input->variable->equalstyle(var_id))
    error->all(FLERR,"Fix bond/react: Variable {} is not equal-style",myarg);
}
