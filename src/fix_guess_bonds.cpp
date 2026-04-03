// clang-format off
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

#include "fix_guess_bonds.h"

#include "atom.h"
#include "compute_guess_bonds.h"
#include "error.h"
#include "force.h"
#include "label_map.h"
#include "modify.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGuessBonds::FixGuessBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), cgb(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR,"fix guess_bonds", error);
  dynamic_group_allow = 1;
  if (strcmp(arg[3],"radii") != 0) error->all(FLERR,"Unknown fix guess_bonds keyword {}", arg[3]);
  prefactor = utils::numeric(FLERR, arg[4], false, lmp);
  int ntypes = atom->ntypes;
  radii.resize(ntypes);
  cutsq.resize(ntypes);
  for (auto &row : cutsq) row.resize(ntypes);

  int iarg = 5;
  int mytype;
  std::string radii_list;
  for (int i = 4; i < narg; i++)
    radii_list += fmt::format("{} ", arg[i]);

  while (iarg < narg) {
    std::string typestr = utils::utf8_subst(arg[iarg]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        mytype = utils::inumeric(FLERR, typestr, false, lmp);
      }
      case 1: {    // type label
        if (!atom->labelmapflag)
          error->all(FLERR, "Invalid atom type {} in fix guess_bonds", typestr);
        mytype = atom->lmap->find_type(typestr, Atom::ATOM);
        if (mytype == -1)
          error->all(FLERR, "Unknown atom type {} in fix guess_bonds", typestr);
        break;
      }
      default:    // invalid
        error->all(FLERR, "Invalid format in fix guess_bonds");
        break;
    }
    radii[mytype-1] = utils::numeric(FLERR, arg[iarg+1], false, lmp);
    iarg += 2;
  }

  for (auto radius : radii)
    if (radius <= 0.0)
      error->all(FLERR, "Fix guess_bonds: A positive radius must be provided for every atom type");

  for (int i = 0; i < atom->ntypes; i++) {
    for (int j = 0; j < atom->ntypes; j++) {
      cutsq[i][j] = prefactor*(radii[i]+radii[j]);
      cutsq[i][j] *= cutsq[i][j];
    }
  }

  // create instances of compute guess_bonds

  char *groupid = arg[1];
  std::string computeid = fmt::format("{}_compute_guess_bonds", id);
  std::string check = fmt::format("{} {} guess_bonds radii {}", computeid, groupid, radii_list);
  cgb = dynamic_cast<ComputeGuessBonds *>(modify->add_compute(
        fmt::format("{} {} guess_bonds radii {}", computeid, groupid, radii_list)));
}

/* ---------------------------------------------------------------------- */

int FixGuessBonds::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

void FixGuessBonds::end_of_step()
{
  int **bond_type = atom->bond_type;
  int *num_bond = atom->num_bond;
  tagint **bond_atom = atom->bond_atom;
  cgb->compute_peratom();
  // newton bond considerations?
  for (int i = 0; i < atom->nlocal; i++) {
    num_bond[i] = cgb->array_atom[i][0];
    for (int j = 0; j < num_bond[i]; j++) {
      bond_type[i][j] = 1; // all bonds set to type 1
      bond_atom[i][j] = cgb->array_atom[i][j+1];
    }
  }

  // recount bonds

  bigint nbonds = 0;
  for (int i = 0; i < atom->nlocal; i++) nbonds += num_bond[i];

  MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (!force->newton_bond) atom->nbonds /= 2;
  //force reneighbor??
}
