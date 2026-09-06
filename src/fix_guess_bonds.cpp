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
#include "fix_store_state.h"
#include "force.h"
#include "label_map.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <cstdio>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGuessBonds::FixGuessBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), cgb(nullptr), fss(nullptr)
{
  if (narg < 11) utils::missing_cmd_args(FLERR,"fix guess_bonds", error);
  dynamic_group_allow = 1;
  int ntypes = atom->ntypes;
  radii.resize(ntypes);
  cutsq.resize(ntypes);
  for (auto &row : cutsq) row.resize(ntypes);

  nevery_history = utils::inumeric(FLERR, arg[3], false, lmp);
  nrepeat_history = utils::inumeric(FLERR, arg[4], false, lmp);
  nfreq_history = utils::inumeric(FLERR, arg[5], false, lmp);
  bond_order_cutoff = utils::numeric(FLERR, arg[6], false, lmp);
  //bond_order_cutoff = 0.3; // default bond order cutoff

  if (strcmp(arg[7],"radii") != 0) error->all(FLERR,"Unknown fix guess_bonds keyword {}", arg[6]);

  prefactor = utils::numeric(FLERR, arg[8], false, lmp);

  int iarg = 8;
  int mytype;
  //std::string radii_list;
  for (int i = iarg; i < narg; i++)
    radii_list += fmt::format("{} ", arg[i]);

  iarg++;

  while (iarg < narg) {
    std::string typestr = utils::utf8_subst(arg[iarg]);
    switch (utils::is_type(typestr)) {
      case 0: {    // numeric
        mytype = utils::inumeric(FLERR, typestr, false, lmp);
        break;
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
        error->all(FLERR, "Invalid keyword {} in fix guess_bonds", typestr);
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

  groupid = arg[1];
  //std::string computeid = fmt::format("{}_compute_guess_bonds", id);
  //std::string check = fmt::format("{} {} guess_bonds radii {}", computeid, groupid, radii_list);
  //cgb = dynamic_cast<ComputeGuessBonds *>(modify->add_compute(
  //      fmt::format("{} {} guess_bonds radii {}", computeid, groupid, radii_list)));

  //return;

  std::string fss_fixid = fmt::format("{}_fix_store_state", id);

  //fss = dynamic_cast<FixStoreState *>(modify->get_fix_by_id(fss_fixid));
  //if (!fss)
  //  fss = dynamic_cast<FixStoreState *>(modify->add_fix(
  //        fmt::format("{} {} store/state 0 c_{}[*] history {} {} {}", fss_fixid, groupid, computeid,
  //                    nevery_history, nrepeat_history, nfreq_history)));

}

void FixGuessBonds::post_constructor()
{
  // create instances of compute guess_bonds

  //char *groupid = arg[1];
  std::string computeid = fmt::format("{}_compute_guess_bonds", id);
  std::string check = fmt::format("{} {} guess_bonds radii {}", computeid, groupid, radii_list);
  cgb = dynamic_cast<ComputeGuessBonds *>(modify->add_compute(
        fmt::format("{} {} guess_bonds radii {}", computeid, groupid, radii_list)));

  std::string fss_fixid = fmt::format("{}_fix_store_state", id);

  fss = dynamic_cast<FixStoreState *>(modify->get_fix_by_id(fss_fixid));
  if (!fss)
    fss = dynamic_cast<FixStoreState *>(modify->add_fix(
          fmt::format("{} {} store/state 0 c_{}[*] history {} {} {}", fss_fixid, groupid, computeid,
                      nevery_history, nrepeat_history, nfreq_history)));
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
  //cgb->compute_peratom();

  int dim;
  int flag_history = *((int *) fss->extract("flag_history",dim));
  if (!flag_history)
    error->all(FLERR,"Fix {} store/state fix does not store history", style);

  int *count_history_ptr = (int *) fss->extract("count_history",dim);
  int n_valid_history = *count_history_ptr;
  if (n_valid_history > nrepeat_history) n_valid_history = nrepeat_history;

  int *most_recent_index_ptr = (int *) fss->extract("most_recent_index",dim);
  int most_recent_index = *most_recent_index_ptr;

  double ***history = (double ***) fss->extract("history",dim);

  int size_peratom_cols = atom->bond_per_atom + 1; // should extract from compute_guess_bonds instead

  std::vector<std::vector<int>> ave_bond_atoms(atom->nlocal, std::vector<int>(size_peratom_cols, 0));
  std::vector<std::vector<double>> ave_bond_persistence(atom->nlocal, std::vector<double>(size_peratom_cols, 0));

  if (n_valid_history > 0) {
  int hslot = most_recent_index;
  for (int n = 0; n < n_valid_history; n++) {
    int i = hslot;
    for (int j = 0; j < atom->nlocal; j++) {
      int num_bond = history[i][j][0];
      for (int k = 0; k < num_bond; k++) {
        int num_ave_bonds = ave_bond_atoms[j][0];
        int found_bond = 0;
        for (int kk = 0; kk < num_ave_bonds; kk++) {
          if (history[i][j][k+1] == ave_bond_atoms[j][kk+1]) {
            ave_bond_persistence[j][kk+1]++;
            found_bond = 1;
            break;
          }
        }
        if (!found_bond) {
          if (ave_bond_atoms[j][0] >= size_peratom_cols-1)
            error->one(FLERR,"Fix guess_bonds: too many bonds per atom, increase bonds per atom");


          ave_bond_atoms[j][0]++;
          ave_bond_atoms[j][ave_bond_atoms[j][0]] = history[i][j][k+1];
          ave_bond_persistence[j][ave_bond_atoms[j][0]] = 1;
        }
      }
    }
    hslot--;
    if (hslot < 0) hslot += nrepeat_history;
  }

  for (int i = 0; i < atom->nlocal; i++) {
    num_bond[i] = 0;
    int num_ave_bonds = ave_bond_atoms[i][0];
    for (int j = 0; j < num_ave_bonds; j++) {
      ave_bond_persistence[i][j+1] /= n_valid_history;

      if (ave_bond_persistence[i][j+1] > 1.0 + 1e-9) {
        error->warning(FLERR, "guess_bonds: persistence {} > 1 for atom {} -- "
                        "unexpected, please report", ave_bond_persistence[i][j+1],
                        atom->tag[i]);
        //printf("what %d %d\n",n_valid_history,num_ave_bonds);
        //for (int m = 0; m < num_ave_bonds; m++) {
        //  printf("that %d\n",ave_bond_atoms[i][m+1]);
        //}
      }

      if (ave_bond_persistence[i][j+1] > bond_order_cutoff) {
        tagint tag_j = ave_bond_atoms[i][j+1];

        if (force->newton_bond && atom->tag[i] > tag_j) continue;

        bond_type[i][num_bond[i]] = 1; // all bonds set to type 1
        bond_atom[i][num_bond[i]] = tag_j;
        num_bond[i]++;
      }
    }
  }
  } else {
    for (int i = 0; i < atom->nlocal; i++) num_bond[i] = 0;
  } // end if (n_valid_history > 0)



  // use Pauling bond order definition?
  // option to weight more recent timesteps higher?
  // newton bond considerations?
  //for (int i = 0; i < atom->nlocal; i++) {
  //  num_bond[i] = cgb->array_atom[i][0];
  //  for (int j = 0; j < num_bond[i]; j++) {
  //    bond_type[i][j] = 1; // all bonds set to type 1
  //    bond_atom[i][j] = cgb->array_atom[i][j+1];
  //  }
  //}

  // recount bonds

  bigint nbonds = 0;
  for (int i = 0; i < atom->nlocal; i++) nbonds += num_bond[i];

  MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (!force->newton_bond) atom->nbonds /= 2;
  //force reneighbor??
}
