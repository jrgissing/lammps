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
#include "comm.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "label_map.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "respa.h"
#include "update.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGuessBonds::FixGuessBonds(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), choose(nullptr), clist(nullptr), chooseghost(nullptr), bufcopy(nullptr)
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

  maxbufcopy = 0;
  maxlocal = 0;
}

/* ---------------------------------------------------------------------- */

int FixGuessBonds::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::init()
{
  if (utils::strmatch(update->integrate_style, "^respa"))
    nlevels_respa = (dynamic_cast<Respa *>(update->integrate))->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::setup(int vflag)
{
  if (utils::strmatch(update->integrate_style, "^verlet"))
    post_force(vflag);
  else {
    (dynamic_cast<Respa *>(update->integrate))->copy_flevel_f(nlevels_respa - 1);
    post_force_respa(vflag, nlevels_respa - 1, 0);
    (dynamic_cast<Respa *>(update->integrate))->copy_f_flevel(nlevels_respa - 1);
  }
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::post_force(int /*vflag*/)
{
  int i,j,m,atom1,atom2;

  auto *list = neighbor->get_best_pair_list();
  if (!list)
    if (comm->me == 0)
      error->warning(FLERR, "No suitable existing neighbor list for fix guess_bonds found");
  double **x = atom->x;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *type = atom->type;
  const int nlocal = atom->nlocal;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  // grow choose array if needed

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;

    //if (choose)
    memory->destroy(choose);
    //if (clist)
    memory->destroy(clist);
    memory->create(choose,maxlocal,"dump:choose");
    memory->create(clist,maxlocal,"dump:clist");
  }

  // remove all bonds
  int *num_bond = atom->num_bond;
  for (i = 0; i < nlocal; i++)
    num_bond[i] = 0;

  // choose all local atoms for output

  for (i = 0; i < nlocal; i++) choose[i] = 1;

  // un-choose if not in group

  if (igroup) {
    int *mask = atom->mask;
    for (i = 0; i < nlocal; i++)
      if (!(mask[i] & groupbit))
        choose[i] = 0;
  }

  // compress choose flags into clist
  // nchoose = # of selected atoms
  // clist[i] = local index of each selected atom

  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  // grab some suitable neighbor list, if available

  comm_forward = 1;
    //int nlocal = atom->nlocal;

    // communicate choose flag for ghost atoms to know if they are selected
    // if bcolor/bdiam = ATOM, setup bufcopy to comm atom color/diam attributes

  if (atom->nmax > maxbufcopy) {
    maxbufcopy = atom->nmax;
    memory->destroy(chooseghost);
    memory->create(chooseghost,maxbufcopy,"dump:chooseghost");
    if (comm_forward == 3) {
      memory->destroy(bufcopy);
      memory->create(bufcopy,maxbufcopy,2,"dump:bufcopy");
    }
  }
  for (i = 0; i < nlocal; i++) chooseghost[i] = choose[i];

  //if (comm_forward == 3) {
  //  for (i = 0; i < nlocal; i++) bufcopy[i][0] = bufcopy[i][1] = 0.0;
  //  m = 0;
  //  for (i = 0; i < nchoose; i++) {
  //    j = clist[i];
  //    bufcopy[j][0] = buf[m];
  //    bufcopy[j][1] = buf[m+1];
  //    m += size_one;
  //  }
  //}

  comm->forward_comm(this);
  //bool checkmass = false;

  // check for element by mass only for units "real" or "metal"
  //if ((strcmp(update->unit_style, "real") == 0) || (strcmp(update->unit_style, "metal") == 0))
  //  checkmass = true;

  for (int ii = 0; ii < nchoose; ii++) {
    atom1 = clist[ii];
    const double xtmp = x[atom1][0];
    const double ytmp = x[atom1][1];
    const double ztmp = x[atom1][2];

    // loop over neighbors
    auto *jlist = firstneigh[atom1];
    const int jnum = numneigh[atom1];
    for (int jj = 0; jj < jnum; ++jj) {
      atom2 = jlist[jj] & NEIGHMASK;
      if (!chooseghost[atom2]) continue;
      if ((newton_pair == 0) && (tag[atom1] > tag[atom2])) continue;
      // skip hydrogen-hydrogen bonds for units real or metal based on their mass
      // this is primarily for water and alkyl groups
      //if (checkmass) {
      //  if (rmass) {
      //    if ((rmass[atom1] < 3.0) && (rmass[atom2] < 3.0)) continue;
      //  } else {
      //    if ((mass[type[atom1]] < 3.0) && (mass[type[atom2]] < 3.0)) continue;
      //  }
      //}
      double dx = x[atom2][0] - xtmp;
      double dy = x[atom2][1] - ytmp;
      double dz = x[atom2][2] - ztmp;
      const double rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < cutsq[type[atom1]-1][type[atom2]-1]) {
        if (num_bond[atom1] == atom->bond_per_atom)
          error->one(FLERR, "New bond exceeded bonds per atom limit of {} in fix guess_bonds",
                     atom->bond_per_atom);
        bond_type[atom1][num_bond[atom1]] = 1; // all bonds set to type 1
        bond_atom[atom1][num_bond[atom1]] = tag[atom2];
        num_bond[atom1]++;
      }
    }
  }

  // recount bonds

  bigint nbonds = 0;
  for (i = 0; i < nlocal; i++) nbonds += num_bond[i];

  MPI_Allreduce(&nbonds, &atom->nbonds, 1, MPI_LMP_BIGINT, MPI_SUM, world);
  if (!force->newton_bond) atom->nbonds /= 2;
  //force reneighbor?
}

/* ---------------------------------------------------------------------- */

int FixGuessBonds::pack_forward_comm(int n, int *list, double *buf,
                                 int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;

  if (comm_forward == 1) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = chooseghost[j];
    }
  } else {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = chooseghost[j];
      buf[m++] = bufcopy[j][0];
      buf[m++] = bufcopy[j][1];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  if (comm_forward == 1)
    for (i = first; i < last; i++) chooseghost[i] = static_cast<int>(buf[m++]);
  else {
    for (i = first; i < last; i++) {
      chooseghost[i] = static_cast<int>(buf[m++]);
      bufcopy[i][0] = buf[m++];
      bufcopy[i][1] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::post_force_respa(int vflag, int ilevel, int /*iloop*/)
{
  if (ilevel == nlevels_respa - 1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixGuessBonds::min_post_force(int vflag)
{
  post_force(vflag);
}
