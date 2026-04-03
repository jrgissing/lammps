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

#include "compute_guess_bonds.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "label_map.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"

#include <cstring>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeGuessBonds::ComputeGuessBonds(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), choose(nullptr), clist(nullptr), chooseghost(nullptr), bufcopy(nullptr),  carray(nullptr)
{
  if (narg < 7) utils::missing_cmd_args(FLERR,"compute guess_bonds", error);
  dynamic_group_allow = 1;
  ncol = atom->bond_per_atom + 1;
  if (strcmp(arg[3],"radii") != 0) error->all(FLERR,"Unknown compute guess_bonds keyword {}", arg[3]);
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
          error->all(FLERR, "Invalid atom type {} in compute guess_bonds", typestr);
        mytype = atom->lmap->find_type(typestr, Atom::ATOM);
        if (mytype == -1)
          error->all(FLERR, "Unknown atom type {} in compute guess_bonds", typestr);
        break;
      }
      default:    // invalid
        error->all(FLERR, "Invalid format in compute guess_bonds");
        break;
    }
    radii[mytype-1] = utils::numeric(FLERR, arg[iarg+1], false, lmp);
    iarg += 2;
  }

  for (auto radius : radii)
    if (radius <= 0.0)
      error->all(FLERR, "Compute guess_bonds: A positive radius must be provided for every atom type");

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

ComputeGuessBonds::~ComputeGuessBonds()
{
  memory->destroy(carray);
  //need to destroy more stuff
}

/* ---------------------------------------------------------------------- */

void ComputeGuessBonds::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeGuessBonds::init()
{
  auto *req = neighbor->add_request(this, NeighConst::REQ_OCCASIONAL);
}

/* ---------------------------------------------------------------------- */

void ComputeGuessBonds::compute_peratom()
{
  int i,j,m,atom1,atom2;

  neighbor->build_one(list);
  if (!list)
    if (comm->me == 0)
      error->warning(FLERR, "No suitable existing neighbor list for compute guess_bonds found");
  double **x = atom->x;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *type = atom->type;
  const int nlocal = atom->nlocal;

  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  // grow choose array if needed

  if (atom->nmax > maxlocal) {
    maxlocal = atom->nmax;
    memory->destroy(carray);
    memory->create(carray, maxlocal, ncol, "guess_bonds:carray");
    array_atom = carray;

    //if (choose)
    memory->destroy(choose);
    //if (clist)
    memory->destroy(clist);
    memory->create(choose,maxlocal,"dump:choose");
    memory->create(clist,maxlocal,"dump:clist");
  }

  // remove all bonds

  for (i = 0; i < nlocal; i++)
    carray[i][0] = 0;

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
  nchoose = 0;
  for (i = 0; i < nlocal; i++)
    if (choose[i]) clist[nchoose++] = i;

  // communicate choose flag for ghost atoms to know if they are selected

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

  comm->forward_comm(this);

  // mostly borrowed from fix graphics autobond
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

      double dx = x[atom2][0] - xtmp;
      double dy = x[atom2][1] - ytmp;
      double dz = x[atom2][2] - ztmp;
      //domain->minimum_image(FLERR, dx, dy, dz);
      const double rsq = dx*dx + dy*dy + dz*dz;

      if (rsq < cutsq[type[atom1]-1][type[atom2]-1]) {
        if (carray[atom1][0] == atom->bond_per_atom)
          error->one(FLERR, "New bond exceeded bonds per atom limit of {} in compute guess_bonds",
                     atom->bond_per_atom);
        carray[atom1][(int) ++carray[atom1][0]] = tag[atom2];
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int ComputeGuessBonds::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/,
                                        int * /*pbc*/)
{
  int i,j,m;

  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = chooseghost[j];
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void ComputeGuessBonds::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;

  for (i = first; i < last; i++) chooseghost[i] = static_cast<int>(buf[m++]);
}
