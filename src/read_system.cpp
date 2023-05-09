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

#include "read_system.h"

#include "comm.h"
#include "read_data.h"
#include "error.h"

using namespace LAMMPS_NS;

// clang-format on
/* ---------------------------------------------------------------------- */
ReadSystem::ReadSystem(LAMMPS *_lmp) : Command(_lmp)
{
  MPI_Comm_rank(world, &me);
}

/* ---------------------------------------------------------------------- */

ReadSystem::~ReadSystem()
{
}

// clang format off

/* ---------------------------------------------------------------------- */

void ReadSystem::command(int narg, char **arg)
{
  if (narg < 1) utils::missing_cmd_args(FLERR, "read_system", error);

  MPI_Barrier(world);
  double time1 = platform::walltime();

  if (!(strcmp(arg[0],"read_data") == 0 || strcmp(arg[0],"read_molecule") == 0))
    error->all(FLERR,"Invalid syntax for read_system command");

  // let's find number of files to be read
  std::vector<int> vnarg;
  //vnarg.push_back(0);
  int nfiles = 0;
  int argc = 1;
  for (int i = 0; i < narg; i++) {
    if (strcmp(arg[i],"read_data") == 0 || strcmp(arg[i],"read_molecule") == 0) {
      //vnarg.push_back(argc);
      nfiles++;
      vnarg.push_back(i);
      argc = 1;
    }
    argc++;
  }
  vnarg.push_back(narg);

  printf("nfiles %d %d\n",nfiles,narg);
  for (int i = 0; i < nfiles; i++)
  {

    read_data = new ReadData(lmp);
    read_data->internal_extractflag = 1;
    read_data->command(vnarg[i+1]-vnarg[i],&arg[vnarg[i]+1]);
    printf("why %d\n",vnarg[i]);
    printf("pastfirst %d %s\n",vnarg[i+1]-vnarg[i],arg[vnarg[i]+1]);

    // full read data routine after determing max topology
    read_data->internal_extractflag = 0;
    //read_data->addflag = ReadData::MERGE;
    //read_data->command(narg,arg);
    read_data->command(vnarg[i+1]-vnarg[i],&arg[vnarg[i]+1]);
    delete [] read_data;
  }

  // total time

  //atom->readsysflag = 1;
  //int index = 1;
  //molecule = new Molecule(lmp,narg,arg,index);

  MPI_Barrier(world);

  if (comm->me == 0)
    utils::logmesg(lmp, "  read_system CPU = {:.3f} seconds\n", platform::walltime() - time1);

  //delete [] read_data;
  //delete [] molecule;
}
