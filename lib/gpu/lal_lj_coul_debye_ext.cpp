/***************************************************************************
                             lj_coul_debye_ext.cpp
                             -------------------
                            Trung Dac Nguyen (ORNL)

  Functions for LAMMPS access to lj/cut/coul/debye acceleration routines.

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                : 
    email                : nguyentd@ornl.gov
 ***************************************************************************/

#include <iostream>
#include <cassert>
#include <math.h>

#include "lal_lj_coul_debye.h"

using namespace std;
using namespace LAMMPS_AL;

static LJCoulDebye<PRECISION,ACC_PRECISION> LJCDMF;

// ---------------------------------------------------------------------------
// Allocate memory on host and device and copy constants to device
// ---------------------------------------------------------------------------
int ljcd_gpu_init(const int ntypes, double **cutsq, double **host_lj1,
                  double **host_lj2, double **host_lj3, double **host_lj4,
                  double **offset, double *special_lj, const int inum,
                  const int nall, const int max_nbors, const int maxspecial,
                  const double cell_size, int &gpu_mode, FILE *screen,
                  double **host_cut_ljsq, double **host_cut_coulsq,
                  double *host_special_coul, const double qqrd2e, 
                  const double kappa) {
  LJCDMF.clear();
  gpu_mode=LJCDMF.device->gpu_mode();
  double gpu_split=LJCDMF.device->particle_split();
  int first_gpu=LJCDMF.device->first_device();
  int last_gpu=LJCDMF.device->last_device();
  int world_me=LJCDMF.device->world_me();
  int gpu_rank=LJCDMF.device->gpu_rank();
  int procs_per_gpu=LJCDMF.device->procs_per_gpu();

  LJCDMF.device->init_message(screen,"lj/cut/coul/debye",first_gpu,last_gpu);

  bool message=false;
  if (LJCDMF.device->replica_me()==0 && screen)
    message=true;

  if (message) {
    fprintf(screen,"Initializing Device and compiling on process 0...");
    fflush(screen);
  }

  int init_ok=0;
  if (world_me==0)
    init_ok=LJCDMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3,
                        host_lj4, offset, special_lj, inum, nall, 300,
                        maxspecial, cell_size, gpu_split, screen, host_cut_ljsq,
                        host_cut_coulsq, host_special_coul, qqrd2e, kappa);

  LJCDMF.device->world_barrier();
  if (message)
    fprintf(screen,"Done.\n");

  for (int i=0; i<procs_per_gpu; i++) {
    if (message) {
      if (last_gpu-first_gpu==0)
        fprintf(screen,"Initializing Device %d on core %d...",first_gpu,i);
      else
        fprintf(screen,"Initializing Devices %d-%d on core %d...",first_gpu,
                last_gpu,i);
      fflush(screen);
    }
    if (gpu_rank==i && world_me!=0)
      init_ok=LJCDMF.init(ntypes, cutsq, host_lj1, host_lj2, host_lj3, host_lj4,
                          offset, special_lj, inum, nall, 300, maxspecial,
                          cell_size, gpu_split, screen, host_cut_ljsq,
                          host_cut_coulsq, host_special_coul, qqrd2e, kappa);

    LJCDMF.device->gpu_barrier();
    if (message) 
      fprintf(screen,"Done.\n");
  }
  if (message)
    fprintf(screen,"\n");

  if (init_ok==0)
    LJCDMF.estimate_gpu_overhead();
  return init_ok;
}

void ljcd_gpu_clear() {
  LJCDMF.clear();
}

int** ljcd_gpu_compute_n(const int ago, const int inum_full,
                         const int nall, double **host_x, int *host_type,
                         double *sublo, double *subhi, int *tag, int **nspecial, 
                         int **special, const bool eflag, const bool vflag,
                         const bool eatom, const bool vatom, int &host_start,
                         int **ilist, int **jnum, const double cpu_time,
                         bool &success, double *host_q, double *boxlo,
                         double *prd) {
  return LJCDMF.compute(ago, inum_full, nall, host_x, host_type, sublo,
                        subhi, tag, nspecial, special, eflag, vflag, eatom,
                        vatom, host_start, ilist, jnum, cpu_time, success,
                        host_q, boxlo, prd);
}  
			
void ljcd_gpu_compute(const int ago, const int inum_full, const int nall,
                     double **host_x, int *host_type, int *ilist, int *numj,
                     int **firstneigh, const bool eflag, const bool vflag,
                     const bool eatom, const bool vatom, int &host_start,
                     const double cpu_time, bool &success, double *host_q,
                     const int nlocal, double *boxlo, double *prd) {
  LJCDMF.compute(ago,inum_full,nall,host_x,host_type,ilist,numj,firstneigh,eflag,
                vflag,eatom,vatom,host_start,cpu_time,success,host_q,
                nlocal,boxlo,prd);
}

double ljcd_gpu_bytes() {
  return LJCDMF.host_memory_usage();
}


