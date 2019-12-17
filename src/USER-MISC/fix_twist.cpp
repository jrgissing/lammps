/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Andrew Jewett (jewett@caltech.edu)
   adapted from fix_restrain.cpp (by Craig Tenney and Andres Jaramillo-Botero),
   and fix_ave_atom.cpp
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include <cmath>
#include <set>
using namespace std;
#include "fix_twist.h"
#include "atom.h"
#include "force.h"
#include "update.h"
#include "domain.h"
#include "comm.h"
#include "respa.h"
#include "input.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"

#include "math_const.h"
#include "math_extra.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;
using namespace MathExtra;

static const int DELTA=1;

/* ---------------------------------------------------------------------- */

FixTwist::FixTwist(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 4)
    error->all(FLERR,"Illegal fix twist command");

  scalar_flag = 1;
  global_freq = 1;
  extscalar = 1;
  peratom_flag = 1;
  create_attribute = 1;
  phi_prev = NULL;
  grow_arrays(atom->nmax);  // this should allocate phi_prev[]
  atom->add_callback(0);    // ? (I copied this from fix_ave_atom.cpp -Andrew)
  for (tagint i = 0; i < atom->nlocal; i++) // now fill with impossible values
    phi_prev[i] = 0.0;

  // parse args

  ntwist = maxtwist = 0;
  ids = NULL;
  k = NULL;
  param_start = param_stop = NULL;


  int iarg = 3;
  while (iarg < narg) {
    if (ntwist == maxtwist) {
      maxtwist += DELTA;
      memory->grow(ids,maxtwist,sizeof(tagint),"twist:ids");
      memory->grow(k,maxtwist,"twist:k");
      memory->grow(param_start,maxtwist,"twist:param_start");
      memory->grow(param_stop,maxtwist,"twist:param_stop");
    }

    if ((strcmp(arg[iarg],"torque") == 0) ||
	(strcmp(arg[iarg],"energy_per_turn") == 0)) {
      if (iarg+7 > narg) 
          error->all(FLERR,"Illegal fix twist command");
      ids[ntwist][0] = force->inumeric(FLERR,arg[iarg+1]);
      ids[ntwist][1] = force->inumeric(FLERR,arg[iarg+2]);
      ids[ntwist][2] = force->inumeric(FLERR,arg[iarg+3]);
      ids[ntwist][3] = force->inumeric(FLERR,arg[iarg+4]);

      param_start[ntwist]= force->numeric(FLERR,arg[iarg+5]);//energy per radian
      if (strcmp(arg[iarg],"energy_per_turn") == 0) 
	param_start[ntwist] /= MY_2PI;

      //phi_prev[ntwist] = MY_2PI*force->inumeric(FLERR,arg[iarg+6]);
      double phi_prev_in_radians = force->numeric(FLERR,arg[iarg+6]); 
      phi_prev[ntwist] = MY_2PI*floor((phi_prev_in_radians+MY_PI) / MY_2PI);
      //param_stop[ntwist]  = force->numeric(FLERR,arg[iarg+6]);
      param_stop[ntwist] = param_start[ntwist]; //(gets ignored later anyway)
      k[ntwist] = -1000.0; //any negative value will do

      tagint i_central_atom = ids[ntwist][1];

      set<tagint> reserved_atoms;
      if (reserved_atoms.find(i_central_atom) != reserved_atoms.end()) {
        char err_msg[1024];
        //stringstream err_msg("Illegal fix twist command: \n");
        //err_msg << 
        //  "    No atom can appear as the central atom (2nd atom) in the\n"
        //  "    list of twist interactions more than once.\n"
        //  "    (However atom " <<i_central_atom<< " has already been used.\n" 
        //  "     in the " << phi_prev[i_central_atom] << "th interaction.)\n";
        //error->all(FLERR,err_msg.str().c_str());
        sprintf(err_msg,
                "Illegal fix twist command:\n"
                "    No atom can appear as the central atom (2nd atom) in the\n"
                "    list of twist interactions more than once.\n"
                "    However atom " TAGINT_FORMAT " has already been used\n"
                "    before the %dth interaction.)\n",
                i_central_atom, ntwist);
        error->all(FLERR,err_msg);
      }
      reserved_atoms.insert(i_central_atom);

      iarg += 7;
    }
    else if (strcmp(arg[iarg],"rate") == 0) {
      if (iarg+9 > narg) 
        error->all(FLERR,"Illegal fix twist command");
      ids[ntwist][0] = force->inumeric(FLERR,arg[iarg+1]);
      ids[ntwist][1] = force->inumeric(FLERR,arg[iarg+2]);
      ids[ntwist][2] = force->inumeric(FLERR,arg[iarg+3]);
      ids[ntwist][3] = force->inumeric(FLERR,arg[iarg+4]);
      k[ntwist] = force->numeric(FLERR,arg[iarg+5]);
      if (k[ntwist] < 0.0) {
	error->all(FLERR,"Illegal fix twist command: k must be >= 0.0");
      }
      param_start[ntwist] = force->numeric(FLERR,arg[iarg+6]);
      param_stop[ntwist]  = force->numeric(FLERR,arg[iarg+7]);
      param_start[ntwist] *= MY_PI/180.0;
      param_stop[ntwist] *= MY_PI/180.0;
      //phi_prev[ntwist] = MY_2PI*force->inumeric(FLERR,arg[iarg+8]);
      double phi_prev_in_radians = force->numeric(FLERR,arg[iarg+8]); 
      phi_prev[ntwist] = MY_2PI*floor((phi_prev_in_radians+MY_PI) / MY_2PI);
      iarg += 9;
    }
    else
      error->all(FLERR,"Illegal fix twist command");

    ntwist++;
  }

  // require atom map to lookup atom IDs

  if (atom->map_style == 0)
    error->all(FLERR,"Fix twist requires an atom map, see atom_modify");
}

/* ---------------------------------------------------------------------- */

FixTwist::~FixTwist()
{
  memory->destroy(ids);
  memory->destroy(param_start);
  memory->destroy(param_stop);
  memory->destroy(k);
}

/* ---------------------------------------------------------------------- */

int FixTwist::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixTwist::init()
{
  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixTwist::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixTwist::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTwist::post_force(int vflag)
{
  energy = 0.0;

  for (tagint m = 0; m < ntwist; m++)
    compute_force(m);
}

/* ---------------------------------------------------------------------- */

void FixTwist::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixTwist::min_post_force(int vflag)
{
  post_force(vflag);
}


// --------------------------------------------
// ------- Calculate the dihedral angle -------
// --------------------------------------------
static const int g_dim=3;

double 
FixTwist::Phi(double const *x1, //array holding x,y,z coords atom 1
              double const *x2, // :       :      :      :        2
              double const *x3, // :       :      :      :        3
              double const *x4, // :       :      :      :        4
              Domain *domain, //<-periodic boundary information
              double *vb12, // will store x2-x1
              double *vb23, // will store x3-x2
              double *vb34, // will store x4-x3
              double *n123, // will store normal to plane x1,x2,x3
              double *n234) // will store normal to plane x2,x3,x4
{

  for (int d=0; d < g_dim; ++d) {
    vb12[d] = x2[d] - x1[d]; // 1st bond
    vb23[d] = x3[d] - x2[d]; // 2nd bond
    vb34[d] = x4[d] - x3[d]; // 3rd bond
  }

  //Consider periodic boundary conditions:
  domain->minimum_image(vb12[0],vb12[1],vb12[2]);
  domain->minimum_image(vb23[0],vb23[1],vb23[2]);
  domain->minimum_image(vb34[0],vb34[1],vb34[2]);

  //--- Compute the normal to the planes formed by atoms 1,2,3 and 2,3,4 ---

  cross3(vb23, vb12, n123);        // <- n123=vb23 x vb12
  cross3(vb23, vb34, n234);        // <- n234=vb23 x vb34

  norm3(n123);
  norm3(n234);

  double cos_phi = -dot3(n123, n234);

  if (cos_phi > 1.0)
    cos_phi = 1.0;
  else if (cos_phi < -1.0)
    cos_phi = -1.0;

  double phi = acos(cos_phi);

  if (dot3(n123, vb34) > 0.0) {
    phi = -phi;   //(Note: Negative dihedral angles are possible only in 3-D.)
    phi += MY_2PI; //<- This insure phi is always in the range 0 to 2*PI
  }
  return phi;
} // Phi()



/* ----------------------------------------------------------------------
   apply torque to a dihedral angle defined by 4 atoms (selected by m)
---------------------------------------------------------------------- */

void FixTwist::compute_force(tagint m)
{
  double f1[3],f2[3],f3[3],f4[3];

  double **x = atom->x;
  double **f = atom->f;
  tagint nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  double frac_time = (update->endstep - update->beginstep);
  if (frac_time != 0.0)
    frac_time = (update->ntimestep - update->beginstep) / frac_time;

  tagint i1 = atom->map(ids[m][0]);
  tagint i2 = atom->map(ids[m][1]);
  tagint i3 = atom->map(ids[m][2]);
  tagint i4 = atom->map(ids[m][3]);
  tagint icentral_atom = i2;

  // newton_bond on: only processor owning i2 computes force
  // newton_bond off: only processors owning any of i1-i4 computes force

  if (newton_bond) {
    if (i2 == -1 || i2 >= nlocal)
      return;
    if (i1 == -1 || i3 == -1 || i4 == -1) {
      char err_msg[1024];
      sprintf(err_msg, "Twist atoms " 
             TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT
              " missing on proc %d at step " BIGINT_FORMAT,
              ids[m][0],ids[m][1],ids[m][2],ids[m][3],
              comm->me,update->ntimestep);
      error->one(FLERR,err_msg);
    }
  }
  else {
    if ((i1 == -1 || i1 >= nlocal) && (i2 == -1 || i2 >= nlocal) &&
        (i3 == -1 || i3 >= nlocal) && (i4 == -1 || i3 >= nlocal))
      return;
    if (i1 == -1 || i2 == -1 || i3 == -1 || i4 == -1) {
      char err_msg[1024];
      sprintf(err_msg, "Twist atoms " 
             TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT " " TAGINT_FORMAT
              " missing on proc %d at step " BIGINT_FORMAT,
              ids[m][0],ids[m][1],ids[m][2],ids[m][3],
              comm->me,update->ntimestep);
      error->one(FLERR,err_msg);
    }
  }


  // The dihedral angle "phi" is the angle between n123 and n234
  // the planes defined by atoms i1,i2,i3, and i2,i3,i4.
  //
  // Definitions of vectors: vb12, vb23, vb34, perp12on23
  //                         proj12on23, perp43on32, proj43on32
  //
  //  Note: The positions of the 4 atoms are labeled x[i1], x[i2], x[i3], x[i4]
  //        (which are also vectors)
  //
  //             proj12on23                          proj34on23
  //             --------->                         ----------->
  //                           .
  //                          .
  //                         .
  //                  x[i2] .                       x[i3]
  //    .                __@----------vb23-------->@ . . . .           .
  //   /|\                /|                        \                  |
  //    |                /                           \                 |
  //    |               /                             \                |
  // perp12vs23        /                               \               |
  //    |             /                                 \          perp34vs23
  //    |          vb12                                  \             |
  //    |           /                                   vb34           |
  //    |          /                                       \           |
  //    |         /                                         \          |
  //    |        /                                           \         |
  //            @                                             \        |
  //                                                          _\|     \|/
  //         x[i1]                                              @
  //
  //                                                           x[i4]
  //

  // 1st bond
  double vb12[g_dim]; // displacement vector from atom i1 towards atom i2
  //     vb12[d]       = x[i2][d] - x[i1][d]      (for d=0,1,2)

  // 2nd bond:
  double vb23[g_dim]; // displacement vector from atom i2 towards atom i3
  //     vb23[d]       = x[i3][d] - x[i2][d]      (for d=0,1,2)

  // 3rd bond:
  double vb34[g_dim]; // displacement vector from atom i3 towards atom i4
  //     vb34[d]       = x[i4][d] - x[i3][d]      (for d=0,1,2)

  //  n123 & n234: These two unit vectors are normal to the planes
  //               defined by atoms 1,2,3 and 2,3,4.
  double n123[g_dim]; //n123=vb23 x vb12 / |vb23 x vb12|  ("x" is cross product)
  double n234[g_dim]; //n234=vb23 x vb34 / |vb23 x vb34|  ("x" is cross product)

  double proj12on23[g_dim];
  //    proj12on23[d] = (vb23[d]/|vb23|) * dot3(vb12,vb23)/|vb12|*|vb23|
  double proj34on23[g_dim];
  //    proj34on23[d] = (vb34[d]/|vb23|) * dot3(vb34,vb23)/|vb34|*|vb23|
  double perp12on23[g_dim];
  //    perp12on23[d] = v12[d] - proj12on23[d]
  double perp34on23[g_dim];
  //    perp34on23[d] = v34[d] - proj34on23[d]



  // ------ Step 1: Compute the dihedral angle "phi" ------
  //

  // Phi() calculates the dihedral angle.
  // This function also calculates the vectors:
  // vb12, vb23, vb34, n123, and n234, which we will need later.



  double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain,
                   vb12, vb23, vb34, n123, n234);

  // ------ Step 2: Compute the gradient of phi with atomic position: ------
  //
  // Gradient variables:
  //
  // dphi_dx1, dphi_dx2, dphi_dx3, dphi_dx4 are the gradients of phi with
  // respect to the atomic positions of atoms i1, i2, i3, i4, respectively.
  // As an example, consider dphi_dx1.  The d'th element is:
  double dphi_dx1[g_dim]; //                 d phi
  double dphi_dx2[g_dim]; // dphi_dx1[d] = ----------    (partial derivatives)
  double dphi_dx3[g_dim]; //               d x[i1][d]
  double dphi_dx4[g_dim]; //where d=0,1,2 corresponds to x,y,z  (if g_dim==3)

  double dot123     = dot3(vb12, vb23);
  double dot234     = dot3(vb23, vb34);
  double L23sqr     = dot3(vb23, vb23);
  double L23        = sqrt(L23sqr);   // (central bond length)
  double inv_L23sqr = 0.0;
  double inv_L23    = 0.0;
  if (L23sqr != 0.0) {
    inv_L23sqr = 1.0 / L23sqr;
    inv_L23 = 1.0 / L23;
  }
  double neg_inv_L23        = -inv_L23;
  double dot123_over_L23sqr = dot123 * inv_L23sqr;
  double dot234_over_L23sqr = dot234 * inv_L23sqr;

  for (int d=0; d < g_dim; ++d) {
    // See figure above for a visual definitions of these vectors:
    proj12on23[d] = vb23[d] * dot123_over_L23sqr;
    proj34on23[d] = vb23[d] * dot234_over_L23sqr;
    perp12on23[d] = vb12[d] - proj12on23[d];
    perp34on23[d] = vb34[d] - proj34on23[d];
  }


  // --- Compute the gradient vectors dphi/dx1 and dphi/dx4: ---

  // These two gradients point in the direction of n123 and n234,
  // and are scaled by the distances of atoms 1 and 4 from the central axis.
  // Distance of atom 1 to central axis:
  double perp12on23_len = sqrt(dot3(perp12on23, perp12on23));
  // Distance of atom 4 to central axis:
  double perp34on23_len = sqrt(dot3(perp34on23, perp34on23));

  double inv_perp12on23 = 0.0;
  if (perp12on23_len != 0.0) inv_perp12on23 = 1.0 / perp12on23_len;
  double inv_perp34on23 = 0.0;
  if (perp34on23_len != 0.0) inv_perp34on23 = 1.0 / perp34on23_len;

  for (int d=0; d < g_dim; ++d) {
    dphi_dx1[d] = n123[d] * inv_perp12on23;
    dphi_dx4[d] = n234[d] * inv_perp34on23;
  }

  // --- Compute the gradient vectors dphi/dx2 and dphi/dx3: ---
  //
  // This is more tricky because atoms 2 and 3 are shared by both planes
  // 123 and 234 (the angle between which defines "phi").  Moving either
  // one of these atoms effects both the 123 and 234 planes
  // Both the 123 and 234 planes intersect with the plane perpendicular to the
  // central bond axis (vb23).  The two lines where these intersections occur
  // will shift when you move either atom 2 or atom 3.  The angle between
  // these lines is the dihedral angle, phi.  We can define four quantities:
  // dphi123_dx2 is the change in "phi" due to the movement of the 123 plane
  //             ...as a result of moving atom 2.
  // dphi234_dx2 is the change in "phi" due to the movement of the 234 plane
  //             ...as a result of moving atom 2.
  // dphi123_dx3 is the change in "phi" due to the movement of the 123 plane
  //             ...as a result of moving atom 3.
  // dphi234_dx3 is the change in "phi" due to the movement of the 234 plane
  //             ...as a result of moving atom 3.

  double proj12on23_len = dot123 * inv_L23;
  double proj34on23_len = dot234 * inv_L23;
  // Interpretation:
  //The magnitude of "proj12on23_len" is the length of the proj12on23 vector.
  //The sign is positive if it points in the same direction as the central
  //bond (vb23).  Otherwise it is negative.  The same goes for "proj34on23".
  //(In the example figure in the comment above, both variables are positive.)

  // The forumula used in the 8 lines below explained separately:
  //   "supporting_information/doc/gradient_formula_explanation/"
  double dphi123_dx2_coef = neg_inv_L23 * (L23 + proj12on23_len);
  double dphi234_dx2_coef = inv_L23 * proj34on23_len;

  double dphi234_dx3_coef = neg_inv_L23 * (L23 + proj34on23_len);
  double dphi123_dx3_coef = inv_L23 * proj12on23_len;

  for (int d=0; d < g_dim; ++d) {
    // Recall that the n123 and n234 plane normal vectors are proportional to
    // the dphi/dx1 and dphi/dx2 gradients vectors
    // It turns out we can save slightly more CPU cycles by expressing
    // dphi/dx2 and dphi/dx3 as linear combinations of dphi/dx1 and dphi/dx2
    // which we computed already (instead of n123 & n234).
    dphi_dx2[d] = dphi123_dx2_coef*dphi_dx1[d] + dphi234_dx2_coef*dphi_dx4[d];
    dphi_dx3[d] = dphi123_dx3_coef*dphi_dx1[d] + dphi234_dx3_coef*dphi_dx4[d];
  }

  // ----- Step 3: Figure out any additional offsets to phi (n*2*pi) -----
  // -----         and correct for discontinuous jumps in phi        -----

  // The function Phi() returns numbers between -Pi and Pi (-180 and 180)
  // If our previous phi angle 179.5 degrees, between that step and this 
  // one, phi increased by 1 degree, then our new phi angle will
  // be -179.5 degrees.  If torque is positive, then will appear as 
  // though the energy increased by almost torque*2PI because the angle
  // suddenly decreased by almost 2PI.  However in reality the angle 
  // increased by a tiny amount.  So we check the magnitude of the change 
  // in angle between invocations to try and correct for this issue.

  double phi_prev_ntwists = floor(phi_prev[icentral_atom] / MY_2PI);
  double phi_prev_ntwists_x_2PI = phi_prev_ntwists * MY_2PI;
  double phi_prev_0_to_2pi = //phi_prev periodically mapped between 0 2pi.
    phi_prev[icentral_atom] - phi_prev_ntwists_x_2PI;
                            //this number should be similar to the current phi
  double delta_phi = phi - phi_prev_0_to_2pi;
  phi += phi_prev_ntwists_x_2PI;
  if (delta_phi < -MY_PI)
    phi += MY_2PI;
  else if (delta_phi > MY_PI)
    phi -= MY_2PI;

  phi_prev[icentral_atom] = phi;


  // ----- Step 4: Calculate the energy and force in the phi direction -----

  double U, dUdphi;

  if (k[m] > 0.0) {
    // If k[m] >= 0, then we assume the user wants to spin these atoms
    // at a constant rate determined by param_start[m] and param_stop[m]
    double phi0 = param_start[m] + frac_time * (param_stop[m]-param_start[m]);
    // (In most usage scenarios, I expect that param_start[m] != param_stop[m]
    //  Otherwise, most users would probably just use "fix restrain" instead.)

    //U = 0.5*k[m]*(1.0-cos(phi - phi0));
    //dUdphi = 0.5*k[m]*sin(phi - phi0);
    double Dphi = phi-phi0;
    U = 0.5*k[m]*Dphi*Dphi;
    dUdphi = k[m]*Dphi;

    energy += U; //Note:Energy is not conserved if phi0 is changing, since you
                 //     are essentially forcing the atoms to spin at some rate.
                 //     It's not clear that calculating the energy is meaningful
  }
  else {
    // If k[m] < 0.0, then we assume the user wants to to apply a constant
    // torque to these atoms. 

    double torque = param_start[m];

    dUdphi = -torque;

    energy += -torque*phi;
  }

  // ----- Step 4: Calculate the force direction in real space -----

  //          d U          d U      d phi
  // -f  =   -----   =    -----  *  -----
  //          d x         d phi      d x
  for(int d=0; d < g_dim; ++d) {
    f1[d] = -dUdphi * dphi_dx1[d];
    f2[d] = -dUdphi * dphi_dx2[d];
    f3[d] = -dUdphi * dphi_dx3[d];
    f4[d] = -dUdphi * dphi_dx4[d];
  }

  // apply force to each of 4 atoms

  if (newton_bond || i1 < nlocal) {
    f[i1][0] += f1[0];
    f[i1][1] += f1[1];
    f[i1][2] += f1[2];
  }

  if (newton_bond || i2 < nlocal) {
    f[i2][0] += f2[0];
    f[i2][1] += f2[1];
    f[i2][2] += f2[2];
  }

  if (newton_bond || i3 < nlocal) {
    f[i3][0] += f3[0];
    f[i3][1] += f3[1];
    f[i3][2] += f3[2];
  }

  if (newton_bond || i4 < nlocal) {
    f[i4][0] += f4[0];
    f[i4][1] += f4[1];
    f[i4][2] += f4[2];
  }

  //I suppose the "ev_tally()" function is not used in fixes, but I keep it here
  //to remind me I should some day try to figure out how this fix effects the 
  //virial. (You could possibly have a lot of dihedrals controlled by this fix.)
  //if (evflag) {
  //  ev_tally(i1,i2,i3,i4,
  //           nlocal,
  //           newton_bond,
  //           energy,
  //           0,
  //           f1,f3,f4,
  //           vb12[0],vb12[1],vb12[2],
  //           vb23[0],vb23[1],vb23[2],
  //           vb34[0],vb34[1],vb34[2]);
  //}

} //FixTwist::compute_force()



/* ----------------------------------------------------------------------
   potential energy of added force
------------------------------------------------------------------------- */

double FixTwist::compute_scalar()
{
  MPI_Allreduce(&energy,&energy_all,1,MPI_DOUBLE,MPI_SUM,world);
  return energy_all;
}





/* ---------------------------------------------------------------------- */

void FixTwist::setup_pre_force(int vflag) {
  FixTwist::calc_phi_prev();
}

void FixTwist::setup_min_pre_force(int vflag) {
  FixTwist::calc_phi_prev();
}

// The next two functions need to be defined so that LAMMPS does not ignore
// the setup_pre_force() and setup_min_pre_force() functions.  (IE. so that
// "n_pre_force" and "n_min_pre_force" are both > 0.  See modify.cpp)
void FixTwist::pre_force(int vflag) {
  //do nothing
}
void FixTwist::min_pre_force(int vflag) {
  //do nothing
}


void FixTwist::calc_phi_prev() {
  double **x = atom->x;
  tagint nlocal = atom->nlocal;
  tagint newton_bond = force->newton_bond;
  for (tagint m = 0; m < ntwist; m++) {
    tagint i1 = atom->map(ids[m][0]);
    tagint i2 = atom->map(ids[m][1]);
    tagint i3 = atom->map(ids[m][2]);
    tagint i4 = atom->map(ids[m][3]);
    tagint icentral_atom = i2;
    // newton_bond on: only processor owning i2 computes force
    // newton_bond off: only processors owning any of i1-i4 computes force
    if (newton_bond) {
      if (i2 == -1 || i2 >= nlocal) return;
      if (i1 == -1 || i3 == -1 || i4 == -1) return;
    }
    else {
      if ((i1 == -1 || i1 >= nlocal) && (i2 == -1 || i2 >= nlocal) &&
          (i3 == -1 || i3 >= nlocal) && (i4 == -1 || i3 >= nlocal))
        return;
      if (i1 == -1 || i2 == -1 || i3 == -1 || i4 == -1)  return;
    }
    double vb12[g_dim]; // displacement vector from atom i1 towards atom i2
    double vb23[g_dim]; // displacement vector from atom i2 towards atom i3
    double vb34[g_dim]; // displacement vector from atom i3 towards atom i4
    double n123[g_dim]; //n123=vb23 x vb12/|vb23 x vb12|  ("x" is cross product)
    double n234[g_dim]; //n234=vb23 x vb34/|vb23 x vb34|  ("x" is cross product)
    double phi = Phi(x[i1], x[i2], x[i3], x[i4], domain,
                     vb12, vb23, vb34, n123, n234);
    phi_prev[icentral_atom] = phi;
  } //for (tagint m = 0; m < ntwist; m++)
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixTwist::memory_usage()
{
  double bytes;
  bytes = atom->nmax * sizeof(double);   //for the phi_prev[] array
  bytes += maxtwist * sizeof(tagint);    //for the ids[] array (could be big)
  bytes += maxtwist * 3*sizeof(double);  //for param_start[], param_stop[], k[]

  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixTwist::grow_arrays(int nmax)
{
  memory->grow(phi_prev,nmax,"fix_twist:phi_prev");
  //Alternately, if phi_prev was of type std::vector<double>, we could try:
  //phi_prev.resize(nmax, 0);  (instead of using memory->grow() ?)
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixTwist::copy_arrays(int i, int j, int delflag)
{
  //The following more general code comes from fix_ave_atom.cpp:
  //for (int m = 0; m < nvalues; m++)
  //  array[j][m] = array[i][m];
  //However, in this simple fix, we only need to store 1 number per atom:
  phi_prev[j] = phi_prev[i];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixTwist::pack_exchange(int i, double *buf)
{
  //Debugging:
  //char err_msg[2048];
  //sprintf(err_msg, "Good news! FixTwist::pack_exchange(%d,%p) invoked.\n"
  //        "Exiting...\n", i, buf);
  //error->one(FLERR, err_msg);

  buf[0] = phi_prev[i];
  //Incidentally, the following more general code comes from fix_ave_atom.cpp:
  //for (int m = 0; m < nvalues; m++)
  //  buf[m] = array[i][m];
  //return nvalues;
  //However, in this simple fix, we only need to store 1 number per atom:

  return 1;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixTwist::unpack_exchange(int nlocal, double *buf)
{
  //Debugging:
  //char err_msg[2048];
  //sprintf(err_msg, "Good news! FixTwist::unpack_exchange(%d,%p) invoked.\n"
  //        "Exiting...\n", nlocal, buf);
  //error->one(FLERR,err_msg);

  phi_prev[nlocal] = buf[0];
  //The following more general code comes from fix_ave_atom.cpp:
  //for (int m = 0; m < nvalues; m++)
  //  array[nlocal][m] = buf[m];
  //return nvalues;
  //However, in this simple fix, we only need to store 1 number per atom:

  return 1;
}
