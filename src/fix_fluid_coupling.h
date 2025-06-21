/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

    DEM simulation engine, released by
    DCS Computing Gmbh, Linz, Austria
    http://www.dcs-computing.com, office@dcs-computing.com

    LIGGGHTS® is part of CFDEM®project:
    http://www.liggghts.com | http://www.cfdem.com

    Core developer and main author:
    Christoph Kloss, christoph.kloss@dcs-computing.com

    LIGGGHTS® is open-source, distributed under the terms of the GNU Public
    License, version 2 or later. It is distributed in the hope that it will
    be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. You should have
    received a copy of the GNU General Public License along with LIGGGHTS®.
    If not, see http://www.gnu.org/licenses . See also top-level README
    and LICENSE files.

    LIGGGHTS® and CFDEM® are registered trade marks of DCS Computing GmbH,
    the producer of the LIGGGHTS® software and the CFDEM®coupling software
    See http://www.cfdem.com/terms-trademark-policy for details.

-------------------------------------------------------------------------
    Contributing author and copyright for this file:
    Robin Walloner
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(fluid_coupling,FixFluidCoupling)

#else

#ifndef LMP_FIX_FLUID_COUPLING_H
#define LMP_FIX_FLUID_COUPLING_H

#include "fix.h"
#include "compute.h"

namespace LAMMPS_NS {

class FixFluidCoupling : public Fix {
 public:
  FixFluidCoupling(class LAMMPS *, int, char **);
  ~FixFluidCoupling();
  int setmask();
  void init();
  void setup(int);
  void post_force(int);
  void final_integrate();

  void read_fluid_velocity(double);
  void write_particle_force();
  void compute_volume_fraction();
  void compute_force();

  double compute_array(int,int);
  void compute_array_atom();
  void grow_arrays(int);
  void set_arrays(int);
  void copy_arrays(int, int,int);
  double memory_usage();  

 private:
  int drag_law;
  int coupling_type;
  int gravity_flag;
  int buoyancy_flag;

  double rho_fluid;
  double mu_fluid;
  double *gravity;

  Compute *c_voronoi;

  double nmax;
  double **v_fluid;
  double *volume;
  double *volfrac;
  double *drag_coeff;
  double *reynolds;
  double **f_drag;
  double **f_buoyancy;
  double **f_gravity;
  double **f_total;
  double **expl_momentum;
  double *impl_momentum;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Region ID for fix addforce does not exist

Self-explanatory.

E: Variable name for fix addforce does not exist

Self-explanatory.

E: Variable for fix addforce is invalid style

Self-explanatory.

E: Cannot use variable energy with constant force in fix addforce

This is because for constant force, LAMMPS can compute the change
in energy directly.

E: Must use variable energy with fix addforce

Must define an energy vartiable when applyting a dynamic
force during minimization.

*/
