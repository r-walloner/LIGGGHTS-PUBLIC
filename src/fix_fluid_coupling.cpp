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

#include <cassert>
#include <math.h>
#include <precice/preciceC.h>
#include "fix_fluid_coupling.h"
#include "fix_gravity.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

enum
{
  DRAG_STOKES,
  DRAG_GIDASPOW,
  DRAG_KOCH_HILL,
  DRAG_ZHAO_SHAN,
};

enum
{
  COUPLING_MOMENTUM_SEMI_IMPLICIT,
  COUPLING_FORCE,
};

/* ----------------------------------------------------------------------
   construct fluid_coupling fix and parse arguments from input script
------------------------------------------------------------------------- */

FixFluidCoupling::FixFluidCoupling(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp, narg, arg)
{
  nmax = 0;
  v_fluid = NULL;
  volume = NULL;
  volfrac = NULL;
  drag_coeff = NULL;
  reynolds = NULL;
  f_drag = NULL;
  f_buoyancy = NULL;
  f_gravity = NULL;
  f_fluid_total = NULL;
  expl_momentum = NULL;
  impl_momentum = NULL;

  memory->create(gravity, 3, "FixFluidCoupling:gravity");

  peratom_flag = 1;
  peratom_freq = 1;
  size_peratom_cols = 27;
  array_atom = NULL;

  // Parse arguments
  if (narg != 9)
    error->all(FLERR, "Illegal fix fluid_coupling command");

  drag_law = -1;
  if (strcmp(arg[3], "stokes") == 0)
    drag_law = DRAG_STOKES;
  else if (strcmp(arg[3], "gidaspow") == 0)
    drag_law = DRAG_GIDASPOW;
  else if (strcmp(arg[3], "koch_hill") == 0)
    drag_law = DRAG_KOCH_HILL;
  else if (strcmp(arg[3], "zhao_shan") == 0)
    drag_law = DRAG_ZHAO_SHAN;
  if (drag_law == -1)
    error->all(FLERR, "Illegal fluid drag law specified");

  coupling_type = -1;
  if (strcmp(arg[4], "momentum_semi_implicit") == 0)
    coupling_type = COUPLING_MOMENTUM_SEMI_IMPLICIT;
  else if (strcmp(arg[4], "force") == 0)
    coupling_type = COUPLING_FORCE;
  if (coupling_type == -1)
    error->all(FLERR, "Illegal coupling type specified");

  rho_fluid = force->numeric(FLERR, arg[5]);
  mu_fluid = force->numeric(FLERR, arg[6]);
  if (rho_fluid <= 0.0)
    error->all(FLERR, "Fluid density must be positive");
  if (mu_fluid <= 0.0)
    error->all(FLERR, "Fluid viscosity must be positive");

  gravity_flag = force->inumeric(FLERR, arg[7]);
  buoyancy_flag = force->inumeric(FLERR, arg[8]);
}

/* ----------------------------------------------------------------------
    destruct fluid_coupling fix and free allocated memory
------------------------------------------------------------------------- */

FixFluidCoupling::~FixFluidCoupling()
{
  memory->destroy(array_atom);
  memory->destroy(v_fluid);
  memory->destroy(volume);
  memory->destroy(volfrac);
  memory->destroy(drag_coeff);
  memory->destroy(reynolds);
  memory->destroy(f_drag);
  memory->destroy(f_buoyancy);
  memory->destroy(f_gravity);
  memory->destroy(f_fluid_total);
  memory->destroy(expl_momentum);
  memory->destroy(impl_momentum);
  memory->destroy(gravity);
}

/* ----------------------------------------------------------------------
   determines when the fix is called during the timestep
------------------------------------------------------------------------- */

int FixFluidCoupling::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= FINAL_INTEGRATE;
  return mask;
}

/* ----------------------------------------------------------------------
   initialization before a run
------------------------------------------------------------------------- */

void FixFluidCoupling::init()
{
  // Find the computes for the voronoi tessellation and gravity
  c_voronoi = modify->find_compute_style_strict("voronoi/atom", 0);
  if (!c_voronoi)
    error->all(FLERR, "Cannot find voronoi/atom compute");

  if (buoyancy_flag || gravity_flag)
  {
    // TODO get gravity vector from input script
    gravity[0] = 0.0;
    gravity[1] = -9.81;
    gravity[2] = 0.0;
  }
  else
  {
    gravity[0] = 0.0;
    gravity[1] = 0.0;
    gravity[2] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   called after pair & molecular forces are computed and communicated
------------------------------------------------------------------------- */

void FixFluidCoupling::post_force(int vflag)
{
  if (atom->nlocal > nmax)
    grow_arrays(atom->nlocal);

  read_fluid_velocity(0);
  compute_volume_fraction();
  compute_force();
  apply_force();
  compute_array_atom();
}

/* ----------------------------------------------------------------------
   called at end of each timestep
------------------------------------------------------------------------- */

void FixFluidCoupling::final_integrate()
{
  read_fluid_velocity(update->dt);
  compute_volume_fraction();
  compute_force();
  write_particle_force();
}

/* ----------------------------------------------------------------------
   read the fluid velocity from preCICE into v_fluid array
------------------------------------------------------------------------- */

void FixFluidCoupling::read_fluid_velocity(double relative_read_time)
{
  precicec_mapAndReadData(
      "Fluid-Mesh", "Velocity",
      atom->nlocal, *atom->x, relative_read_time, *v_fluid);
}

/* ----------------------------------------------------------------------
   write the particle force / momentum contribution and the particle volume
   fraction to preCICE
------------------------------------------------------------------------- */

void FixFluidCoupling::write_particle_force()
{
  if (coupling_type == COUPLING_MOMENTUM_SEMI_IMPLICIT)
  {
    precicec_writeAndMapData(
        "Fluid-Mesh", "ExplicitMomentum",
        atom->nlocal, *atom->x, *expl_momentum);
    precicec_writeAndMapData(
        "Fluid-Mesh", "ImplicitMomentum",
        atom->nlocal, *atom->x, impl_momentum);
  }

  else if (coupling_type == COUPLING_FORCE)
  {
    precicec_writeAndMapData(
        "Fluid-Mesh", "DragForce",
        atom->nlocal, *atom->x, *f_drag);
    precicec_writeAndMapData(
        "Fluid-Mesh", "Alpha",
        atom->nlocal, *atom->x, volume);
  }
}

/* ----------------------------------------------------------------------
   compute the particle volume fraction into volfrac array
------------------------------------------------------------------------- */

void FixFluidCoupling::compute_volume_fraction()
{
  precicec_startProfilingSection("solver.advance.compute_volume_fraction");

  // Access per-atom voronoi cell volume from voronoi compute
  if (!(c_voronoi->invoked_flag & INVOKED_PERATOM))
  {
    c_voronoi->compute_peratom();
    c_voronoi->invoked_flag |= INVOKED_PERATOM;
  }
  double **voro_tess = c_voronoi->array_atom;
  if (!voro_tess)
    error->all(FLERR, "voronoi/atom compute does not provide a per-atom array");

  // Compute volume fraction
  for (int i = 0; i < atom->nlocal; i++)
  {
    if (!(atom->mask[i] & groupbit))
      continue;

    volume[i] = (4.0 / 3.0) * M_PI * pow(atom->radius[i], 3);
    volfrac[i] = volume[i] / voro_tess[i][0];

    assert(volume[i] >= 0.0 && "Particle volume is negative");
    assert(volfrac[i] >= 0.0 && "Particle volume fraction is negative");
    assert(volfrac[i] <= 1.0 && "Particle volume fraction exceeds 1.0");
  }

  precicec_stopLastProfilingSection();
}

/* ----------------------------------------------------------------------
   compute the interaction force and its momentum contribution to the fluid
   into f_drag, expl_momentum and impl_momentum arrays
------------------------------------------------------------------------- */

void FixFluidCoupling::compute_force()
{
  precicec_startProfilingSection("solver.advance.compute_force");

  double *v_rel = (double *)malloc(3 * sizeof(double));
  double mag_v_rel;
  double beta, F_0, F_3, X, diameter, volfrac_p, volfrac_f;

  for (int i = 0; i < atom->nlocal; i++)
  {
    if (!(atom->mask[i] & groupbit))
      continue;

    // relative particle velocity
    sub3(v_fluid[i], atom->v[i], v_rel);
    mag_v_rel = len3(v_rel);

    // particle diameter
    diameter = 2 * atom->radius[i];

    // volume fractions
    volfrac_p = volfrac[i];
    volfrac_f = 1 - volfrac[i];

    // drag force (using different drag laws)
    // --------------------------------------

    if (drag_law == DRAG_STOKES)
    {
      for (int d = 0; d < 3; d++)
        f_drag[i][d] = 3 * M_PI * diameter *
                       mu_fluid * volfrac_f * v_rel[d];
      // note: this does not support semi-implicit momentum coupling yet
    }

    else if (drag_law == DRAG_ZHAO_SHAN)
    {
      reynolds[i] = rho_fluid * diameter * mag_v_rel / mu_fluid;

      if (reynolds[i] == 0)
        drag_coeff[i] = 0;
      else
        drag_coeff[i] = .2924 * pow(9.06 / sqrt(reynolds[i]) + 1, 2);

      if (reynolds[i] == 0)
        X = 0; // drag does not matter if ralative velocity is zero
      else
        X = 3.7 - 0.65 * exp(-0.5 * pow(1.5 - log10(reynolds[i]), 2));

      for (int d = 0; d < 3; d++)
        f_drag[i][d] = .125 * rho_fluid * M_PI * pow(diameter, 2) *
                       drag_coeff[i] * pow(volfrac_f, -(X + 1)) * v_rel[d] * mag_v_rel;

      if (coupling_type == COUPLING_MOMENTUM_SEMI_IMPLICIT)
      {
        impl_momentum[i] = .125 * rho_fluid * M_PI * pow(diameter, 2) *
                           drag_coeff[i] * pow(volfrac_f, -(X + 1)) * mag_v_rel;
        for (int d = 0; d < 3; d++)
          expl_momentum[i][d] = impl_momentum[i] * atom->v[i][d];
      }
    }

    else if (drag_law == DRAG_GIDASPOW || drag_law == DRAG_KOCH_HILL)
    {
      reynolds[i] = volfrac_f * rho_fluid * mag_v_rel * diameter /
                 mu_fluid;

      if (drag_law == DRAG_GIDASPOW)
      {
        // drag coefficient
        if (reynolds[i] == 0)
          drag_coeff[i] = 0; // drag does not matter if ralative velocity is zero
        else if (reynolds[i] < 1000)
          drag_coeff[i] = 24.0 * (1.0 + 0.15 * pow(reynolds[i], 0.687)) / reynolds[i];
        else
          drag_coeff[i] = 0.44;

        // inter-phase momentum exchange coefficient beta
        if (volfrac_f < 0.8)
          beta = 150.0 * pow(volfrac_p, 2) * mu_fluid /
                     (volfrac_f * pow(diameter, 2)) +
                 1.75 * volfrac_p * rho_fluid * mag_v_rel / diameter;
        else
          beta = 3.0 * drag_coeff[i] * volfrac_f * pow(volfrac_p, 2) *
                 rho_fluid * mag_v_rel * pow(volfrac_f, -2.65) /
                 (2.0 * diameter);
      }

      else if (drag_law == DRAG_KOCH_HILL)
      {
        if (volfrac_p < 0.4)
          F_0 = (1 + 3 * sqrt(volfrac_p / 2) + 2.109 * volfrac_p * log(volfrac_p) + 16.14 * volfrac_p) /
                (1 + 0.681 * volfrac_p - 8.48 * pow(volfrac_p, 2) + 8.16 * pow(volfrac_p, 3));
        else
          F_0 = 10 * volfrac_p / pow(volfrac_f, 3);

        F_3 = 0.0673 + 0.0212 * volfrac_p + (0.0232 / pow(volfrac_f, 5));

        beta = 18 * mu_fluid * pow(volfrac_f, 2) * volfrac_p *
               (F_0 + 0.5 * F_3 * reynolds[i]) /
               pow(diameter, 2);
      }

      // drag force
      for (int d = 0; d < 3; d++)
        f_drag[i][d] = beta * volume[i] * v_rel[d] / volfrac_p;

      if (coupling_type == COUPLING_MOMENTUM_SEMI_IMPLICIT)
      {
        // momentum contribution to fluid
        impl_momentum[i] = beta * volume[i] / (volfrac_p * volfrac_f);
        for (int d = 0; d < 3; d++)
          expl_momentum[i][d] = impl_momentum[i] * atom->v[i][d];
      }
    }

    // gravity and buoyancy forces
    // ---------------------------
    for (int d = 0; d < 3; d++)
    {
      f_gravity[i][d] = gravity_flag * volume[i] * atom->density[i] * gravity[d];
      f_buoyancy[i][d] = buoyancy_flag * volume[i] * rho_fluid * -gravity[d];
    }
  }

  free(v_rel);

  precicec_stopLastProfilingSection();
}

/* ----------------------------------------------------------------------
   apply the computed forces to the particles
------------------------------------------------------------------------- */

void FixFluidCoupling::apply_force()
{
  for (int i = 0; i < atom->nlocal; i++)
  {
    for (int d = 0; d < 3; d++)
    {
      f_fluid_total[i][d] = f_drag[i][d] + f_buoyancy[i][d];
      atom->f[i][d] += f_fluid_total[i][d] + f_gravity[i][d];
    }
  }
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixFluidCoupling::grow_arrays(int nmax_new)
{
  memory->destroy(array_atom);
  memory->create(array_atom, nmax_new, size_peratom_cols, "FixFluidCoupling:array_atom");
  memory->grow(v_fluid, nmax_new, 3, "FixFluidCoupling:v_fluid");
  memory->grow(volume, nmax_new, "FixFluidCoupling:volume");
  memory->grow(volfrac, nmax_new, "FixFluidCoupling:volfrac");
  memory->grow(drag_coeff, nmax_new, "FixFluidCoupling:drag_coeff");
  memory->grow(reynolds, nmax_new, "FixFluidCoupling:reynolds");
  memory->grow(f_drag, nmax_new, 3, "FixFluidCoupling:f_drag");
  memory->grow(f_buoyancy, nmax_new, 3, "FixFluidCoupling:f_buoyancy");
  memory->grow(f_gravity, nmax_new, 3, "FixFluidCoupling:f_gravity");
  memory->grow(f_fluid_total, nmax_new, 3, "FixFluidCoupling:f_fluid_total");
  memory->grow(expl_momentum, nmax_new, 3, "FixFluidCoupling:expl_momentum");
  memory->grow(impl_momentum, nmax_new, "FixFluidCoupling:impl_momentum");
  nmax = nmax_new;
}

/* ----------------------------------------------------------------------
   return I,J array value for output to file
------------------------------------------------------------------------- */

double FixFluidCoupling::compute_array(int i, int j)
{
  if (j == 1) return v_fluid[i][0];
  else if (j == 2) return v_fluid[i][1];
  else if (j == 3) return v_fluid[i][2];
  else if (j == 4) return volume[i];
  else if (j == 5) return volfrac[i];
  else if (j == 6) return f_drag[i][0];
  else if (j == 7) return f_drag[i][1];
  else if (j == 8) return f_drag[i][2];
  else if (j == 9) return expl_momentum[i][0];
  else if (j == 10) return expl_momentum[i][1];
  else if (j == 11) return expl_momentum[i][2];
  else if (j == 12) return impl_momentum[i];
  else if (j == 13) return drag_coeff[i];
  else if (j == 14) return v_fluid[i][0] - atom->v[i][0];
  else if (j == 15) return v_fluid[i][1] - atom->v[i][1];
  else if (j == 16) return v_fluid[i][2] - atom->v[i][2];
  else if (j == 17) return 1 - volfrac[i];
  else if (j == 18) return reynolds[i];
  else if (j == 19) return f_buoyancy[i][0];
  else if (j == 20) return f_buoyancy[i][1];
  else if (j == 21) return f_buoyancy[i][2];
  else if (j == 22) return f_fluid_total[i][0];
  else if (j == 23) return f_fluid_total[i][1];
  else if (j == 24) return f_fluid_total[i][2];
  else if (j == 25) return f_gravity[i][0];
  else if (j == 26) return f_gravity[i][1];
  else if (j == 27) return f_gravity[i][2];
  else
    error->all(FLERR, "Invalid column index in FixFluidCoupling::compute_array");
  return 0.0; // unreachable, but avoids compiler warning
}

/* ----------------------------------------------------------------------
   prepare array_atom for output to file
------------------------------------------------------------------------- */

void FixFluidCoupling::compute_array_atom()
{
  for (int i = 0; i < nmax; i++)
    for (int j = 0; j < size_peratom_cols; j++)
      array_atom[i][j] = compute_array(i, j + 1);
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixFluidCoupling::set_arrays(int i)
{
  v_fluid[i][0] = 0.0;
  v_fluid[i][1] = 0.0;
  v_fluid[i][2] = 0.0;
  volume[i] = 0.0;
  volfrac[i] = 0.0;
  drag_coeff[i] = 0.0;
  reynolds[i] = 0.0;
  f_drag[i][0] = 0.0;
  f_drag[i][1] = 0.0;
  f_drag[i][2] = 0.0;
  f_buoyancy[i][0] = 0.0;
  f_buoyancy[i][1] = 0.0;
  f_buoyancy[i][2] = 0.0;
  f_gravity[i][0] = 0.0;
  f_gravity[i][1] = 0.0;
  f_gravity[i][2] = 0.0;
  f_fluid_total[i][0] = 0.0;
  f_fluid_total[i][1] = 0.0;
  f_fluid_total[i][2] = 0.0;
  expl_momentum[i][0] = 0.0;
  expl_momentum[i][1] = 0.0;
  expl_momentum[i][2] = 0.0;
  impl_momentum[i] = 0.0;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixFluidCoupling::copy_arrays(int i, int j, int delflag)
{
  v_fluid[j][0] = v_fluid[i][0];
  v_fluid[j][1] = v_fluid[i][1];
  v_fluid[j][2] = v_fluid[i][2];
  volume[j] = volume[i];
  volfrac[j] = volfrac[i];
  drag_coeff[j] = drag_coeff[i];
  reynolds[j] = reynolds[i];
  f_drag[j][0] = f_drag[i][0];
  f_drag[j][1] = f_drag[i][1];
  f_drag[j][2] = f_drag[i][2];
  f_buoyancy[j][0] = f_buoyancy[i][0];
  f_buoyancy[j][1] = f_buoyancy[i][1];
  f_buoyancy[j][2] = f_buoyancy[i][2];
  f_gravity[j][0] = f_gravity[i][0];
  f_gravity[j][1] = f_gravity[i][1];
  f_gravity[j][2] = f_gravity[i][2];
  f_fluid_total[j][0] = f_fluid_total[i][0];
  f_fluid_total[j][1] = f_fluid_total[i][1];
  f_fluid_total[j][2] = f_fluid_total[i][2];
  expl_momentum[j][0] = expl_momentum[i][0];
  expl_momentum[j][1] = expl_momentum[i][1];
  expl_momentum[j][2] = expl_momentum[i][2];
  impl_momentum[j] = impl_momentum[i];
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixFluidCoupling::memory_usage()
{
  double bytes = (23 + size_peratom_cols) * atom->nmax * sizeof(double);
  return bytes;
}
