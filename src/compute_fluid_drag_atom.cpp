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
#include <cmath>
#include <precice/preciceC.h>
#include "compute_fluid_drag_atom.h"
#include "atom.h"
#include "error.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace MathExtra;

enum
{
  DRAG_STOKES,
  DRAG_XIAO_SUN,
  DRAG_KOCH_HILL
};

/* ---------------------------------------------------------------------- */

ComputeFluidDragAtom::ComputeFluidDragAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) : Compute(lmp, iarg, narg, arg)
{
  if (iarg + 1 != narg)
    error->all(FLERR, "Illegal compute fluid_drag/atom command");

  peratom_flag = 1;
  size_peratom_cols = 3;

  c_v_fluid = NULL;
  c_vol_frac = NULL;

  drag_law = -1;
  if (strcmp(arg[iarg], "stokes") == 0)
    drag_law = DRAG_STOKES;
  else if (strcmp(arg[iarg], "xiao_sun") == 0)
    drag_law = DRAG_XIAO_SUN;
  else if (strcmp(arg[iarg], "koch_hill") == 0)
    drag_law = DRAG_KOCH_HILL;
  if (drag_law == -1)
    error->all(FLERR, "Illegal fluid drag law specified");
  iarg++;

  f_drag = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeFluidDragAtom::~ComputeFluidDragAtom()
{
  memory->destroy(f_drag);
}

/* ---------------------------------------------------------------------- */

void ComputeFluidDragAtom::init()
{
  // Find required computes
  c_v_fluid = modify->find_compute_id("v_fluid");
  if (!c_v_fluid)
    error->all(FLERR, "Cannot find compute v_fluid");

  c_vol_frac = modify->find_compute_style_strict("volume_fraction/atom", 0);
  if (!c_vol_frac)
    error->all(FLERR, "Cannot find volume_fraction/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeFluidDragAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow data array if necessary
  if (atom->nlocal > nmax)
  {
    memory->destroy(f_drag);
    nmax = atom->nlocal;
    memory->create(f_drag, nmax, size_peratom_cols, "compute:fluid_drag/atom");
    array_atom = f_drag;
  }

  // Access per-atom fluid velocity and particle volume fraction from compute
  if (!(c_v_fluid->invoked_flag & INVOKED_PERATOM))
  {
    c_v_fluid->compute_peratom();
    c_v_fluid->invoked_flag |= INVOKED_PERATOM;
  }
  double **v_fluid = c_v_fluid->array_atom;
  if (!v_fluid)
    error->all(FLERR, "Compute v_fluid does not provide a per-atom array");

  if (!(c_vol_frac->invoked_flag & INVOKED_PERATOM))
  {
    c_vol_frac->compute_peratom();
    c_vol_frac->invoked_flag |= INVOKED_PERATOM;
  }
  double *vol_frac = c_vol_frac->vector_atom;
  if (!vol_frac)
    error->all(FLERR, "volume_fraction/atom compute does not provide a per-atom array");

  // Assume that all particles have the same radius
  // TODO refactor this into its own compute to reduce code duplication
  //   - maybe we can use the compute property/atom command
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++)
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  // Compute particle diameter and volume
  const double diameter = 2 * radius;
  const double volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  // TODO get constants from global properties
  // air:
  const double rho_fluid = 1; // kg/m^3
  const double mu_fluid = 1.8e-5; // Pa.s
  // water:
  // const double rho_fluid = 1000; // kg/m^3
  // const double mu_fluid = 1e-3;  // Pa.s

  double beta, drag_coeff, reynolds, F_0, F_3;
  double *v_rel = (double *)malloc(3 * sizeof(double));

  double *expl_momentum = (double *)malloc(3 * sizeof(double));
  double *impl_momentum = (double *)malloc(sizeof(double));

  // Calculate and apply drag force
  for (int i = 0; i < atom->nlocal; i++)
  {
    if (!(atom->mask[i] & groupbit))
      continue;

    // particle velocity relative to fluid
    sub3(v_fluid[i], atom->v[i], v_rel);
    const double mag_v_rel = len3(v_rel);

    if (drag_law == DRAG_XIAO_SUN || drag_law == DRAG_KOCH_HILL)
    {
      // reynolds number
      reynolds = (1 - vol_frac[i]) * rho_fluid * mag_v_rel * diameter /
                 mu_fluid;

      if (drag_law == DRAG_XIAO_SUN)
      {
        // drag coefficient
        if (reynolds == 0)
          drag_coeff = 0; // drag is zero for stationary particles, so drag_coeff does not matter
        else if (reynolds < 1000)
          drag_coeff = 24.0 * (1.0 + 0.15 * pow(reynolds, 0.687)) / reynolds;
        else
          drag_coeff = 0.44;

        // beta
        if (1 - vol_frac[i] < 0.8)
          beta = 150.0 * pow(vol_frac[i], 2) * mu_fluid /
                    ((1 - vol_frac[i]) * pow(diameter, 2)) +
                1.75 * vol_frac[i] * rho_fluid * mag_v_rel / diameter;
        else
          beta = 3.0 * drag_coeff * (1 - vol_frac[i]) * pow(vol_frac[i], 2) *
                rho_fluid * mag_v_rel * pow(1 - vol_frac[i], -2.65) /
                (2.0 * diameter);
      }

      else if (drag_law == DRAG_KOCH_HILL)
      {
        reynolds = (1 - vol_frac[i]) * rho_fluid * mag_v_rel * diameter /
                  mu_fluid;

        if (vol_frac[i] < 0.4)
          F_0 = (1 + 3 * sqrt(vol_frac[i] / 2) + 2.109 * vol_frac[i] * log(vol_frac[i]) + 16.14 * vol_frac[i]) /
                (1 + 0.681 * vol_frac[i] - 8.48 * pow(vol_frac[i], 2) + 8.16 * pow(vol_frac[i], 3));
        else
          F_0 = 10 * vol_frac[i] / pow(1 - vol_frac[i], 3);

        F_3 = 0.0673 + 0.0212 * vol_frac[i] + (0.0232 / pow(1 - vol_frac[i], 5));

        beta = 18 * mu_fluid * pow(1 - vol_frac[i], 2) * vol_frac[i] *
              (F_0 + 0.5 * F_3 * reynolds) /
              pow(diameter, 2);
      }

      // drag force
      for (int d = 0; d < size_peratom_cols; d++)
        f_drag[i][d] = beta * volume * v_rel[d] / vol_frac[i];

      // TODO move this to its own compute
      *impl_momentum = beta * volume / (vol_frac[i] * (1 - vol_frac[i]));
      precicec_writeAndMapData("Fluid-Mesh", "ImplicitMomentum", 1, atom->x[i], impl_momentum);

      for (int d = 0; d < 3; d++)
        expl_momentum[d] = *impl_momentum * atom->v[i][d];
      precicec_writeAndMapData("Fluid-Mesh", "ExplicitMomentum", 1, atom->x[i], expl_momentum);

      precicec_writeAndMapData("Fluid-Mesh", "DragForce", 1, atom->x[i], f_drag[i]);
    }

    else if (drag_law == DRAG_STOKES)
    {
      // stokes drag force 
      for (int d = 0; d < 3; d++)
      {
        f_drag[i][d] = 3 * M_PI * mu_fluid * diameter * (1 - vol_frac[i]) * v_rel[d];
      }
    }
  }

  free(expl_momentum);
  free(impl_momentum);
  free(v_rel);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeFluidDragAtom::memory_usage()
{
  double bytes = nmax * size_peratom_cols * sizeof(double);
  return bytes;
}
