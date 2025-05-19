/* ----------------------------------------------------------------------
    This is the

    ██╗     ██╗ ██████╗  ██████╗  ██████╗ ██╗  ██╗████████╗███████╗
    ██║     ██║██╔════╝ ██╔════╝ ██╔════╝ ██║  ██║╚══██╔══╝██╔════╝
    ██║     ██║██║  ███╗██║  ███╗██║  ███╗███████║   ██║   ███████╗
    ██║     ██║██║   ██║██║   ██║██║   ██║██╔══██║   ██║   ╚════██║
    ███████╗██║╚██████╔╝╚██████╔╝╚██████╔╝██║  ██║   ██║   ███████║
    ╚══════╝╚═╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚══════╝®

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
// TODO clean headers
#include <cassert>
#include <cmath>
#include "fix_fluid_drag.h"
#include "atom.h"
#include "error.h"
#include "math_extra.h"
#include "modify.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathExtra;

/* ----------------------------------------------------------------------
   parse arguments
------------------------------------------------------------------------- */

FixFluidDrag::FixFluidDrag(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  c_v_fluid = NULL;
  c_vol_frac = NULL;
}

/* ----------------------------------------------------------------------
   free allocated memory
------------------------------------------------------------------------- */

FixFluidDrag::~FixFluidDrag()
{
}

/* ----------------------------------------------------------------------
   set mask for when to call this fix
-------------------------------------------------------------------------- */

int FixFluidDrag::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::init()
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

void FixFluidDrag::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::post_force(int vflag)
{
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
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++)
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  // Compute particle diameter and volume
  const double diameter = 2 * radius;
  const double volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  // TODO get constants from global properties
  const double rho_fluid = 1000.0; // kg/m^3
  const double mu_fluid = 0.001;   // Pa.s

  double beta, drag_coeff, reynolds;
  double *v_rel = (double *)malloc(3 * sizeof(double));

  // Calculate and apply drag force
  for (int i = 0; i < atom->nlocal; i++)
  {
    if (!(atom->mask[i] & groupbit))
      continue;

    // particle velocity relative to fluid
    sub3(v_fluid[i], atom->v[i], v_rel);
    const double mag_v_rel = len3(v_rel);

    // reynolds number
    reynolds = (1 - vol_frac[i]) * rho_fluid * mag_v_rel * diameter /
               mu_fluid;

    // drag coefficient
    if (reynolds < 1000)
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

    // drag force
    for (int d = 0; d < 3; d++)
      atom->f[i][d] += beta * volume * v_rel[d] / vol_frac[i];
  }

  free(v_rel);
}

/* ---------------------------------------------------------------------- */

void FixFluidDrag::min_post_force(int vflag)
{
  post_force(vflag);
}
