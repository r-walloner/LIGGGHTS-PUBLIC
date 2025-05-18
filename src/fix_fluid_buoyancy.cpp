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
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include "precice/preciceC.h"
#include "fix_fluid_buoyancy.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"
#include "compute.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   parse arguments
------------------------------------------------------------------------- */

FixFluidBuoyancy::FixFluidBuoyancy(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 3)
    error->all(FLERR, "Illegal fix fluid/buoyancy command");

  
}

/* ----------------------------------------------------------------------
   free allocated memory
------------------------------------------------------------------------- */

FixFluidBuoyancy::~FixFluidBuoyancy()
{
}

/* ----------------------------------------------------------------------
   set mask for when to call this fix
-------------------------------------------------------------------------- */

int FixFluidBuoyancy::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixFluidBuoyancy::init()
{
  // Store the compute index for later use
  icompute_p_grad_fluid = modify->find_compute("p_grad_fluid");
  if (icompute_p_grad_fluid < 0)
    error->all(FLERR, "Compute ID 'p_grad_fluid' not found");
}

/* ---------------------------------------------------------------------- */

void FixFluidBuoyancy::setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluidBuoyancy::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixFluidBuoyancy::post_force(int vflag)
{
  // Assume that all particles have the same radius
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++) {
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  }
  // Compute particle volume
  const double volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  // Access per-atom pressure gradient from compute
  Compute *c = modify->compute[icompute_p_grad_fluid];
  if (!(c->invoked_flag & INVOKED_PERATOM)) {
    c->compute_peratom();
    c->invoked_flag |= INVOKED_PERATOM;
  }
  double **p_grad_fluid = c->array_atom;
  if (!p_grad_fluid)
    error->all(FLERR, "Compute does not provide a per-atom array");

  // Calculate and apply buoyancy force
  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit))
      continue;
    
    for (int d = 0; d < 3; d++)
      atom->f[i][d] += -volume * p_grad_fluid[i][d];
  }
}

/* ---------------------------------------------------------------------- */

void FixFluidBuoyancy::min_post_force(int vflag)
{
  post_force(vflag);
}
