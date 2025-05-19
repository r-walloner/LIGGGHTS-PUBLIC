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
#include "compute_fluid_buoyancy_atom.h"
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeFluidBuoyancyAtom::ComputeFluidBuoyancyAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  peratom_flag = 1;
  size_peratom_cols = 3;

  c_p_grad_fluid = NULL;

  f_buoyancy = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeFluidBuoyancyAtom::~ComputeFluidBuoyancyAtom()
{
  memory->destroy(f_buoyancy);
}

/* ---------------------------------------------------------------------- */

void ComputeFluidBuoyancyAtom::init() {
  // Find the required computes
  c_p_grad_fluid = modify->find_compute_id("p_grad_fluid");
  if (!c_p_grad_fluid)
    error->all(FLERR, "Cannot find compute p_grad_gluid");
}

/* ---------------------------------------------------------------------- */

void ComputeFluidBuoyancyAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow data array if necessary
  if (atom->nlocal > nmax) {
    memory->destroy(f_buoyancy);
    nmax = atom->nlocal;
    memory->create(f_buoyancy, nmax, size_peratom_cols, "compute:fluid_buoyancy/atom");
    array_atom = f_buoyancy;
  }

  // Assume that all particles have the same radius
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++) {
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  }
  // Compute particle volume
  const double volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  // Access per-atom pressure gradient from compute
  if (!(c_p_grad_fluid->invoked_flag & INVOKED_PERATOM)) {
    c_p_grad_fluid->compute_peratom();
    c_p_grad_fluid->invoked_flag |= INVOKED_PERATOM;
  }
  double **p_grad_fluid = c_p_grad_fluid->array_atom;
  if (!p_grad_fluid)
    error->all(FLERR, "Compute does not provide a per-atom array");

  // Calculate and apply buoyancy force
  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit))
      continue;
    
    for (int d = 0; d < size_peratom_cols; d++)
      f_buoyancy[i][d] = -volume * p_grad_fluid[i][d];
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeFluidBuoyancyAtom::memory_usage()
{
  double bytes = nmax * size_peratom_cols * sizeof(double);
  return bytes;
}
