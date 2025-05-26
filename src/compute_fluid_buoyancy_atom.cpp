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
  size_peratom_cols = 0;

  f_buoyancy = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeFluidBuoyancyAtom::~ComputeFluidBuoyancyAtom()
{
  memory->destroy(f_buoyancy);
}

/* ---------------------------------------------------------------------- */

void ComputeFluidBuoyancyAtom::init() {}

/* ---------------------------------------------------------------------- */

void ComputeFluidBuoyancyAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow data array if necessary
  if (atom->nlocal > nmax) {
    memory->destroy(f_buoyancy);
    nmax = atom->nlocal;
    memory->create(f_buoyancy, nmax, "compute:fluid_buoyancy/atom");
    vector_atom = f_buoyancy;
  }

  // Assume that all particles have the same radius
  // TODO refactor this into its own compute to reduce code duplication
  //   - maybe we can use the compute property/atom command
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++) {
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  }
  // Compute particle volume
  const double volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  // TODO get constants from global properties
  const double rho_fluid = 1000; // kg/m^3
  const double g = 9.81; // m/s^2

  // Calculate and apply buoyancy force
  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit))
      continue;
    
    f_buoyancy[i] = volume * rho_fluid * g;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeFluidBuoyancyAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
