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

#include <string.h>
#include "precice/preciceC.h"
#include "compute_precice_read_atom.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "comm.h"
#include "fix_multisphere.h" 
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PreciceReadAtom::PreciceReadAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  // if (narg != iarg) error->all(FLERR,"Illegal compute ke/atom command");

  peratom_flag = 1;
  size_peratom_cols = 3; // TODO does this need to be configured from precice?
  
  nmax = 0;
  precice_data = NULL;
}

/* ---------------------------------------------------------------------- */

PreciceReadAtom::~PreciceReadAtom()
{
  memory->destroy(precice_data);
}

/* ---------------------------------------------------------------------- */

void PreciceReadAtom::init()
{
  // int count = 0;
  // for (int i = 0; i < modify->ncompute; i++)
  //   if (strcmp(modify->compute[i]->style,"ke/atom") == 0) count++;
  // if (count > 1 && comm->me == 0)
  //   error->warning(FLERR,"More than one compute ke/atom");
}

/* ---------------------------------------------------------------------- */

void PreciceReadAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow precice_data vector if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(precice_data);
    nmax = atom->nlocal;
    memory->create(precice_data,nmax,3,"precice");
    array_atom = precice_data;
  }

  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit))
      continue;
    
    precicec_mapAndReadData("Fluid-Mesh", "Force", 1, atom->x[i], update->dt, precice_data[i]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double PreciceReadAtom::memory_usage()
{
  double bytes = nmax * domain->dimension * sizeof(double);
  return bytes;
}
