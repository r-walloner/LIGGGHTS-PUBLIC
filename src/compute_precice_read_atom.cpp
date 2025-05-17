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
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PreciceReadAtom::PreciceReadAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  if ((narg - iarg < 2) || (narg - iarg > 3))
    error->all(FLERR,"Illegal compute precice_read/atom command");

  // parse arguments

  mesh_name = strdup(arg[iarg++]);
  data_name = strdup(arg[iarg++]);

  if (iarg < narg)
    relative_read_time = force->numeric(FLERR,arg[iarg++]);
  else
    relative_read_time = 1.0;
  relative_read_time *= update->dt;

  // set  up compute

  peratom_flag = 1;
  size_peratom_cols = precicec_getDataDimensions(mesh_name, data_name);

  nmax = 0;
  data = NULL;
}

/* ---------------------------------------------------------------------- */

PreciceReadAtom::~PreciceReadAtom()
{
  free(mesh_name);
  free(data_name);
  memory->destroy(data);
}

/* ---------------------------------------------------------------------- */

void PreciceReadAtom::init() {}

/* ---------------------------------------------------------------------- */

void PreciceReadAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow data array if necessary

  if (atom->nlocal > nmax) {
    memory->destroy(data);
    nmax = atom->nlocal;
    char mem_name[256];
    snprintf(mem_name, sizeof(mem_name), "precice%s%s", mesh_name, data_name);
    memory->create(data, nmax, size_peratom_cols, mem_name);
    array_atom = data;
  }

  // read data

  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit))
      continue;
    
    precicec_mapAndReadData(
      mesh_name, data_name, 1, atom->x[i], relative_read_time, data[i]);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double PreciceReadAtom::memory_usage()
{
  double bytes = nmax * size_peratom_cols * sizeof(double);
  return bytes;
}
