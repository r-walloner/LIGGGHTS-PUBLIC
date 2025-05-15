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

#include <stdlib.h>
#include <string.h>
#include "precice/preciceC.h"
#include "fix_map_and_read_data.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ----------------------------------------------------------------------
   parse arguments
------------------------------------------------------------------------- */

FixMapAndReadData::FixMapAndReadData(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 7)
    error->all(FLERR, "Illegal fix map_and_read_data command");

  // parse arguments and copy strings
  mesh_name = strdup(arg[3]);
  data_name = strdup(arg[4]);
  relative_read_time = force->numeric(FLERR, arg[5]);
  target_property = strdup(arg[6]);
}

/* ----------------------------------------------------------------------
   free allocated memory
------------------------------------------------------------------------- */

FixMapAndReadData::~FixMapAndReadData()
{
  free(mesh_name);
  free(data_name);
  free(target_property);
}

/* ----------------------------------------------------------------------
   set mask for when to call this fix
-------------------------------------------------------------------------- */

int FixMapAndReadData::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixMapAndReadData::init()
{
  // error checks on coarsegraining
  // TODO do we need this?
  if (force->cg_active())
    error->cg(FLERR, this->style);
}

/* ---------------------------------------------------------------------- */

void FixMapAndReadData::setup(int vflag)
{
  if (!strstr(update->integrate_style, "verlet"))
    error->all(FLERR, "fix map_and_read_data is only implemented for verlet integration.");

  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMapAndReadData::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixMapAndReadData::initial_integrate(int vflag)
{
  if (precicec_requiresWritingCheckpoint())
  {
    // TODO write checkpoint
  }
}

/* ---------------------------------------------------------------------- */

void FixMapAndReadData::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  // constant force
  // potential energy = - x dot f in unwrapped coords

  double read_forces[3];
  for (int i = 0; i < nlocal; i++)
  {
    if (!(mask[i] & groupbit))
      continue;

    // TODO relative read time
    relative_read_time = update->dt;
    precicec_mapAndReadData(mesh_name, data_name, 1, x[i], relative_read_time, read_forces);

    f[i][0] += read_forces[0];
    f[i][1] += read_forces[1];
    f[i][2] += read_forces[2];
  }
}

/* ---------------------------------------------------------------------- */

void FixMapAndReadData::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */
void FixMapAndReadData::end_of_step()
{
}
