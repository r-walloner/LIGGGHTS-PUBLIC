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

#include <stdlib.h>
#include <string.h>
#include "precice/preciceC.h"
#include "precice_initialize.h"
#include "domain.h"
#include "error.h"
#include "mpi_liggghts.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PreciceInitialize::PreciceInitialize(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void PreciceInitialize::command(int narg, char **arg)
{
  if (narg < 3)
    error->all(FLERR, "Illegal precice_initialize command");

  if (domain->box_exist == 0)
    error->all(FLERR, "Precice initialize command before simulation box is defined");

  // parse arguments
  const char *participant_name = arg[0];
  const char *config_file_name = arg[1];
  const int n_meshes = narg - 2;
  const char **mesh_names = (const char **)(arg + 2);
  
  // create participant
  int mpi_rank = -1;
  int mpi_size = -1;
  MPI_Comm_rank(world, &mpi_rank);
  MPI_Comm_size(world, &mpi_size);
  precicec_createParticipant_withCommunicator(
      participant_name, config_file_name, mpi_rank, mpi_size, &world);

  // set mesh access region
  const double mesh_access_region[6] = {
      domain->sublo[0], domain->subhi[0],
      domain->sublo[1], domain->subhi[1],
      domain->sublo[2], domain->subhi[2]};
  for (int i = 0; i < n_meshes; i++)
    precicec_setMeshAccessRegion(mesh_names[i], mesh_access_region);

  // initialize
  precicec_initialize();
}
