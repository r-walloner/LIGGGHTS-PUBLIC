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
#include "compute_volume_fraction_atom.h"
#include "atom.h"
#include "error.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeVolumeFractionAtom::ComputeVolumeFractionAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) :
  Compute(lmp, iarg, narg, arg)
{
  peratom_flag = 1;
  size_peratom_cols = 0;

  c_voronoi = NULL;

  v_frac = NULL;
  nmax = 0;
}

/* ---------------------------------------------------------------------- */

ComputeVolumeFractionAtom::~ComputeVolumeFractionAtom()
{
  memory->destroy(v_frac);
}

/* ---------------------------------------------------------------------- */

void ComputeVolumeFractionAtom::init() {
  // Find the ID of the voronoi tessellation compute
  c_voronoi = modify->find_compute_style_strict("voronoi/atom", 0);
  if (!c_voronoi)
    error->all(FLERR, "Cannot find voronoi/atom compute");
}

/* ---------------------------------------------------------------------- */

void ComputeVolumeFractionAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // Assume that all particles have the same radius
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++) {
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  }
  // Compute particle volume
  const double particle_volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  // Access per-atom voronoi cell volume from compute
  if (!(c_voronoi->invoked_flag & INVOKED_PERATOM)) {
    c_voronoi->compute_peratom();
    c_voronoi->invoked_flag |= INVOKED_PERATOM;
  }
  double **voro_tess = c_voronoi->array_atom;
  if (!voro_tess)
    error->all(FLERR, "voronoi/atom compute does not provide a per-atom array");

  // grow data array if necessary
  if (atom->nlocal > nmax) {
    memory->destroy(v_frac);
    nmax = atom->nlocal;
    memory->create(v_frac, nmax, "compute:volume_fraction/atom");
    vector_atom = v_frac;
  }

  // Calculate volume fraction
  for (int i = 0; i < atom->nlocal; i++) {
    if (!(atom->mask[i] & groupbit)) {
      v_frac[i] = 0.0;
      continue;
    }

    // Volume fraction = particle volume / voronoi cell volume
    v_frac[i] = particle_volume / voro_tess[i][0];
    
    assert(v_frac[i] >= 0.0 && "Particle volume fraction is negative");
    assert(v_frac[i] <= 1.0 && "Particle volume fraction is greater than 1.0");

    // Write volumes to precice
    // TODO move this to a separate fix or compute
    precicec_writeAndMapData("Fluid-Mesh", "Alpha", 1, atom->x[i], &particle_volume);
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeVolumeFractionAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
