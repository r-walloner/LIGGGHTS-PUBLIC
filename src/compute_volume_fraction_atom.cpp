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
#include <string.h>
#include "compute_volume_fraction_atom.h"
#include "atom.h"
#include "error.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "update.h"
#include "vector_liggghts.h"

using namespace LAMMPS_NS;

enum
{
  VOLUME_FRACTION_VORONOI,
  VOLUME_FRACTION_KERNEL
};

/* ---------------------------------------------------------------------- */

ComputeVolumeFractionAtom::ComputeVolumeFractionAtom(LAMMPS *lmp, int &iarg, int narg, char **arg) : Compute(lmp, iarg, narg, arg)
{
  if (narg < iarg + 1)
    error->all(FLERR, "Illegal compute volume_fraction/atom command");

  method = -1;
  if (strcmp(arg[iarg], "voronoi") == 0)
    method = VOLUME_FRACTION_VORONOI;
  else if (strcmp(arg[iarg], "kernel") == 0)
    method = VOLUME_FRACTION_KERNEL;
  if (method == -1)
    error->all(FLERR, "Illegal compute volume_fraction/atom command");
  iarg++;

  if (method == VOLUME_FRACTION_KERNEL)
  {
    if (narg < iarg + 1)
      error->all(FLERR, "Illegal compute volume_fraction/atom command");
    kernel_radius = force->numeric(FLERR, arg[iarg++]);
    kernel_radius_sq = kernel_radius * kernel_radius;
    if (kernel_radius <= 0.0)
      error->all(FLERR, "kernel_radius > 0 required");
  }

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

void ComputeVolumeFractionAtom::init()
{
  if (method == VOLUME_FRACTION_VORONOI)
  {
    // Find the ID of the voronoi tessellation compute
    c_voronoi = modify->find_compute_style_strict("voronoi/atom", 0);
    if (!c_voronoi)
      error->all(FLERR, "Cannot find voronoi/atom compute");
  }
}

/* ---------------------------------------------------------------------- */

void ComputeVolumeFractionAtom::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */
// Gaussian weighting function for kernel method
double ComputeVolumeFractionAtom::weightingFunction(const double r)
{
  // w = exp(-r^2/kernel_sqRadius_*9)
  // the constant 9 is coming from 3^2 which is the cutoff of the gaussian (at 3 sigma)
  const double w = exp(-r * r * 9.0 / kernel_radius_sq);
  // phi = w(r)/[4 pi \int_0^inf w(s) s^2 ds]
  // constant below is (pi (e^9 sqrt(pi) erf(3)-6))/(27 e^9)
  const double normalization = kernel_radius_sq * kernel_radius_sq * 0.20614365813686698838756664006357;
  return w / normalization;
}

/* ---------------------------------------------------------------------- */

void ComputeVolumeFractionAtom::compute_peratom()
{
  invoked_peratom = update->ntimestep;

  // grow data array if necessary
  if (atom->nlocal > nmax)
  {
    memory->destroy(v_frac);
    nmax = atom->nlocal;
    memory->create(v_frac, nmax, "compute:volume_fraction/atom");
    vector_atom = v_frac;
  }

  // Assume that all particles have the same radius
  const double radius = atom->radius[0];
  for (int i = 1; i < atom->nlocal; i++)
  {
    assert(atom->radius[i] == radius && "All particles must have the same radius");
  }
  // Compute particle volume
  const double particle_volume = (4.0 / 3.0) * M_PI * pow(radius, 3);

  if (method == VOLUME_FRACTION_VORONOI)
  {
    // Access per-atom voronoi cell volume from compute
    if (!(c_voronoi->invoked_flag & INVOKED_PERATOM))
    {
      c_voronoi->compute_peratom();
      c_voronoi->invoked_flag |= INVOKED_PERATOM;
    }
    double **voro_tess = c_voronoi->array_atom;
    if (!voro_tess)
      error->all(FLERR, "voronoi/atom compute does not provide a per-atom array");

    // Calculate volume fraction
    for (int i = 0; i < atom->nlocal; i++)
    {
      if (!(atom->mask[i] & groupbit))
      {
        v_frac[i] = 0.0;
        continue;
      }

      // Volume fraction = particle volume / voronoi cell volume
      v_frac[i] = particle_volume / voro_tess[i][0];

      assert(v_frac[i] >= 0.0 && "Particle volume fraction is negative");
      assert(v_frac[i] <= 1.0 && "Particle volume fraction is greater than 1.0");
    }
  }

  // TODO does not work jet
  else if (method == VOLUME_FRACTION_KERNEL)
  {
    for (int i = 0; i < atom->nlocal; i++)
    {
      if (!(atom->mask[i] & groupbit))
      {
        v_frac[i] = 0.0;
        continue;
      }

      v_frac[i] = .1;
    }
  }
  // {
  //   int i, j;

  //   for (int ii = 0; ii < list->inum; ii++)
  //   {
  //     i = list->ilist[ii];
  //     if (!(atom->mask[i] & groupbit))
  //     {
  //       v_frac[i] = 0.0;
  //       continue;
  //     }

  //     for (int jj = 0; jj <= list->numneigh[i]; jj++)
  //     {
  //       if (jj == list->numneigh[i])
  //         j = i; // last j is i itself
  //       else
  //         j = list->firstneigh[i][jj];

  //       if (!(atom->mask[j] & groupbit))
  //         continue;

  //       // compute particle distance
  //       const double xij[3] = {atom->x[i][0] - atom->x[j][0],
  //                              atom->x[i][1] - atom->x[j][1],
  //                              atom->x[i][2] - atom->x[j][2]};
  //       const double dist_ij_sq = vectorMag3DSquared(xij);
  //       const double dist_ij = sqrt(dist_ij_sq);

  //       if (dist_ij_sq > kernel_radius_sq)
  //         continue; // skip if outside kernel radius

  //       v_frac[i] += weightingFunction(dist_ij) * particle_volume;
  //     }
  //   }
  // }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeVolumeFractionAtom::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
