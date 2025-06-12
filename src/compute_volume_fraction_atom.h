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

#ifdef COMPUTE_CLASS

ComputeStyle(volume_fraction/atom,ComputeVolumeFractionAtom)

#else

#ifndef LMP_COMPUTE_VOLUME_FRACTION_ATOM_H
#define LMP_COMPUTE_VOLUME_FRACTION_ATOM_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeVolumeFractionAtom : public Compute {
 public:
  ComputeVolumeFractionAtom(class LAMMPS *, int &iarg, int, char **);
  ~ComputeVolumeFractionAtom();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double weightingFunction(const double r);
  double memory_usage();

 private:
  int method; // VOLUME_FRACTION_VORONOI or VOLUME_FRACTION_KERNEL
  Compute *c_voronoi;
  NeighList *list;
  double kernel_radius;
  double kernel_radius_sq;
  double *v_frac;
  int nmax;
};

}

#endif
#endif

/* ERROR/WARNING messages:

// TODO

*/
