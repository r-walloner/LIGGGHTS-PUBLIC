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
#include "precice_advance.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "memory.h"
#include "modify.h"
#include "neighbor.h"
#include "read_restart.h"
#include "signal_handling.h"
#include "universe.h"
#include "update.h"
#include "write_restart.h"

using namespace LAMMPS_NS;

// same as write_restart.cpp

enum{VERSION,SMALLINT,TAGINT,BIGINT,
       UNITS,NTIMESTEP,DIMENSION,NPROCS,PROCGRID_0,PROCGRID_1,PROCGRID_2,
       NEWTON_PAIR,NEWTON_BOND,XPERIODIC,YPERIODIC,ZPERIODIC,
       BOUNDARY_00,BOUNDARY_01,BOUNDARY_10,BOUNDARY_11,BOUNDARY_20,BOUNDARY_21,
       ATOM_STYLE,NATOMS,NTYPES,
       NBONDS,NBONDTYPES,BOND_PER_ATOM,
       NANGLES,NANGLETYPES,ANGLE_PER_ATOM,
       NDIHEDRALS,NDIHEDRALTYPES,DIHEDRAL_PER_ATOM,
       NIMPROPERS,NIMPROPERTYPES,IMPROPER_PER_ATOM,
       BOXLO_0,BOXHI_0,BOXLO_1,BOXHI_1,BOXLO_2,BOXHI_2,
       SPECIAL_LJ_1,SPECIAL_LJ_2,SPECIAL_LJ_3,
       SPECIAL_COUL_1,SPECIAL_COUL_2,SPECIAL_COUL_3,
       XY,XZ,YZ};
enum{MASS};
enum{PAIR,BOND,ANGLE,DIHEDRAL,IMPROPER};

/* ---------------------------------------------------------------------- */

PreciceAdvance::PreciceAdvance(LAMMPS *lmp) :
    ReadRestart(lmp),
    checkpoint_file(strdup("precice_checkpoint.%"))
{}

/* ---------------------------------------------------------------------- */

PreciceAdvance::~PreciceAdvance()
{
  free(checkpoint_file);
}

/* ---------------------------------------------------------------------- */

void PreciceAdvance::command(int narg, char **arg)
{
  precicec_advance(update->dt);

  // // revert to last checkpoint if needed
  // if (precicec_requiresReadingCheckpoint())
  //   read_checkpoint();

  // end simulation if coupling is finished
  if (!precicec_isCouplingOngoing()) {
    precicec_finalize();
    SignalHandler::request_quit = true;
    return;
  }

  // // write checkpoint if needed
  // if (precicec_requiresWritingCheckpoint())
  //   write_checkpoint();
}

/* ---------------------------------------------------------------------- */

// void PreciceAdvance::write_checkpoint()
// {
//   WriteRestart *restart = new WriteRestart(lmp);
//   restart->write(checkpoint_file);
//   delete restart;
// }

/* ---------------------------------------------------------------------- */

// void PreciceAdvance::read_checkpoint()
// {
//   MPI_Comm_rank(world, &me);
//   MPI_Comm_size(world, &nprocs);

//   // open base file
//   if (me == 0) {
//     if (screen) fprintf(screen,"Reading checkpoint file ...\n");
//     char *hfile;
//     hfile = new char[strlen(checkpoint_file) + 16];
//     char *ptr = strchr(checkpoint_file,'%');
//     *ptr = '\0';
//     sprintf(hfile,"%s%s%s",checkpoint_file,"base",ptr+1);
//     *ptr = '%';
//     fp = fopen(hfile,"rb");
//     if (fp == NULL) {
//       char str[512];
//       sprintf(str,"Cannot open checkpoint file %s",hfile);
//       error->one(FLERR,str);
//     }
//     delete [] hfile;
//   }

//   // read header info and create atom style and simulation box

//   read_header();

//   // Create a new AtomVec instance with the specified style
//   const char *style = strdup(atom->atom_style);
//   atom->create_avec(style, 0, nullptr, nullptr); // TODO dont pass nullptrs
//   atom->nlocal = 0;
//   AtomVec *avec = atom->avec;

//   // problem setup using info from header
//   // TODO do we need this?
//   int n;
//   if (nprocs == 1) n = static_cast<int> (atom->natoms);
//   else n = static_cast<int> (1.1 * atom->natoms / nprocs); // TODO check LB_FACTOR

//   atom->allocate_type_arrays();
//   atom->avec->grow(n);
//   n = atom->nmax;

//   // read group names
//   group->read_restart(fp);

//   // read atom type arrays
//   type_arrays();

//   // read in force fields
//   // TODO we probably dont need this. However, other usecases might.
//   force_fields();

//   int nextra = modify->read_restart(fp);
//   // TODO we probably dont need this. However, other usecases might.
//   atom->nextra_store = nextra;
//   memory->create(atom->extra,atom->nmax,nextra,"atom:extra");

//   // read atom data
//   // each proc reads its own file, keeping all atoms in the files
//   // perform irregular comm to migrate atoms to correct procs
//   // close restart file when done

//   int maxbuf = 0;
//   double *buf = NULL;
  
//   if (me == 0) fclose(fp);
//   char *perproc = new char[strlen(checkpoint_file) + 16];
//   char *ptr = strchr(checkpoint_file,'%');
  
//   *ptr = '\0';
//   sprintf(perproc,"%s%d%s",checkpoint_file,me,ptr+1);
//   *ptr = '%';
//   fp = fopen(perproc,"rb");
//   if (fp == NULL) {
//     char str[512];
//     sprintf(str,"Cannot open checkpoint file %s",perproc);
//     error->one(FLERR,str);
//   }
  
//   nread_int(&n,1,fp);
//   if (n > maxbuf) {
//     maxbuf = n;
//     memory->destroy(buf);
//     memory->create(buf,maxbuf,"read_restart:buf");
//   }
//   if (n > 0) nread_double(buf,n,fp);

//   int m = 0;
//   while (m < n)
//     m += avec->unpack_restart(&buf[m]);
//   fclose(fp);

//   delete [] perproc;

//   // // create a temporary fix to hold and migrate extra atom info
//   // // necessary b/c irregular will migrate atoms
//   // // TODO do we need this? 

//   // if (nextra) {
//   //   char cextra[8],fixextra[8];
//   //   sprintf(cextra,"%d",nextra);
//   //   sprintf(fixextra,"%d",modify->nfix_restart_peratom);
//   //   char **newarg = new char*[5];
//   //   newarg[0] = (char *) "_read_restart";
//   //   newarg[1] = (char *) "all";
//   //   newarg[2] = (char *) "READ_RESTART";
//   //   newarg[3] = cextra;
//   //   newarg[4] = fixextra;
//   //   modify->add_fix(5,newarg);
//   //   delete [] newarg;
//   // }

//   // // move atoms to new processors via irregular()
//   // // in case read by different proc than wrote restart file
//   // // first do map_init() since irregular->migrate_atoms() will do map_clear()

//   // if (atom->map_style) atom->map_init();
//   // if (domain->triclinic) domain->x2lamda(atom->nlocal);
//   // Irregular *irregular = new Irregular(lmp);
//   // irregular->migrate_atoms();
//   // delete irregular;
//   // if (domain->triclinic) domain->lamda2x(atom->nlocal);

//   // // put extra atom info held by fix back into atom->extra
//   // // destroy temporary fix

//   // if (nextra) {
//   //   memory->destroy(atom->extra);
//   //   memory->create(atom->extra,atom->nmax,nextra,"atom:extra");
//   //   int ifix = modify->find_fix("_read_restart");
//   //   FixReadRestart *fix = (FixReadRestart *) modify->fix[ifix];
//   //   int *count = fix->count;
//   //   double **extra = fix->extra;
//   //   double **atom_extra = atom->extra;
//   //   int nlocal = atom->nlocal;
//   //   for (int i = 0; i < nlocal; i++)
//   //     for (int j = 0; j < count[i]; j++)
//   //       atom_extra[i][j] = extra[i][j];
//   //   modify->delete_fix("_read_restart");
//   // }

//   // clean-up memory

//   memory->destroy(buf);

//   // check that all atoms were assigned to procs

//   bigint natoms;
//   bigint nblocal = atom->nlocal;
//   MPI_Allreduce(&nblocal,&natoms,1,MPI_LMP_BIGINT,MPI_SUM,world);

//   if (me == 0) {
//     if (screen) fprintf(screen,"  " BIGINT_FORMAT " atoms\n",natoms);
//     if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " atoms\n",natoms);
//   }

//   if (natoms != atom->natoms)
//     error->all(FLERR,"Did not assign all atoms correctly");

//   if (me == 0) {
//     if (atom->nbonds) {
//       if (screen) fprintf(screen,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
//       if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " bonds\n",atom->nbonds);
//     }
//     if (atom->nangles) {
//       if (screen) fprintf(screen,"  " BIGINT_FORMAT " angles\n",
//                           atom->nangles);
//       if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " angles\n",
//                            atom->nangles);
//     }
//     if (atom->ndihedrals) {
//       if (screen) fprintf(screen,"  " BIGINT_FORMAT " dihedrals\n",
//                           atom->ndihedrals);
//       if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " dihedrals\n",
//                            atom->ndihedrals);
//     }
//     if (atom->nimpropers) {
//       if (screen) fprintf(screen,"  " BIGINT_FORMAT " impropers\n",
//                           atom->nimpropers);
//       if (logfile) fprintf(logfile,"  " BIGINT_FORMAT " impropers\n",
//                            atom->nimpropers);
//     }
//   }

//   // // check if tags are being used
//   // // create global mapping and bond topology now that system is defined
//   // // TODO We probably dont need this. However, other usecases might.
//   // int flag = 0;
//   // for (int i = 0; i < atom->nlocal; i++)
//   //   if (atom->tag[i] > 0) flag = 1;
//   // int flag_all;
//   // MPI_Allreduce(&flag,&flag_all,1,MPI_INT,MPI_MAX,world);
//   // if (atom->natoms > 0 && flag_all == 0) atom->tag_enable = 0;

//   // if (atom->map_style) {
//   //   atom->map_init();
//   //   atom->map_set();
//   // }
//   // if (atom->molecular) {
//   //   Special special(lmp);
//   //   special.build();
//   // }

// }

/* ----------------------------------------------------------------------
   read header of checkpoint file
------------------------------------------------------------------------- */

// void PreciceAdvance::read_header()
// {
//   int px=0,py=0,pz=0;
//   int xperiodic=0,yperiodic=0,zperiodic=0;
//   int boundary[3][2] = {};
//   double boxlo[3],boxhi[3];
//   double special_lj[4],special_coul[4];
//   double triclinic_xy=NULL,triclinic_xz=NULL,triclinic_yz=NULL;

//   // read flags and values until flag = -1

//   int flag = read_int();
//   while (flag >= 0) {

//     // check checkpoint file version, error if different

//     if (flag == VERSION) {
//       char *version = read_char();
//       if (strcmp(version,universe->version) != 0 && me == 0) {
//         error->all(FLERR,
//                    "Checkpoint file version does not match LIGGGHTS version"); 
//         if (screen) fprintf(screen,"   --> restart file = %s\n   --> LIGGGHTS = %s\n", 
//                             version,universe->version);
//       }
//       // parse version number
//       // version format is:
//       // Version LIGGGHTS-REPOSITORY-NAME MAJOR.MINOR.[....]
//       // MAJOR and MINOR are integers
//       std::string ver = std::string(version);
//       std::size_t space1 = ver.find(' ');
//       std::size_t space2 = ver.find(' ', space1+1);
//       std::size_t dot1 = ver.find('.', space2+1);
//       std::size_t dot2 = ver.find('.', dot1+1);
//       if (space1 != std::string::npos &&
//           space2 != std::string::npos &&
//           dot1 != std::string::npos &&
//           dot2 != std::string::npos)
//       {
//           std::string ver_major = ver.substr(space2+1, dot1-space2-1);
//           std::string ver_minor = ver.substr(dot1+1, dot2-dot1-1);
//           restart_major = atoi(ver_major.c_str());
//           restart_minor = atoi(ver_minor.c_str());
//           printf("version %d %d\n", restart_major, restart_minor);
//       }
//       delete [] version;

//     // check lmptype.h sizes, error if different

//     } else if (flag == SMALLINT) {
//       int size = read_int();
//       if (size != sizeof(smallint))
//         error->all(FLERR,"Smallint setting in lmptype.h is not compatible");
//     } else if (flag == TAGINT) {
//       int size = read_int();
//       if (size != sizeof(tagint))
//         error->all(FLERR,"Tagint setting in lmptype.h is not compatible");
//     } else if (flag == BIGINT) {
//       int size = read_int();
//       if (size != sizeof(bigint))
//         error->all(FLERR,"Bigint setting in lmptype.h is not compatible");

//     // read unit_style from checkpoint file, error if different
      
//      } else if (flag == UNITS) {
//       char *style = read_char();
//       if (strcmp(style,update->unit_style) != 0 && me == 0)
//         error->all(FLERR,"Checkpoint file used different units");
//       delete [] style;

//     // set timestep from checkpoint file

//     } else if (flag == NTIMESTEP) {
//       update->ntimestep = read_bigint();

//     // read dimension from checkpoint file, error if different

//     } else if (flag == DIMENSION) {
//       int dimension = read_int();
//       if (dimension != domain->dimension && me == 0)
//         error->all(FLERR,"Checkpoint file used different dimension");

//     // read nprocs from checkpoint file, error if different

//     } else if (flag == NPROCS) {
//       nprocs_file = read_int();
//       if (nprocs_file != comm->nprocs && me == 0)
//         error->all(FLERR,"Checkpoint file used different # of processors");

//     // read procgrid from checkpoint file, error if different

//     } else if (flag == PROCGRID_0) {
//       px = read_int();
//     } else if (flag == PROCGRID_1) {
//       py = read_int();
//     } else if (flag == PROCGRID_2) {
//       pz = read_int();
//       if (comm->user_procgrid[0] != 0 &&
//           (px != comm->user_procgrid[0] || py != comm->user_procgrid[1] ||
//            pz != comm->user_procgrid[2]) && me == 0)
//         error->all(FLERR,"Checkpoint file used different 3d processor grid");

//     // read newton_pair from checkpoint file, error if different

//     } else if (flag == NEWTON_PAIR) {
//       int newton_pair_file = read_int();
//       if (newton_pair_file != force->newton_pair && me == 0)
//         error->all(FLERR,"Checkpoint file used different newton pair setting");
    
//     // read newton_bond from checkpoint file, error if different

//     } else if (flag == NEWTON_BOND) {
//       int newton_bond_file = read_int();
//       if (newton_bond_file != force->newton_bond && me == 0)
//         error->all(FLERR,"Checkpoint file used different newton bond setting");

//     // read boundary settings from checkpoint file, error if different

//     } else if (flag == XPERIODIC) {
//       xperiodic = read_int();
//     } else if (flag == YPERIODIC) {
//       yperiodic = read_int();
//     } else if (flag == ZPERIODIC) {
//       zperiodic = read_int();
//     } else if (flag == BOUNDARY_00) {
//       boundary[0][0] = read_int();
//     } else if (flag == BOUNDARY_01) {
//       boundary[0][1] = read_int();
//     } else if (flag == BOUNDARY_10) {
//       boundary[1][0] = read_int();
//     } else if (flag == BOUNDARY_11) {
//       boundary[1][1] = read_int();
//     } else if (flag == BOUNDARY_20) {
//       boundary[2][0] = read_int();
//     } else if (flag == BOUNDARY_21) {
//       boundary[2][1] = read_int();

//       if (xperiodic != domain->periodicity[0] ||
//           yperiodic != domain->periodicity[1] ||
//           zperiodic != domain->periodicity[2]) {
//         if (me == 0)
//           error->all(FLERR,"Checkpoint file used different domain periodicity");
//       }

//       if (boundary[0][0] != domain->boundary[0][0] ||
//           boundary[0][1] != domain->boundary[0][1] ||
//           boundary[1][0] != domain->boundary[1][0] ||
//           boundary[1][1] != domain->boundary[1][1] ||
//           boundary[2][0] != domain->boundary[2][0] ||
//           boundary[2][1] != domain->boundary[2][1]) {
//         if (me == 0)
//           error->all(FLERR,"Checkpoint file used different boundary settings");
//       }

//     // read atom_style from checkpoint file, error if different

//     } else if (flag == ATOM_STYLE) {
//       char *style = read_char();

//       int nwords = 0;
//       char **words = NULL;

//       if (strcmp(style,"hybrid") == 0) {
//         nwords = read_int();
//         words = new char*[nwords];
//         for (int i = 0; i < nwords; i++) words[i] = read_char();
//       }

//       if (strcmp(style,atom->atom_style) != 0 && me == 0) {
//         error->all(FLERR,"Checkpoint file used different atom_style");
//       }

//       // TODO do we need to do something different here?
//       // atom->create_avec(style,nwords,words);
//       // atom->avec->read_restart_settings(fp);

//       for (int i = 0; i < nwords; i++) delete [] words[i];
//       delete [] words;
//       delete [] style;

//     // read counts from checkpoint file, warn/error if different
    
//     } else if (flag == NATOMS) {
//       atom->natoms = read_bigint();
//     } else if (flag == NTYPES) {
//       atom->ntypes = read_int();
//     } else if (flag == NBONDS) {
//       atom->nbonds = read_bigint();
//     } else if (flag == NBONDTYPES) {
//       atom->nbondtypes = read_int();
//     } else if (flag == BOND_PER_ATOM) {
//       atom->bond_per_atom = read_int();
//     } else if (flag == NANGLES) {
//       atom->nangles = read_bigint();
//     } else if (flag == NANGLETYPES) {
//       atom->nangletypes = read_int();
//     } else if (flag == ANGLE_PER_ATOM) {
//       atom->angle_per_atom = read_int();
//     } else if (flag == NDIHEDRALS) {
//       atom->ndihedrals = read_bigint();
//     } else if (flag == NDIHEDRALTYPES) {
//       atom->ndihedraltypes = read_int();
//     } else if (flag == DIHEDRAL_PER_ATOM) {
//       atom->dihedral_per_atom = read_int();
//     } else if (flag == NIMPROPERS) {
//       atom->nimpropers = read_bigint();
//     } else if (flag == NIMPROPERTYPES) {
//       atom->nimpropertypes = read_int();
//     } else if (flag == IMPROPER_PER_ATOM) {
//       atom->improper_per_atom = read_int();

//     // read simulation box from checkpoint file, warn if different

//     } else if (flag == BOXLO_0) {
//       boxlo[0] = read_double();
//     } else if (flag == BOXHI_0) {
//       boxhi[0] = read_double();
//     } else if (flag == BOXLO_1) {
//       boxlo[1] = read_double();
//     } else if (flag == BOXHI_1) {
//       boxhi[1] = read_double();
//     } else if (flag == BOXLO_2) {
//       boxlo[2] = read_double();
//     } else if (flag == BOXHI_2) {
//       boxhi[2] = read_double();

//       if (boxlo[0] != domain->boxlo[0] || boxhi[0] != domain->boxhi[0] ||
//           boxlo[1] != domain->boxlo[1] || boxhi[1] != domain->boxhi[1] ||
//           boxlo[2] != domain->boxlo[2] || boxhi[2] != domain->boxhi[2]) {
//         if (me == 0)
//           error->warning(FLERR,"Checkpoint file has different simulation box");
//       }

//     // read special lj/coul prefactors from checkpoint file, warn if different

//     } else if (flag == SPECIAL_LJ_1) {
//       special_lj[1] = read_double();
//     } else if (flag == SPECIAL_LJ_2) {
//       special_lj[2] = read_double();
//     } else if (flag == SPECIAL_LJ_3) {
//       special_lj[3] = read_double();
//     } else if (flag == SPECIAL_COUL_1) {
//       special_coul[1] = read_double();
//     } else if (flag == SPECIAL_COUL_2) {
//       special_coul[2] = read_double();
//     } else if (flag == SPECIAL_COUL_3) {
//       special_coul[3] = read_double();

//       if (special_lj[1] != force->special_lj[1] ||
//           special_lj[2] != force->special_lj[2] ||
//           special_lj[3] != force->special_lj[3] ||
//           special_coul[1] != force->special_coul[1] ||
//           special_coul[2] != force->special_coul[2] ||
//           special_coul[3] != force->special_coul[3]) {
//         if (me == 0)
//           error->warning(FLERR,"Checkpoint file has different force prefactors");
//       }

//     // read triclinic box angles from checkpoint file, warn if different

//     } else if (flag == XY) {
//       triclinic_xy = read_double();
//     } else if (flag == XZ) {
//       triclinic_xz = read_double();
//     } else if (flag == YZ) {
//       triclinic_yz = read_double();

//       if (triclinic_xy != NULL || triclinic_xz != NULL || triclinic_yz != NULL) {
//         if (domain->triclinic != 1 && me == 0)
//           error->all(FLERR,"Checkpoint file has triclinic box");

//         if (triclinic_xy != domain->xy || triclinic_xz != domain->xz ||
//             triclinic_yz != domain->yz) {
//           if (me == 0)
//             error->warning(FLERR,"Checkpoint file has different triclinic angles");
//         }
//       }

//     } else error->all(FLERR,"Invalid flag in header section of checkpoint file");

//     flag = read_int();
//   }
// }
