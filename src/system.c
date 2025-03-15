#include <stdlib.h>
#include "system.h"
#include "utils.h"

void init_system(System *system)
{
  system->box[0] = 4.2; // 4.2 nm
  system->box[1] = 4.2; // 4.2 nm
  system->box[2] = 4.2; // 4.2 nm
  system->temperature = 300;
  system->dt = 0.001; // 1 fs
  system->nsteps = 1000;
  system->natoms = 1000;
  system->cut_off = 2.0; // 2.0 nm
  system->atoms = (Atom *)malloc(system->natoms * sizeof(Atom));
  system->f_trajectory = fopen("trajectory.gro", "w");
  Atom *atoms = system->atoms;

  double binx = system->box[0] / 10;
  double biny = system->box[1] / 10;
  double binz = system->box[2] / 10;

  // printf("binx = %f, biny = %f, binz = %f\n", binx, biny, binz);

  for (int idx = 0; idx < 10; idx++)
  {
    for (int idy = 0; idy < 10; idy++)
    {
      for (int idz = 0; idz < 10; idz++)
      {
        int id = idx * 100 + idy * 10 + idz;
        atoms[id].q[0] = idx * binx;
        atoms[id].q[1] = idy * biny;
        atoms[id].q[2] = idz * binz;
        atoms[id].p[0] = 0.0;
        atoms[id].p[1] = 0.0;
        atoms[id].p[2] = 0.0;
        if (id == 0)
        {
          atoms[id].p[0] = 1.0;
          atoms[id].p[1] = 1.0;
          atoms[id].p[2] = 1.0;
        }
      }
    }
  }

  output_gro(system, 0, system->f_trajectory);

  return;
}