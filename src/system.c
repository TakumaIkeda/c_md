#include <stdlib.h>
#include <math.h>
#include "system.h"
#include "utils.h"

void init_system(System *system)
{
  system->box[0] = 4.2; // 4.2 nm
  system->box[1] = 4.2; // 4.2 nm
  system->box[2] = 4.2; // 4.2 nm
  system->temperature = 50;
  system->dt = 0.001; // 1 fs
  system->nsteps = 10000;
  system->natoms = 1000;
  system->cut_off = 2.0; // 2.0 nm
  system->atoms = (Atom *)malloc(system->natoms * sizeof(Atom));
  system->f_trajectory = fopen("trajectory.gro", "w");
  system->f_energy = fopen("energy.xvg", "w");
  system->f_temperature = fopen("temperature.xvg", "w");
  Atom *atoms = system->atoms;

  double binx = system->box[0] / 10;
  double biny = system->box[1] / 10;
  double binz = system->box[2] / 10;

  srand(42);

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
      }
    }
  }

  double total_px = 0.0;
  double total_py = 0.0;
  double total_pz = 0.0;
  for (int i = 0; i < system->natoms; i++)
  {
    double sigma = sqrt(system->temperature);
    atoms[i].p[0] = sigma * (rand() / (double)RAND_MAX * 2 - 1);
    atoms[i].p[1] = sigma * (rand() / (double)RAND_MAX * 2 - 1);
    atoms[i].p[2] = sigma * (rand() / (double)RAND_MAX * 2 - 1);

    total_px += atoms[i].p[0];
    total_py += atoms[i].p[1];
    total_pz += atoms[i].p[2];
  }

  for (int i = 0; i < system->natoms; i++)
  {
    atoms[i].p[0] -= total_px / system->natoms;
    atoms[i].p[1] -= total_py / system->natoms;
    atoms[i].p[2] -= total_pz / system->natoms;
  }

  output_gro(system, 0, system->f_trajectory);

  return;
}