#include <stdlib.h>

#include "src/system.h"
#include "src/md.h"

int main()
{
  System *system;
  system = (System *)malloc(sizeof(System));

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

  run_md(system);

  free(system->atoms);
  free(system);
  fclose(system->f_trajectory);

  return 0;
}