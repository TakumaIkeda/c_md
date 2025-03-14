#include <stdio.h>
#include "system.h"
#include "utils.h"

void calc_force(System *system)
{
  Atom *atoms = system->atoms;
  double dx[3];
  double eps = 1.0;
  double sigma = 0.398; // nm
  double sigma2 = sigma * sigma;
  double sigma6 = sigma2 * sigma2 * sigma2;
  double sigma12 = sigma6 * sigma6;

  for (int i = 0; i < system->natoms; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      atoms[i].f[j] = 0.0;
    }
  }

  for (int i = 0; i < system->natoms; i++)
  {
    for (int j = 0; j < system->natoms; j++)
    {
      if (i == j)
      {
        continue;
      }

      pbc_dx(&atoms[i], &atoms[j], dx, system->box);
      double r = norm(dx);
      double r2 = r * r;
      double r6 = r2 * r2 * r2;
      double r12 = r6 * r6;

      // printf("r = %f, r2 = %f, r6 = %f, r12 = %f\n", r, r2, r6, r12);
      if (r < system->cut_off)
      {
        double f = 4 * eps * (12 * sigma12 / r12 - 6 * sigma6 / r6) / r2;
        for (int k = 0; k < 3; k++)
        {
          atoms[i].f[k] += f * dx[k] / r;
        }
      }
    }
  }

  return;
}

void calc_momentum(System *system)
{
  Atom *atoms = system->atoms;

  for (int i = 0; i < system->natoms; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      atoms[i].p[j] += 0.5 * system->dt * atoms[i].f[j];
    }
  }

  return;
}

void calc_position(System *system)
{
  Atom *atoms = system->atoms;

  for (int i = 0; i < system->natoms; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      atoms[i].q[j] += system->dt * atoms[i].p[j];
      pbc(&atoms[i], system->box);
    }
  }

  return;
}

void run_md(System *system)
{
  // init system
  init_system(system);
  calc_force(system);
  // run md
  for (int step = 0; step < system->nsteps; step++)
  {
    calc_momentum(system);
    calc_position(system);
    calc_force(system);
    if (step % 10 == 0)
    {
      printf("Step %d\n", step);
      output_gro(system, step, system->f_trajectory);
    }
  }

  // output results

  return;
}
