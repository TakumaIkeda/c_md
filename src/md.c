#include <stdio.h>
#include <stdlib.h>
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

double calc_potential(System *system)
{
  Atom *atoms = system->atoms;
  double dx[3];
  double eps = 1.0;
  double sigma = 0.398; // nm
  double sigma2 = sigma * sigma;
  double sigma6 = sigma2 * sigma2 * sigma2;
  double sigma12 = sigma6 * sigma6;
  double potential = 0.0;

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

      if (r < system->cut_off)
      {
        potential += 4 * eps * (sigma12 / r12 - sigma6 / r6);
      }
    }
  }

  return potential;
}

double calc_kinetic(System *system)
{
  Atom *atoms = system->atoms;
  double kinetic = 0.0;

  for (int i = 0; i < system->natoms; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      kinetic += 0.5 * atoms[i].p[j] * atoms[i].p[j];
    }
  }

  return kinetic;
}

void run_md()
{
  System system;

  init_system(&system);
  calc_force(&system);
  for (int step = 0; step < system.nsteps; step++)
  {
    calc_momentum(&system);
    calc_position(&system);
    calc_force(&system);
    calc_momentum(&system);
    if (step % 10 == 0)
    {
      double potential = calc_potential(&system);
      double kinetic = calc_kinetic(&system);
      printf("Step %d : Potential= %f Kinetic= %f Total= %f\n", step, potential, kinetic, potential + kinetic);
      output_gro(&system, step, system.f_trajectory);
    }
  }

  free(system.atoms);
  fclose(system.f_trajectory);

  return;
}
