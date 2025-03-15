#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "system.h"
#include "utils.h"

void calc_force(System *system, int rank, int size)
{
  Atom *atoms = system->atoms;
  double dx[3];
  double eps = 1.0;
  double sigma = 0.398; // nm
  double sigma2 = sigma * sigma;
  double sigma6 = sigma2 * sigma2 * sigma2;
  double sigma12 = sigma6 * sigma6;
  int nranks = 6;
  double *f_system = (double *)malloc(system->natoms * 3 * sizeof(double));

  for (int i = 0; i < system->natoms * 3; i++)
  {
    f_system[i] = 0.0;
  }

  int start = rank * system->natoms / nranks;
  int end = (rank == size - 1) ? system->natoms : (rank + 1) * system->natoms / nranks;

  for (int i = start; i < end; i++)
  {
    for (int j = i + 1; j < system->natoms; j++)
    {
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
          f_system[i * 3 + k] += f * dx[k] / r;
          f_system[j * 3 + k] -= f * dx[k] / r;
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE, f_system, system->natoms * 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for (int i = 0; i < system->natoms; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      atoms[i].f[j] = f_system[i * 3 + j];
    }
  }

  return;
}

void calc_momentum(System *system, int rank, int size)
{
  Atom *atoms = system->atoms;

  if (rank == 0)
  {
    for (int i = 0; i < system->natoms; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        atoms[i].p[j] += 0.5 * system->dt * atoms[i].f[j];
      }
    }
  }

  MPI_Bcast(atoms, system->natoms * sizeof(Atom), MPI_BYTE, 0, MPI_COMM_WORLD);

  return;
}

void calc_position(System *system, int rank, int size)
{
  Atom *atoms = system->atoms;

  if (rank == 0)
  {
    for (int i = 0; i < system->natoms; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        atoms[i].q[j] += system->dt * atoms[i].p[j];
        pbc(&atoms[i], system->box);
      }
    }
  }

  MPI_Bcast(atoms, system->natoms * sizeof(Atom), MPI_BYTE, 0, MPI_COMM_WORLD);

  return;
}

void run_md()
{
  System system;
  int rank;
  int size;

  // init MPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("MPI rank %d, size %d\n", rank, size);

  // init system
  init_system(&system);
  calc_force(&system, rank, size);
  // run md
  for (int step = 0; step < system.nsteps; step++)
  {
    calc_momentum(&system, rank, size);
    calc_position(&system, rank, size);
    calc_force(&system, rank, size);
    if (step % 100 == 0 && rank == 0)
    {
      printf("Step %d\n", step);
    }

    if (step % 10 == 0 && rank == 0)
    {
      output_gro(&system, step + 1, system.f_trajectory);
    }
  }

  // output results
  free(system.atoms);
  fclose(system.f_trajectory);
  MPI_Finalize();

  return;
}
