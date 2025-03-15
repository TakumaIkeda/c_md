#ifndef SYSTEM_H
#define SYSTEM_H

#include <stdio.h>
#include "atom.h"

typedef struct t_system
{
  double box[3];
  int natoms;
  double temperature;
  double dt;
  int nsteps;
  double cut_off;
  Atom *atoms;
  FILE *f_trajectory;
  FILE *f_energy;
  FILE *f_temperature;
} System;

void init_system(System *system);

#endif