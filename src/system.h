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
} System;

void init_system(System *system);
void calc_force(System *system);

#endif