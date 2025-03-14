#ifndef MD_H
#define MD_H

#include "system.h"

void run_md(System *system);
void calc_momentum(System *system);
void calc_position(System *system);

#endif