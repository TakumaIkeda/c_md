#ifndef UTILS_H
#define UTILS_H

double norm(double *v);
void pbc(Atom *atom, double *box);
void output_gro(System *system, int step, FILE *file);

#endif