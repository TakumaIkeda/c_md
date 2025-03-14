#ifndef ATOM_H
#define ATOM_H

typedef struct t_atom
{
  double q[3], p[3], f[3];
} Atom;

void pbc_dx(Atom *atom1, Atom *atom2, double *dx, double *box);

#endif