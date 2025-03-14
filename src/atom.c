#include "atom.h"

void pbc_dx(Atom *atom1, Atom *atom2, double *dx, double *box)
{
  for (int i = 0; i < 3; i++)
  {
    dx[i] = atom1->q[i] - atom2->q[i];
    while (dx[i] > 0.5 * box[i])
    {
      dx[i] -= box[i];
    }

    while (dx[i] < -0.5 * box[i])
    {
      dx[i] += box[i];
    }
  }

  return;
}
