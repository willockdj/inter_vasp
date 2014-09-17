/*********************************************************************/
/* move_molecule.c ; re-position any molecule fragment               */
/*********************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"


void move_molecule(atom *p_molecule, int num_atoms, double *move_vec)

{
int i; 
double *vec_comp;
atom *p_atom;

p_atom = p_molecule;

for (i=0;i<= num_atoms; i++)
    {

	vec_comp = move_vec;
	p_atom->x += *vec_comp;
	vec_comp++;
	p_atom->y += *vec_comp;
	vec_comp++;
	p_atom->z += *vec_comp;
        p_atom++;
    }

return;
}
