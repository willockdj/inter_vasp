/*********************************************************************/
/* move_selected_molecule.c ; re-position any molecule fragment      */
/*                            according to indices supplied          */
/*                            Dave Willock, begun Oct 04             */
/*********************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"


void move_selected_atoms(atom *p_molecule, int num_chosen, 
                            int *p_chosen_indices, double *move_vec)

{
int i; 
double *vec_comp;
atom *p_atom;

for (i=0;i<= num_chosen; i++)
    {
        p_atom = p_molecule + *(p_chosen_indices+i);

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
