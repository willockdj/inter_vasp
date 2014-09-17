/***************************************************************************/
/* move_molecule_with_flags.c ; re-position any flagged part of a molecule */
/***************************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"


void move_molecule_with_flags(atom *p_molecule, int *p_flag_list, 
                                           int num_atoms, double *move_vec)

{
int i, *p_flag; 
double *vec_comp;
atom *p_atom;

p_atom = p_molecule;
p_flag = p_flag_list;

for (i=0;i<= num_atoms; i++)
    {
      if (*p_flag) 
        {
	   vec_comp = move_vec;
	   p_atom->x += *vec_comp;
	   vec_comp++;
	   p_atom->y += *vec_comp;
	   vec_comp++;
	   p_atom->z += *vec_comp;
        }
        p_atom++;
        p_flag++;
    }

return;
}
