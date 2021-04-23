/*********************************************************************/
/* inter_atom_vector.c simply calcuates the vector between two atoms */
/*                     This will be the vector from p_atom to p_atom2*/
/* added April 2019, Dave Willock                                    */
/*********************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"


void inter_atom_vector(atom *p_atom, atom *p_atom2, double *p_vec)
{
    *p_vec =  p_atom2->x-p_atom->x; p_vec++;
    *p_vec =  p_atom2->y-p_atom->y; p_vec++;
    *p_vec =  p_atom2->z-p_atom->z; 

return;
}
