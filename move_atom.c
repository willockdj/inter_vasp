/*********************************************************************/
/* move_atom.c              ; re-position a single ation             */
/*                            Dave Willock, added April 2019         */
/*********************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"


void move_atom(atom *p_atom, double scale, double *move_vec)

{
	p_atom->x += scale* *move_vec;
	move_vec++;
	p_atom->y += scale* *move_vec;
	move_vec++;
	p_atom->z += scale* *move_vec;

return;
}
