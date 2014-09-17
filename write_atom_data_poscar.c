
/* routine to write out an atom data line in a biosym .car format */
#include <stdio.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"


void write_atom_data_poscar(atom *p_atom, double *p_fract, int use_direct)
{

   if (use_direct)
	   {
      printf("%14.9f %14.9f %14.9f\n" , p_atom->x, p_atom->y, p_atom->z);
           }
   else
   {
      printf("%14.9f %14.9f %14.9f \n", *p_fract, *(p_fract+1), *(p_fract+2));  
   }

return;
}

