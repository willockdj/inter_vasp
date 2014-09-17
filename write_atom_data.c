
/* routine to write out an atom data line in a biosym .car format */
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

void write_atom_data(FILE *fp, atom *p_atom, double scale_factor, coord_flags *p_fix_flags)
{
int have_fixed;

      have_fixed = FALSE;

      if (p_fix_flags->fx || p_fix_flags->fy || p_fix_flags->fz) 
        {
          strcpy(p_atom->group, "F");
          have_fixed=TRUE;
        }

      if (have_fixed)
        {
          if (p_fix_flags->fx) strcat(p_atom->group, "F"); else strcat(p_atom->group, "T"); 
          if (p_fix_flags->fy) strcat(p_atom->group, "F"); else strcat(p_atom->group, "T"); 
          if (p_fix_flags->fz) strcat(p_atom->group, "F"); else strcat(p_atom->group, "T"); 
        }

      fprintf(fp, "%-5s %14.9f %14.9f %14.9f %-4s %-7s%-7s %-2s %6.3f\n" ,
               p_atom->label,
               p_atom->x * scale_factor, 
               p_atom->y * scale_factor, 
               p_atom->z * scale_factor,
               p_atom->group,
               p_atom->group_no,              
               p_atom->pot,
               p_atom->elem,
               p_atom->part_chge);

return;
}

