
/* routine to write out an atom data line in a gulp format */
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"


void write_atom_data_gulp(FILE *fp, atom *p_atom, int is_core)
{
char *p_Al_elem;
char *p_Si_elem;
char *p_O_elem;
char *p_Na_elem;
double occupancy, part_charge;

p_Al_elem= "Al";
p_Si_elem= "Si";
p_O_elem= "O";
p_Na_elem= "Na";

if (is_core)
  {
    if (strncmp(p_Al_elem, &(p_atom->elem[0]), 2) == 0)
      {
        occupancy= 1.0;
        part_charge= 3.0;
      }
    else if (strncmp(p_Si_elem, &(p_atom->elem[0]), 2) == 0)
      {
        occupancy= 1.0;
        part_charge= 4.0;
      }
    else if (strncmp(p_O_elem, &(p_atom->elem[0]), 1) == 0)
      {
        occupancy= 1.0;
        part_charge= -2.0;
      }
    else if (strncmp(p_Na_elem, &(p_atom->elem[0]), 1) == 0)
      {
        occupancy= 1.0;
        part_charge= 1.0;
      }
  }
else
  {
    if (strncmp(p_O_elem, &(p_atom->elem[0]), 1) == 0)
      {
        occupancy= 1.0;
        part_charge= 0.0;
      }
  }
 
if (is_core)
  {
    fprintf(fp, "%-5s core  %14.9f %14.9f %14.9f\n" ,
                  p_atom->elem,
                  p_atom->x, p_atom->y, p_atom->z);
  }
else
  {
    fprintf(fp, "%-5s shel  %14.9f %14.9f %14.9f \n" ,
                  p_atom->elem,
                  p_atom->x, p_atom->y, p_atom->z);
  }

return;
}

