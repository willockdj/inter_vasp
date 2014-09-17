/************************************************************************/
/* sort_by_elem.c : order atoms in list so they form blocks of elements */
/* started Dave Willock July 2003                                       */
/************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"

void sort_by_elem( atom *p_molecule, int num_atoms, atom_number *p_types, int num_types ) 
{
int iatom,iatom2;
int idummy2, idummy1;
int itype,isorted;

atom *p_atom1, *p_atom2;
atom molecule_sorted[MAXATOMS];

atom_number *p_this_type;

/***************************************************************************/
/*** Make copy of atoms in sorted order ************************************/
/***************************************************************************/

isorted=0;
p_this_type=p_types;
for (itype=0; itype <= num_types; itype++)
  {
     p_atom1= p_molecule;
     for (iatom=0; iatom <= num_atoms; iatom++)
       {
          if (!strcmp(&(p_this_type->atom_type[0]), &(p_atom1->elem[0])))
            {
              molecule_sorted[isorted]= *p_atom1;
              isorted++;
            }
          p_atom1++;
       }
     p_this_type++;
  }

p_atom1= p_molecule;
for (iatom=0; iatom <= num_atoms; iatom++)
   {
      *p_atom1 = molecule_sorted[iatom];
      p_atom1++;
   }

  return;
}
