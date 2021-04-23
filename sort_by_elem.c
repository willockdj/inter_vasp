/************************************************************************/
/* sort_by_elem.c : order atoms in list so they form blocks of elements */
/* started Dave Willock July 2003                                       */
/* malloced Dave Willock January 2016                                   */
/* num_atoms is the real number of atoms so top index is one less       */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"

void sort_by_elem( atom *p_molecule, int num_atoms, atom_number *p_types, int num_types ) 
{
int iatom;
int itype,isorted;

atom *p_atom1;
atom *p_molecule_sorted;
atom *p_atom_sorted;

atom_number *p_this_type;

/***************************************************************************/
/*** Make copy of atoms in sorted order ************************************/
/***************************************************************************/

printf("sorting %d atoms\n", num_atoms);

/***************************************************************************/
/** Malloc local array for sorting *****************************************/
/***************************************************************************/
printf("In sort_by_elem.c : mallocing molecule_sorted for %d atoms\n", num_atoms+1);
printf("Have %d types...\n", num_types);

p_molecule_sorted=(atom*)malloc((num_atoms+1)*sizeof(atom));

isorted=0;
p_this_type=p_types;
p_atom_sorted = p_molecule_sorted;
for (itype=0; itype <= num_types; itype++)
  {
     p_atom1= p_molecule;
     for (iatom=0; iatom <= num_atoms; iatom++)
       {
//          printf("comparing type %s with atom %s element %s\n", 
//                        p_this_type->atom_type, p_atom1->label, p_atom1->elem);

          if (strcmp(&(p_this_type->atom_type[0]), &(p_atom1->elem[0]))==0)
            {
              *p_atom_sorted= *p_atom1;
              p_atom_sorted++;
            }
          p_atom1++;
       }
     p_this_type++;
  }

//printf("Through loop....next copying sorted molecule back\n");
p_atom1= p_molecule;
p_atom_sorted = p_molecule_sorted;
for (iatom=0; iatom <= num_atoms; iatom++)
   {
      *p_atom1 = *p_atom_sorted;
      p_atom1++; p_atom_sorted++;
   }

  free(p_molecule_sorted);
  return;
}
