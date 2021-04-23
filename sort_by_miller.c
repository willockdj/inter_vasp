/************************************************************************/
/* sort_by_miller : order atoms by their position along a vector        */
/*                  defined by a Miller index set.                      */
/* started Dave Willock March 2019                                      */
/* num_atoms is the real number of atoms so top index is one less       */
/************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"

void min_image( double *x, double *y, double *z, double *p_recip_latt, double *p_latt_vec);

void sort_by_miller( atom *p_molecule, int num_atoms, atom_number *p_types, int num_types,
                     double *p_latt_vec, double *p_recip_latt, int *p_miller ) 
{
int iatom, jatom;
int itype, start_type, isorted;

double surf_norm[3], x,y,z;
double *p_dots, *p_this_dot, *p_this_dot2, dummy_dot;

atom_number *p_this_type;

atom *p_atom1,*p_atom2;
atom *p_molecule_sorted;
atom *p_atom_sorted;

atom dummy_atom;

/***************************************************************************/
/*** Make copy of atoms in sorted order ************************************/
/***************************************************************************/
printf("In sort_by_miller have lattice vecs:\n");
printf("%10.6f %10.6f %10.6f \n", *p_latt_vec, *(p_latt_vec+1), *(p_latt_vec+2));
printf("%10.6f %10.6f %10.6f \n", *(p_latt_vec+3), *(p_latt_vec+4), *(p_latt_vec+5));
printf("%10.6f %10.6f %10.6f \n", *(p_latt_vec+6), *(p_latt_vec+7), *(p_latt_vec+8));
printf("In sort_by_miller have recip lattice vecs:\n");
printf("%10.6f %10.6f %10.6f \n", *p_recip_latt, *(p_recip_latt+1), *(p_recip_latt+2));
printf("%10.6f %10.6f %10.6f \n", *(p_recip_latt+3), *(p_recip_latt+4), *(p_recip_latt+5));
printf("%10.6f %10.6f %10.6f \n", *(p_recip_latt+6), *(p_recip_latt+7), *(p_recip_latt+8));

printf("sorting %d atoms by Miller index %d %d %d\n", num_atoms, *p_miller, *(p_miller+1), *(p_miller+2));

surf_norm[0] = *p_miller *(*p_recip_latt)    + *(p_miller+1)* *(p_recip_latt+3) + *(p_miller+2)* *(p_recip_latt+6);
surf_norm[1] = *p_miller * *(p_recip_latt+1) + *(p_miller+1)* *(p_recip_latt+4) + *(p_miller+2)* *(p_recip_latt+7);
surf_norm[2] = *p_miller * *(p_recip_latt+2) + *(p_miller+1)* *(p_recip_latt+5) + *(p_miller+2)* *(p_recip_latt+8);

printf("Surface Normal vector: %10.6f  %10.6f  %10.6f\n\n", surf_norm[0], surf_norm[1], surf_norm[2]);

/***************************************************************************/
/** Malloc local array for sorting *****************************************/
/***************************************************************************/
printf("In sort_by_miller.c : mallocing molecule_sorted for %d atoms\n", num_atoms+1);

p_molecule_sorted=(atom*)malloc((num_atoms+2)*sizeof(atom));
p_dots           =(double*)malloc((num_atoms+2)*sizeof(double));


/*** copy atoms and min_image all co-ordinates **/
printf("Doing dots......\n");
p_this_dot=p_dots;
p_atom1= p_molecule; p_atom_sorted=p_molecule_sorted;
for (iatom=0; iatom <= num_atoms; iatom++)
  {
     *p_atom_sorted= *p_atom1;
     min_image( &(p_atom_sorted->x), &(p_atom_sorted->y), &(p_atom_sorted->z), p_recip_latt, p_latt_vec);

     printf("Atom %d %s coords %10.6f  %10.6f  %10.6f \n", iatom+1, p_atom1->label,  
                                                           p_atom1->x, 
                                                           p_atom1->y, 
                                                           p_atom1->z);

/* Work out dot product of this atom with surface normal */
     *p_this_dot= surf_norm[0]* p_atom_sorted->x
                 +surf_norm[1]* p_atom_sorted->y
                 +surf_norm[2]* p_atom_sorted->z;

     p_atom_sorted++; p_atom1++; p_this_dot++;
  }

/*** report initial list **/
printf("Initial dot_list:\n");
p_atom_sorted= p_molecule_sorted; p_this_dot=p_dots;
for (iatom=0; iatom <= num_atoms; iatom++)
  {
     printf("Atom %d %s coords %10.6f  %10.6f  %10.6f dot: %10.6f\n", iatom+1, p_atom_sorted->label,  
                                                                      p_atom_sorted->x, 
                                                                      p_atom_sorted->y, 
                                                                      p_atom_sorted->z,
                                                                      *p_this_dot);
     p_atom_sorted++; p_this_dot++;
  }
/*** Sort atoms within types groups by the dot list ***/
start_type=0;
p_atom1= p_molecule_sorted; p_this_dot=p_dots; p_this_type= p_types;
for (itype=0; itype <= num_types; itype++)
  {
    printf("Starting sorting loop, start_type= %d\n", start_type);
    for (iatom=start_type; iatom <= start_type+p_this_type->num; iatom++)
      { 
        p_atom2=p_molecule_sorted+iatom+1; p_this_dot2=p_dots+iatom+1;
        for (jatom=iatom+1; jatom <= start_type+p_this_type->num; jatom++)
          {
/** Bubble **/
             if ( *p_this_dot > *p_this_dot2 ) 
               {
                  dummy_dot   = *p_this_dot2;
                  *p_this_dot2= *p_this_dot;
                  *p_this_dot = dummy_dot;

                  dummy_atom    = *p_atom2;
                  *p_atom2 = *p_atom1;
                  *p_atom1 = dummy_atom;
               }
             p_this_dot2++; p_atom2++;
          }
        p_this_dot++; p_atom1++;
      }
    start_type+=p_this_type->num+1;
    p_this_type++;
  }
 
printf("After sort:\n");
p_atom_sorted= p_molecule_sorted; p_this_dot=p_dots;
for (iatom=0; iatom <= num_atoms; iatom++)
  {
     printf("Atom %d %s coords %10.6f  %10.6f  %10.6f dot: %10.6f\n", iatom+1, p_atom_sorted->label,  
                                                                      p_atom_sorted->x, 
                                                                      p_atom_sorted->y, 
                                                                      p_atom_sorted->z,
                                                                      *p_this_dot);
     p_atom_sorted++; p_this_dot++;
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
  free(p_dots);
  return;
}
