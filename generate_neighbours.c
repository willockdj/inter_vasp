/**********************************************************************/
/* generate_neighbours.c : search molecule for neighbour list         */
/* NOTE: neighb is referenced from ZERO                               */
/* started Dave and Dewi 23/3/95                                      */
/**********************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "global_values.h"
#include "maxima.h"
#include "structures.h"

double standard_bond( char *atom1, char *atom2 );

double atom_separation_squared(atom *p_A, atom *p_B, int pbc,  double *p_recip_latt_vec, double *p_latt_vec);

void generate_neighbours( atom *p_molecule, int num_atoms, 
                          atom_number *p_types, int *p_num_types,
                          int use_pbc,  double *p_recip_latt_vec, double *p_latt_vec,
                          charge_list *p_spec_charges, int set_labels)
{
int iatom1,iatom2;
int idummy2, idummy1;
int itype;

double standard,actual_dist;

atom *p_atom1, *p_atom2;

atom_number *p_this_type;

if (set_labels)
  printf("generate_neighbours being used to generate atom labels\n");


/********* Work out how many of each atom are present in this molecule *****/
/********* and use types as GULP species                               *****/

*p_num_types= 0;
strcpy( &(p_types->atom_type[0]), &(p_molecule->elem[0]));
strcpy( &(p_spec_charges->label[0]), &(p_molecule->elem[0]));

p_spec_charges->is_core = TRUE;
p_spec_charges->part_chge = p_molecule->part_chge;

iatom1=0;
p_atom1= p_molecule;
/*** Debug Dec 08, set number of type 1 to zero as it is declared ***/
p_types->num = 0;
//
if (set_labels)
  {
    if (strlen(p_atom1->elem) == 1)
       sprintf(p_atom1->label,"%1s%d", p_atom1->elem,p_types->num+1);
    else if (strlen(p_atom1->elem) == 2)
       sprintf(p_atom1->label,"%2s%d", p_atom1->elem,p_types->num+1);
    else
       printf("elem doesn't fit, strlen = %d\n", strlen(p_atom1->elem));
  }

/*** Loop of rest of atoms ****/

for (iatom1 = 1; iatom1 <= num_atoms; iatom1++)
   {

     p_atom1++;
     p_this_type= p_types-1;

     for (itype=0; itype<= *p_num_types; itype++)
       {
         p_this_type++;
         if (!strcmp(&(p_this_type->atom_type[0]), &(p_atom1->elem[0]))) break;
       }

/**** test to see if that was a known element*****/

     if (itype <= *p_num_types)
       {
         (p_this_type->num)++;
         if (set_labels)
           {
              sprintf(p_atom1->label,"%s%d",p_atom1->elem,p_this_type->num+1);
           }
       }
     else
       {
         (*p_num_types)++;
         p_this_type++;
         strcpy( &(p_this_type->atom_type[0]), &(p_atom1->elem[0]));
         p_this_type->num= 0;
//
         if (set_labels)
           {
              sprintf(p_atom1->label,"%s%d",p_atom1->elem,p_this_type->num+1);
           }

         p_spec_charges++;
         strcpy( &(p_spec_charges->label[0]), &(p_atom1->elem[0]));
         p_spec_charges->is_core = TRUE;
         p_spec_charges->part_chge = p_atom1->part_chge;
       }

   }

/**********************************************************************/
/******************* loop over pairs of atoms in molecule *************/
/*******************         looking for bonds            *************/
/**********************************************************************/

/**********************************************************************/
/* Set atoms num_neighs to zero outside of double loop to *************/
/* ensure that neighbour numbers are counted correctly.   *************/
/* Note this routine is out of line with Zebedde at the   *************/
/* moment in the use of a 0 index here.                   *************/
/* Dave Willock, Nov. 2018                                *************/
/**********************************************************************/
for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
    p_atom1= p_molecule+iatom1;
    p_atom1->num_neigh=0;
  }

printf("In generate neighbours with %d atoms...\n", num_atoms);
printf("Have found %d types....................\n", *p_num_types);

for (iatom1 = 0; iatom1 <= num_atoms; iatom1++)
  {
  p_atom1= p_molecule+iatom1;

  for (iatom2=iatom1+1; iatom2 <= num_atoms; iatom2++)
     {
      p_atom2= p_molecule+iatom2;

      standard= standard_bond( &(p_atom1->elem[0]), &(p_atom2->elem[0]));
      if (standard != -1) 
        {
           actual_dist = atom_separation_squared(p_atom1, p_atom2, 
                                                  use_pbc, p_recip_latt_vec, p_latt_vec);
           actual_dist = sqrt(actual_dist);

/* test bond length against standard */

	   if (actual_dist <= standard+BOND_TOL)
             {
                 if (idummy1 < MAX_NEIGHS)
                   {
		      idummy1= (p_atom1->num_neigh)++;
	              p_atom1->neighb[idummy1]= iatom2;
                   }
                 else
                   {
                      printf("Warning: Number of neighbours for atom %s exceeds MAX_NEIGHS, (%d), ignoring rest.\n",
                                   p_atom1->label, MAX_NEIGHS);
                   }

                 if (idummy1 < MAX_NEIGHS)
                   {
		      idummy2= (p_atom2->num_neigh)++;
	              p_atom2->neighb[idummy2]= iatom1; 
                   }
                 else
                   {
                      printf("Warning: Number of neighbours for atom %s exceeds MAX_NEIGHS, (%d), ignoring rest.\n",
                                   p_atom2->label, MAX_NEIGHS);
                   }

//                 printf("iatom1 %d iatom2 %d, idummies %d %d\n", iatom1, iatom2, idummy1, idummy2);
             }

/* end of if (Standard) */
	 }
     }
  } 
  return;
}
