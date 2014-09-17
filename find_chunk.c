/***************************************************************/
/* find_chunk.c : scan template and and flag all atoms in the **/
/******           same chunk of molecule as atom2 but not     **/
/******           atom1                                       **/
/******     started May 95 Dave and Dewi                      **/
/***************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

int find_chunk(atom *p_molecule, int num_atoms, int *p_flag_list,
                int atom1, int atom2)
{
#include "header.h"
int ineigh, iflag, num_ends, atoms_at_ends[MAX_ENDS];
int *p_flag;
int this_neigh, iatom, iends;
atom *p_current_atom, *p_atom;

/********* zero the flags ************/

p_flag= p_flag_list; 

for (iflag=0; iflag<= num_atoms; iflag++)
	{
      *p_flag = FALSE;
      p_flag++;
    }

/******** set atom2 as the first member of the chunk *********/

p_flag =  p_flag_list+  atom2;
*p_flag = 1;

num_ends= 0;
atoms_at_ends[0]= atom2;

/******** set the flags for all end atoms and generate new ends ******/
/******** till we have been everywhere in the chunk             ******/

while (num_ends >= 0)
  {
    p_current_atom= p_molecule + atoms_at_ends[0];

    for (ineigh = 0; ineigh <  (p_current_atom->num_neigh); ineigh++)
      {
        this_neigh=  p_current_atom->neighb[ineigh];
        p_flag=  p_flag_list+ this_neigh;

/***************************************************************************/
/***** If atom1 crops up as a neighbour of someone else in this ************/
/***** chunk then atoms 1 and 2 must be in a ring               ************/
/***************************************************************************/

        if ( this_neigh == atom1 && atoms_at_ends[0] != atom2)
           {
             return FALSE;
           }

/***************************************************************************/
/***** add this neighbour to the ends list if it is not already flagged ****/
/***** dissallow atom1 joining the list to stop the other chunk being   ****/
/***** crossed into                                                     ****/
/***************************************************************************/

        if ( ! *p_flag &&  this_neigh != atom1 )
          { 
            *p_flag= TRUE;
            if (num_ends < MAX_ENDS && (p_molecule+this_neigh)->num_neigh != 1) 
              {
                num_ends++;
                atoms_at_ends[num_ends]= this_neigh;       
              }
            else if (num_ends > MAX_ENDS)
              {
                printf("Run out of ends in find_chunk.c increase MAX_ENDS");
                exit(1);
              }
          }
      } 

/******* shuffle atoms_at_ends list to cover the one we have dealt with ****/

    for (iends= 0; iends < num_ends; iends++)
      {
          atoms_at_ends[iends] =  atoms_at_ends[iends+1];   
      }
    num_ends--;
  }

/**************************************/

return TRUE;
}

		
