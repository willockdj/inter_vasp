/*************************************************************/
/* find_mol.c : scan atom list and assign molecule indices  **/
/******         according to atom groups that are bonded to **/
/******         one another.                                **/
/******                     started June 04 Dave Willock    **/
/****** Last updated April 06, removed problem with an      **/
/****** atom being added to ends list multiple times.       **/
/*************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "debug.h"
#include "structures.h"

int find_mol(atom *p_molecule, int num_atoms)
{
#include "header.h"
int ineigh, iflag, num_ends, atoms_at_ends[MAX_ENDS];
int *p_flag, iline;
int *p_added, *p_this_added;
int this_neigh, iatom, iends, mol_index;
atom *p_current_atom, *p_atom, *p_neigh;

/********* zero the mol indicies *****/

printf("In find_mol have %d atoms\n", num_atoms);

/***** malloc space for array to keep track of atoms added to molecule ***/
p_added = (int*)malloc((num_atoms+1)*sizeof(int));

p_atom= p_molecule; p_this_added=p_added;
for (iatom=0; iatom <= num_atoms; iatom++)
    {
      p_atom->mol = 0;
      *p_this_added=FALSE;
      p_atom++; p_this_added++;
    }

printf("added list updated\n");

/******** set atom2 as the first member of the chunk *********/

mol_index=0;
p_atom= p_molecule; p_this_added=p_added;
for (iatom=0; iatom<= num_atoms; iatom++)
  {
//    printf("iatom=%d\n", iatom);
    if (p_atom->mol == 0 )
      {
        mol_index++;

/***** UNCLOSED ****/

        num_ends= 0;
        atoms_at_ends[0]= iatom;
        *p_this_added=TRUE;

/******** set the mol index for all end atoms    ******/
/******** and generate new ends                  ******/
/******** till we have all members of this molecule ***/

        while (num_ends >= 0)
          {
//            printf("Atom at end is: %d\n", atoms_at_ends[0]);
            p_current_atom= p_molecule + atoms_at_ends[0];
            p_current_atom->mol=mol_index;

            if (debug)
              {
                printf("\n\nCurrent atom %d %s\n", atoms_at_ends[0], p_current_atom->label);
                printf("Has %d neighbours\n", p_current_atom->num_neigh);
              }

            for (ineigh = 0; ineigh <  (p_current_atom->num_neigh); ineigh++)
              {
                this_neigh=  p_current_atom->neighb[ineigh];
                p_neigh= p_molecule+this_neigh; 

                if (debug) printf("Neigh %d is %d %s\n", ineigh, this_neigh, p_neigh->label);

/***************************************************************************/
/***** add this neighbour to the ends list if it is not already done    ****/
/***************************************************************************/

               if ( p_neigh->mol == 0 )
                 { 
                   if (num_ends < MAX_ENDS && p_neigh->num_neigh != 1 && !*(p_added+this_neigh)) 
                     {
                       if (debug) printf("Adding this to the ends list\n");

                       num_ends++;
                       *(p_added+this_neigh)=TRUE;
                       atoms_at_ends[num_ends]= this_neigh;       
                     }
                   else if (num_ends >= MAX_ENDS)
                     {
                       printf("ERROR: Run out of ends in find_mol.c increase MAX_ENDS\n");
                       printf("num_ends = %d\n", num_ends);
                       exit(1);
                     }
                   else
                     {
                       if (debug) printf("Assigning this to mol_index:%d\n", mol_index);

                       p_neigh->mol = mol_index;
                     }
                 }
             } 

/******* shuffle atoms_at_ends list to cover the one we have dealt with ****/

            if (num_ends >= MAX_ENDS)
              {
                printf("ERROR: Run out of ends in find_mol.c increase MAX_ENDS\n");
                printf("num_ends = %d\n", num_ends);
                exit(1);
              }

            if (debug)
              {
                iline=0;
                printf("\nCurrent ends list:\n");
              }

            for (iends= 0; iends < num_ends; iends++)
              {
                  atoms_at_ends[iends] =  atoms_at_ends[iends+1];   

                  if (debug)
                    {
                       printf("%d ",atoms_at_ends[iends]); 

                       if (iline < 10)
                          {
                            iline++;
                          }
                       else
                          {
                            iline=0;
                            printf("\n");
                          }
                    }
              }
            num_ends--;
         }
     }
   p_atom++; p_this_added++;
  }

/* Free malloced arrays ***************/

 free(p_added);
/**************************************/

return mol_index;
}

		
