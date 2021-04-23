/******************************************************************/
/** Routine to count oxygen species.                        *******/
/** Begun May 2019, Dave, Mala, Connie and Sachin           *******/
/******************************************************************/
/** p_molecule is a pointer to the start of the atom list   *******/
/** num_atoms is the total number of atoms in the list      *******/
/******************************************************************/
/*** num_oxy_species will hold the number of each oxygen species found ***/
/*** num_oxy_species[0] is for O lattice ***/
/*** num_oxy_species[1] is for OH        ***/
/*** num_oxy_species[2] is for water O   ***/
/*** num_oxy_species[3] is for O alcohol ***/
/******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"
#include "global_values.h"

/* prototype list for this routine */

/*--------------------------------------------------------------------------------------------*/

void count_oxygen_species(atom *p_molecule, int num_atoms, int *oxy_list_ptrs[],
                          atom_number *p_num_oxy_species, int num_oxy_types, int just_count)
  {
  int iatom, ineigh, neigh_index, isH[2], countH;
  int ispecies; 
  atom_number *p_this_num_species;
  atom *p_atom, *p_neigh;

/* Zero arrays */
  p_this_num_species= p_num_oxy_species;
  for (ispecies=0; ispecies<num_oxy_types; ispecies++)
    {
       p_this_num_species->num = 0;
       p_this_num_species++;
    }

  printf("\n\nCount oxygen species routine.....\n\n");

  p_atom=p_molecule;
  for (iatom=0; iatom< num_atoms; iatom++)
    {
      printf("%s %10.6f %10.6f %10.6f \n", p_atom->label,
                                           p_atom->x,
                                           p_atom->y,
                                           p_atom->z );


      printf("This atom has %d neighbours....\n", p_atom->num_neigh);

      if ( strcmp(p_atom->elem, "O") == 0)
        {
          printf("This is an Oxygen atom\n");

          if ( p_atom->num_neigh == 0 )
            {
              printf("This is a lattice oxygen.....\n");
              if (!just_count) oxy_list_ptrs[0][p_num_oxy_species->num] = iatom; 
              ++(p_num_oxy_species->num);
            }
          else
            {
               if (p_atom->num_neigh == 1 )
                 {
                    neigh_index= p_atom->neighb[0];
                    p_neigh = p_molecule + neigh_index;

                    if ( strcmp(p_neigh->elem, "H") == 0 )
                      {
                        printf("This is a hydroxyl oxygen.\n");
                        if (!just_count) oxy_list_ptrs[1][(p_num_oxy_species+1)->num] = iatom; 
                        ++((p_num_oxy_species+1)->num);
                      }
                    else
                      {
                        printf("ERROR: This is an oxygen with a single neighbour\n");
                        printf("ERROR: but not H, don't know how to classify.....\n");
                        exit(0);
                      }
                 }
               else if (p_atom->num_neigh == 2 )
                 {
                   for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                     {
                       neigh_index= p_atom->neighb[ineigh];
                       p_neigh = p_molecule + neigh_index;

                       isH[ineigh] = strcmp(p_neigh->elem, "H") == 0;

                       printf("looking at neighbour %d label %s elem %s get isH = %d\n",
                                    ineigh, p_neigh->label, p_neigh->elem, isH[ineigh]);
                     }

                   countH = isH[0]+isH[1];

                   if (countH == 0)
                     {
                       printf("ERROR: This is an oxygen with two neighbours\n");
                       printf("ERROR: but neither is H, cannot classify....\n");
                       exit(0);
                     }
                   else if (countH == 1)
                     {
                        if ( isH[0] == 1 )
                          {
                            neigh_index= p_atom->neighb[1];
                            p_neigh = p_molecule + neigh_index;
                          }
                        else
                          {
                            neigh_index= p_atom->neighb[0];
                            p_neigh = p_molecule + neigh_index;
                          }

                        if ( strcmp(p_neigh->elem, "C") == 0 )
                          {
                            printf("This is an oxygen in an alcohol.\n");
                            if (!just_count) oxy_list_ptrs[3][(p_num_oxy_species+3)->num] = iatom; 
                            ++((p_num_oxy_species+3)->num);
                          }
                        else
                          {
                            printf("ERROR: This is an oxygen with two neighbours\n");
                            printf("ERROR: one is H but the other is            \n");
                            printf("ERROR: neither H or C, cannot classify....\n");
                            exit(0);
                          }
                     }
                   else if (countH == 2)
                     {
                       printf("This is an oxygen in water.\n");
                       if (!just_count) oxy_list_ptrs[2][(p_num_oxy_species+2)->num] = iatom; 
                       ++((p_num_oxy_species+2)->num);
                     }
                   else
                     {
                       printf("ERROR: Count H is greater than two.\n");
                       printf("ERROR: This is impossible.............\n");
                       exit(0);
                     }
                 }
               else
                 {
                    printf("ERROR: This is an oxygen with more than two neighbours\n");
                    printf("ERROR: don't know how to classify.....\n");
                    exit(0);
                 }

               for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                 {
                   neigh_index= p_atom->neighb[ineigh];
                   printf("%s elem is %s (%d)",(p_molecule+neigh_index)->label, (p_molecule+neigh_index)->elem, neigh_index);
                 }
               printf("\n");
            }
       }
      else
       {
         printf("This is not an O atom.....\n");
       }

      p_atom++;
    }

  return;
  }
