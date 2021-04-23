/*** Shift atoms that are not in the chosen list by linear interpolation **/
/*** The vector should contain the move vector as an x,y,z ordered list  **/
/*** for all 3 * num_atoms.                                              **/
/*** Routine introduced April 2019, Dave Willock, Ali Nasrallah.         **/
#include <stdio.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void interpolate_all_but_chosen(atom *p_molecule, atom *p_step_molecule, int num_atoms,
                                int *p_chosen, int num_chosen, double scale, double *p_vec )
{
atom *p_atom, *p_step_atom;

int imol, iatom, move_this;

int *p_this_chosen;


p_atom=p_molecule; p_step_atom=p_step_molecule;
for (iatom=0; iatom < num_atoms; iatom++)
  {
/************************************************************/
/**** Check if the atom is a member of the group defined ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
     move_this= TRUE; p_this_chosen=p_chosen;
     for (imol=0; imol<=num_chosen; imol++)
       {
          if ( *p_this_chosen == iatom ) move_this= FALSE;
          p_this_chosen++;
       }

     if ( move_this )
       {
          *p_step_atom = *p_atom;
          p_step_atom->x = p_step_atom->x + scale * *(p_vec+iatom*3);
          p_step_atom->y = p_step_atom->y + scale * *(p_vec+iatom*3+1);
          p_step_atom->z = p_step_atom->z + scale * *(p_vec+iatom*3+2);
       }
     p_atom++;p_step_atom++;
   }

return;
}
