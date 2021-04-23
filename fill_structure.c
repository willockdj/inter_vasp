/**************************************************************************/
/*** Copy over atoms from master to molecule as indicated by the        ***/
/*** distrib array, Dave Willock July 2013                              ***/
/*** **********************************************************************/

#include <stdio.h>
#include "maxima.h"
#include "structures.h"

void fill_structure(atom *p_master, atom *p_molecule, int *p_indices, int *p_distrib, 
                    int num_indices, int *p_num_new_atoms, int num_atoms)
  {
    int iloop, index, iatom;
    int *p_this_index;
    atom *p_atom;

    *p_num_new_atoms=-1;

    for (iloop= 0; iloop <= num_atoms; iloop++)
      {
        if (iloop == *p_indices)
          {
             if (*p_distrib)
               {
                 *p_molecule= *p_master;
                 p_molecule++;
                 ++*p_num_new_atoms;
               }
             p_indices++; p_distrib++;
          }
        else
          {
             *p_molecule= *p_master;
             p_molecule++;
             ++*p_num_new_atoms;
          }
        p_master++;
      }
    --*p_num_new_atoms;

    return;
 }


