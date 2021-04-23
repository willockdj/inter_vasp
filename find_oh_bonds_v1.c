/**********************************************************************/
/* find_OH_bonds looks through an atom list, identifies and measures  */
/* H-bonds involving the O atoms in pre-compiled lists.               */
/* started Dave, Mala, Sachin, Connie, May 2019                       */
/**********************************************************************/
/*** Pre-filled arrays will be: ***************************************/
/*** num_oxy_species will hold the number of each oxygen species found ***/
/*** num_oxy_species[0] is for O lattice ***/
/*** num_oxy_species[1] is for OH        ***/
/*** num_oxy_species[2] is for water O   ***/
/*** num_oxy_species[3] is for O alcohol ***/
/******************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "global_values.h"
#include "maxima.h"
#include "structures.h"

double atom_separation_squared(atom *p_A, atom *p_B, int pbc,  double *p_recip_latt_vec, double *p_latt_vec);

void find_oh_bonds( atom *p_molecule, int num_atoms, int *oxy_list_ptrs[],
                    atom_number *p_num_oxy_species, int num_oxy_types, 
                    int use_pbc, double *p_latt_vec, double *p_recip_latt_vec)
{
int iatom1,iatom2;
int ioxy, itype;

double standard,actual_dist;

atom *p_atom1, *p_atom2;

atom_number *p_this_type;

printf("Arrived in find_oh_bonds\n\n");
/**********************************************************************/
/******************* loop over types of oxygen we have    *************/
/**********************************************************************/
for (itype = 0; itype < num_oxy_types; itype++)
  {
    printf("List of oxygens that are: %s\n", p_num_oxy_species->atom_type); 
    for (ioxy=0; ioxy < p_num_oxy_species->num; ioxy++)
      { 
        p_atom1= p_molecule+oxy_list_ptrs[itype][ioxy];
        printf("%d %s which has %d neighbours.\n", ioxy, p_atom1->label, p_atom1->num_neigh );
      }
    
    p_num_oxy_species++;
  }

return;
}
