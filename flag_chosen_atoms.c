#include <stdio.h>
#include <string.h>

/****************************************************************/
/** flag_chosen_atoms looks for the atoms selected for **********/
/** special treatment during interpolation and defines **********/
/** the rigid body interpolation vector as required    **********/
/****************************************************************/

#include "maxima.h"
#include "global_values.h"
#include "structures.h"

void centre_of_mass_flagged(double *p_c_of_m, double *p_total_mass, 
                            atom *p_molecule, 
                            int *p_chosen_indices, int num_chosen );

void move_selected_atoms(atom *p_molecule, int num_chosen, 
                         int *p_chosen_indices, double *move_vec);

void flag_chosen_atoms( atom *p_molecule, atom *p_end_molecule,
                        int num_atoms, int *p_chosen_indices,
                        int *p_num_chosen_atoms, int mol_ind, 
                        char *p_mol_cnt_lab,
                        double *p_start_cofm, 
                        double *p_end_cofm, double *p_inter_cofm)
{
int iatom, index, icomp;
int *p_ci, follow_centre=TRUE;
int ind_centre;

double vec[3], total_mass;

atom *p_atom, *p_end_atom;

*p_num_chosen_atoms= -1;

printf("Molecule %d is the rigid body\n", mol_ind);

p_ci=p_chosen_indices;
p_atom = p_molecule;
ind_centre=-1;
for (iatom= 0; iatom < num_atoms; iatom++)
  {
    printf("Checking %s\n", p_atom->label);
    if (p_atom->mol == mol_ind)
      {
        (*p_num_chosen_atoms)++;
        *p_ci= iatom;
/***********************************************/
/** Identify index supplied by user ************/
/***********************************************/
        if (follow_centre && strcmp(p_mol_cnt_lab, p_atom->label) == 0)
          {
            ind_centre = iatom;
          }

        p_ci++;
      } 
    p_atom++;
  }

printf("Number of chosen atoms=%d\n",*p_num_chosen_atoms);
p_ci=p_chosen_indices;
for (iatom= 0; iatom <= *p_num_chosen_atoms; iatom++) 
  {
   printf("%d %d\n",iatom,*(p_chosen_indices+iatom));
   p_ci++;
  }

printf("%d Chosen atoms for rigid body interpolation are:\n", 
                                                 *p_num_chosen_atoms);

p_ci=p_chosen_indices;
for (iatom= 0; iatom <= *p_num_chosen_atoms; iatom++)
  {
    p_atom     = p_molecule     + *p_ci;
    p_end_atom = p_end_molecule + *p_ci;
    printf("%d %s corresponding end atom %s", *p_ci, 
                                                p_atom->label,
                                                p_end_atom->label);
    if (follow_centre && *p_ci == ind_centre) 
      {
          printf("  the atom on which rigid body interpolation will be based\n");
      }
    else
      {
          printf("\n");
      }
    p_ci++;
  }

if (follow_centre)
  {
     *p_start_cofm     = (p_molecule+ind_centre)->x;
     *(p_start_cofm+1) = (p_molecule+ind_centre)->y;
     *(p_start_cofm+2) = (p_molecule+ind_centre)->z;
   
     *p_end_cofm     = (p_end_molecule+ind_centre)->x;
     *(p_end_cofm+1) = (p_end_molecule+ind_centre)->y;
     *(p_end_cofm+2) = (p_end_molecule+ind_centre)->z;
  }
else
  {
/**** get centres of mass at start and end ****/

    centre_of_mass_flagged(p_start_cofm, &total_mass, p_molecule,
                           p_chosen_indices, *p_num_chosen_atoms );

    centre_of_mass_flagged(p_end_cofm, &total_mass, p_end_molecule,
                           p_chosen_indices, *p_num_chosen_atoms);
  }

/**** inter_vec between cofms ******/

*p_inter_cofm     = *p_end_cofm     - *p_start_cofm;
*(p_inter_cofm+1) = *(p_end_cofm+1) - *(p_start_cofm+1);
*(p_inter_cofm+2) = *(p_end_cofm+2) - *(p_start_cofm+2);

printf("Centre of mass starts at %10.6f %10.6f %10.6f\n", 
                                           *p_start_cofm, 
                                           *(p_start_cofm+1), 
                                           *(p_start_cofm+2));

printf("Centre of mass ends at   %10.6f %10.6f %10.6f\n", 
                                           *p_end_cofm, 
                                           *(p_end_cofm+1), 
                                           *(p_end_cofm+2));

printf("Inter-centre vector      %10.6f %10.6f %10.6f\n", 
                                           *p_inter_cofm, 
                                           *(p_inter_cofm+1), 
                                           *(p_inter_cofm+2));

/**** show co-ordinates of chosen_ones relative to the c_of_m *****/

for (icomp=0; icomp<3; icomp++) vec[icomp]= -*(p_start_cofm+icomp); 

//move_selected_atoms(p_molecule, *p_num_chosen_atoms, 
//                                     p_chosen_indices, vec );


printf("Initial co-ordinates of selected molecule relative to c_of_m\n");

for (iatom= 0; iatom <= *p_num_chosen_atoms; iatom++)
  {
     index= *(p_chosen_indices+iatom);
     p_atom=p_molecule+index;
     printf("%s %10.6f %10.6f %10.6f\n", p_atom->label, 
                                         p_atom->x + vec[0],
                                         p_atom->y + vec[1],
                                         p_atom->z + vec[2]);
  }

for (icomp=0; icomp<3; icomp++) vec[icomp]= -*(p_end_cofm+icomp); 

//move_selected_atoms(p_end_molecule, *p_num_chosen_atoms, 
//                                          p_chosen_indices, vec );

printf("Final co-ordinates of selected molecule relative to c_of_m\n");

for (iatom= 0; iatom <= *p_num_chosen_atoms; iatom++)
  {
     index= *(p_chosen_indices+iatom);
     p_atom=p_end_molecule+index;
     printf("%s %10.6f %10.6f %10.6f\n", p_atom->label, 
                                         p_atom->x + vec[0],
                                         p_atom->y + vec[1],
                                         p_atom->z + vec[2]);
  }

return;
}
