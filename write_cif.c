/**************************************************************************/
/*** Write a VASP POSCAR file *********************************************/
/*** Code assumes that the num_atoms is the upper index  ******************/
/*** Last update April 2016 Dave Willock                 ******************/
/**************************************************************************/

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"

/* ------Prototype-list---------------------------------------- */

void put_string(FILE *fp, int *p_ichar, int length);

/* ------------------------------------------------------------ */

void write_cif( FILE *fp, atom *p_molecule, double *p_fract_coords,
                atom_number *p_types, int num_types,
                double *p_latt_vec, double *p_abc, double *p_scale_factor, int num_atoms,
                int *p_title_line, char *p_c_title_line, int pbc, int is_fract,
                coord_flags *p_fix_flags, int need_zsort)
{
   int iloop, jloop, this_atom;
   int itype, iatom, ind, start_list;

   double dummy, *p_frac1, *p_frac2;

   coord_flags *p_this_flag;

   atom_number *p_this_type;

   atom temp_atom, *p_start_atom, *p_atom;

/* Check to see if a title was given in the structure file */
/* If not use the title from the input file                */

   printf("Writing cif file.....\n");
 
   if ( *p_title_line > 0)
     {
       fprintf(fp, "\#");
       put_string(fp, p_title_line, 100);
     }
   else
     {
       fprintf(fp, "\#%s\n", p_c_title_line);
     }
   printf("Writing cif.....1.5..\n");

/* periodic boundaries must be required if user needs POSCAR */

   if (!pbc) 
    {
	    printf("ERROR: cif file requested for non-periodic input\n");
	    exit(0);
    }

/* print out lattice vectors */
   fprintf(fp,"data_intervasp\n");
   fprintf(fp,"_symmetry_space_group_name_Hall  \'P1\'\n");
   fprintf(fp,"_symmetry_space_group_name_H-M   \'P1\'\n");
   fprintf(fp,"_cell_angle_alpha                %10.6f\n",*(p_abc+3));
   fprintf(fp,"_cell_angle_beta                 %10.6f\n",*(p_abc+4));
   fprintf(fp,"_cell_angle_gamma                %10.6f\n",*(p_abc+5));
   fprintf(fp,"_cell_length_a                   %10.6f\n",*p_abc);
   fprintf(fp,"_cell_length_b                   %10.6f\n",*(p_abc+1));
   fprintf(fp,"_cell_length_c                   %10.6f\n",*(p_abc+2));


/**********************************************************************/
/* Atoms should already be ordered according to type                ***/
/* So use types information to generate number of each element line ***/
/**********************************************************************/
  
   printf("Writing cif_file.3...\n");
  
/**********************************************************************/
/*** If needed re-order atoms by z, this makes fixing slab layers *****/
/*** easier. Dave Willock, May 2018...............................*****/
/**********************************************************************/

if (need_zsort)
  {
    printf("Sorting each type in z-co-ordinate order....\n");

    p_this_type= p_types; p_start_atom=p_molecule; start_list=0;
    for (itype=0; itype <= num_types; itype++) 
      {
        printf("\nreordering atoms of type %s there are %d of these\n", p_this_type->atom_type, p_this_type->num);

        for (this_atom=0; this_atom <= p_this_type->num; this_atom++)
          {
             for (iatom=this_atom+1; iatom<= p_this_type->num; iatom++)
               {
                  if ( (p_start_atom+this_atom)->z >  (p_start_atom+iatom)->z ) 
                    {
                      temp_atom                = *(p_start_atom+this_atom); 
                      *(p_start_atom+this_atom)= *(p_start_atom+iatom);
                      *(p_start_atom+iatom)    = temp_atom;
//
// Re-order fractionals too
//
                      p_frac1= p_fract_coords + 3*(start_list+this_atom);
                      p_frac2= p_fract_coords + 3*(start_list+iatom);
                 
                      for (ind=0; ind<3; ind++)
                        {
                           dummy= *p_frac1;
                           *p_frac1 = *p_frac2;
                           *p_frac2 = dummy;

                           p_frac1++; p_frac2++;
                        } 
                    } 
               }
          }
        p_start_atom+=(p_this_type->num) + 1; 
        start_list+=(p_this_type->num) + 1; 
        p_this_type++;
      }
  }
/**********************************************************************/
/* Print out atom co-ordinates ****************************************/
/**********************************************************************/

 if (is_fract)
  {
   printf("Writing out %d atoms using fractional co-oridinates...\n", num_atoms);
   fprintf(fp,"loop_\n");
   fprintf(fp,"_atom_site_label\n");
   fprintf(fp,"_atom_site_type_symbol\n");
   fprintf(fp,"_atom_site_fract_x\n");
   fprintf(fp,"_atom_site_fract_y\n");
   fprintf(fp,"_atom_site_fract_z\n");
#
   p_atom=p_molecule; 
   for (this_atom=0; this_atom < num_atoms; this_atom++)
     {
       fprintf(fp,"%s  %s   %14.10f  %14.10f  %14.10f\n", 
                 p_atom->label, p_atom->elem, *p_fract_coords, *(p_fract_coords+1), *(p_fract_coords+2));
       p_fract_coords += 3;
       p_atom++;
     }
  }
else
  {
   printf("ERROR: cif files must write fractional co-ordinates, but is_fract not set\n");
  }

return;
}
