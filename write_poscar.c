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

void write_poscar( FILE *fp, atom *p_molecule, double *p_fract_coords,
                    atom_number *p_types, int num_types,
                   double *p_latt_vec, double *p_scale_factor, int num_atoms,
                   int *p_title_line, char *p_c_title_line, int pbc, int is_fract,
                   coord_flags *p_fix_flags, int need_zsort)
{
   int iloop, jloop, this_atom;
   int itype, iatom, ind, start_list;

   double dummy, *p_frac1, *p_frac2;

   coord_flags temp_flag;

   atom_number *p_this_type;

   atom temp_atom, *p_start_atom;

/* Check to see if a title was given in the structure file */
/* If not use the title from the input file                */

   printf("Writing POSCAR.......\n");
 
   if ( *p_title_line > 0)
     {
       put_string(fp, p_title_line, 100);
     }
   else
     {
       fprintf(fp, "%s\n", p_c_title_line);
     }
   printf("Writing POSCAR..1.5..\n");

/* periodic boundaries must be required if user needs POSCAR */

   if (!pbc) 
    {
	    printf("ERROR: VASP POSCAR file requested for non-periodic input\n");
	    exit(0);
    }

/* make scale factor 1 if none available */

  if (*p_scale_factor < 0.0)
    {
      fprintf( fp,"1.0\n"); 
    }
  else
    {
      fprintf( fp,"%10.6f\n",*p_scale_factor);
    }

/* print out lattice vectors */

  for (iloop=0; iloop < 3; iloop++)
    {
      for (jloop=0; jloop < 3; jloop++)
        {
          fprintf(fp,"%14.10f ",*p_latt_vec);
          p_latt_vec++;
        }
      fprintf(fp,"\n");
    }

/**********************************************************************/
/* Atoms should already be ordered according to type                ***/
/* So use types information to generate number of each element line ***/
/**********************************************************************/
  
p_this_type= p_types;
for (itype=0; itype <= num_types; itype++) 
   {
     fprintf( fp,"%s ", p_this_type->atom_type); 
     p_this_type++;
   }
fprintf( fp,"\n");

p_this_type= p_types;
for (itype=0; itype <= num_types; itype++) 
   {
     fprintf( fp,"%d ", (p_this_type->num) + 1 ); 
     p_this_type++;
   }
fprintf( fp,"\n");

fprintf(fp,"Selective dynamics\n");
 
if (is_fract)
  {
    fprintf( fp,"Direct\n");
  }
else
  {
    fprintf( fp,"cartesian\n");
  }

   printf("Writing POSCAR...3...\n");
  
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
//
// Re-order flags
//
                      printf("Shifting flags, start_list %d this atom %d iatom %d\n", start_list, this_atom, iatom);
                      temp_flag                           = *(p_fix_flags+start_list+this_atom);
                      *(p_fix_flags+start_list+this_atom) = *(p_fix_flags+start_list+iatom);
                      *(p_fix_flags+start_list+iatom)     = temp_flag;
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
   for (this_atom=0; this_atom < num_atoms; this_atom++)
   {
     fprintf(fp,"%14.10f  %14.10f  %14.10f", 
               *p_fract_coords, *(p_fract_coords+1), *(p_fract_coords+2));
     p_fract_coords += 3;

     if (p_fix_flags->fx)
       {
         fprintf(fp, " F");
       }
     else
       {
         fprintf(fp, " T");
       }
     if (p_fix_flags->fy)
       {
         fprintf(fp, " F");
       }
     else
       {
         fprintf(fp, " T");
       }
     if (p_fix_flags->fz)
       {
         fprintf(fp, " F\n");
       }
     else
       {
         fprintf(fp, " T\n");
       }
     p_fix_flags++;
   }
  }
else
  {
   printf("Writing out %d atoms...\n", num_atoms);
   for (this_atom=0; this_atom < num_atoms; this_atom++)
   {
     fprintf(fp,"%14.10f  %14.10f  %14.10f", 
                        p_molecule->x, p_molecule->y, p_molecule->z);

/**********************************************************/
/*** In my routines fix_flags will be TRUE if it should  **/
/*** be fixed so requiring and F flag in the POSCAR file **/
/**********************************************************/
     if (p_fix_flags->fx)
       {
         fprintf(fp, " F");
       }
     else
       {
         fprintf(fp, " T");
       }
     if (p_fix_flags->fy)
       {
         fprintf(fp, " F");
       }
     else
       {
         fprintf(fp, " T");
       }
     if (p_fix_flags->fz)
       {
         fprintf(fp, " F\n");
       }
     else
       {
         fprintf(fp, " T\n");
       }
     p_molecule++;
     p_fix_flags++;
   }
  }

return;
}
