
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
                   coord_flags *p_fix_flags)
{
   int iloop, jloop, mol_current, this_atom;
   int itype, havefixed;

   coord_flags *p_this_flag;

/* Check to see if a title was given in the structure file */
/* If not use the title from the input file                */
 
   if ( *p_title_line > 0)
     {
       put_string(fp, p_title_line, 100);
    /*   fprintf(fp, "\n"); */
     }
   else
     {
       fprintf(fp, "%s\n", p_c_title_line);
     }

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
  
for (itype=0; itype <= num_types; itype++) 
   {
     fprintf( fp,"%d ", (p_types->num) + 1 ); 
     p_types++;
   }
fprintf( fp,"\n");

havefixed=FALSE;
p_this_flag= p_fix_flags;
for (this_atom=0; this_atom < num_atoms; this_atom++)
  {
     if (p_this_flag->fx ||  p_this_flag->fy ||  p_this_flag->fz )  
       {
          havefixed=TRUE;
          break;
       }
  }

if (havefixed)
  {
    fprintf(fp,"Selective dynamics\n");
  }
 
if (is_fract)
  {
    fprintf( fp,"direct\n");
  }
else
  {
    fprintf( fp,"cartesian\n");
  }

  
/**********************************************************************/
/* Print out atom co-ordinates ****************************************/
/**********************************************************************/

 if (is_fract)
  {
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
