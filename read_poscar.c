/***************************************************************************/
/*** Read in a VASP POSCAR and related files *******************************/
/*** Dave Willock June 03 **************************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

int read_line(FILE *fp, int *p_ichar);

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int read_atom_data_vasp( FILE *fp, atom *p_atom, int *p_mol_number, coord_flags *p_fix_flags);
                  

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

double get_double(FILE *input_fp, int skip_lines, int *p_error);

int get_integer( FILE *input_fp, int skip_lines, int *p_error );

void put_string (FILE *fp, int *p_ichar, int length);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                            double *p_abc );
/*---------------------------------------------------------------------------*/

/* read in a poscar vasp formatted file */

int read_poscar( FILE *fp, int *p_title_line, 
                 atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms, 
                 int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number, 
                 double *p_latt_vec, double *p_recip_latt_vec, double *p_abc, int *p_been_before, 
		 int *p_is_fract, int *p_is_cart, int *p_ion_number, int *p_num_types, 
                 double *p_scale_factor, coord_flags *p_fix_flags)

{
  int ichar[LINESIZ];
  int place,itsanum;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int at_end, skip, lower_case;
  int start,start2,start3,start4;
  int error, good_read;
  int *p_this_ion_number;

  double dot, cell_volume;
  double a_cross_b[3], b_cross_c[3], c_cross_a[3];

  char *p_key;
  char *tok;

  atom *p_atom;

/* assume poscar file only ever contains one molecule */

     *p_num_of_mols=1;

/* get title line */

     num_of_chars= read_line ( fp, p_title_line);
     printf("Read title line: \n");
     put_string(stdout, p_title_line, num_of_chars);
     printf("\n");

/* get scale factor */

     skip= TRUE;
     lower_case= TRUE;
     
     *p_scale_factor = get_double( fp, skip, &error);

     if (!error)
       printf("Scale factor %10.6f\n", *p_scale_factor);
     else
       printf("Error while reading scale factor from POSCAR\n");

/* VASP POSCAR is always periodic, work out cell vectors and angles */
/* VASP POSCAR has lattice vectors in array as:   ax ay az          */
/* VASP POSCAR                                    bx by bz          */
/* VASP POSCAR                                    cx cy cz          */

     skip= TRUE;
     *p_latt_vec        = get_double(fp, skip, &error);
     skip= FALSE;
     *(p_latt_vec+1)    = get_double(fp, skip, &error);
     *(p_latt_vec+2)    = get_double(fp, skip, &error);

     skip= TRUE;
     *(p_latt_vec+3)    = get_double(fp, skip, &error);
     skip= FALSE;
     *(p_latt_vec+4)    = get_double(fp, skip, &error);
     *(p_latt_vec+5)    = get_double(fp, skip, &error);

     skip= TRUE;
     *(p_latt_vec+6)    = get_double(fp, skip, &error);
     skip= FALSE;
     *(p_latt_vec+7)    = get_double(fp, skip, &error);
     *(p_latt_vec+8)    = get_double(fp, skip, &error);

    latt_vecs_from_cart( p_latt_vec, p_recip_latt_vec, p_abc );

/* Read in the ion number flags */

      skip = TRUE;
	      
      error = FALSE;
      *p_num_types=0;
      p_this_ion_number =p_ion_number;
      while (!error)
	{
          *p_this_ion_number=get_integer( fp, skip, &error );

          skip = FALSE;
          ++*p_num_types;
          if (*p_num_types < MAXTYPES)
            {
               p_this_ion_number++;
            }
          else
            {
               printf("ERROR: MAXTYPES exceeded while reading VASP POSCAR file\n");
               exit(0);
            }
        }

      --*p_num_types;

      printf("Have %d types of ion with ",*p_num_types);

      *p_num_atoms=0;
      p_this_ion_number =p_ion_number;
      for (iloop=0; iloop < *p_num_types; iloop++)
	{
	   printf("%d, ", *p_this_ion_number);
	   *p_num_atoms += *p_this_ion_number;
           p_this_ion_number++;
	}

      printf(" instances of each, totaling %d\n", *p_num_atoms);

/* Read in the direct/cartesian flag data */

    skip=TRUE;
    tok = tok_get( fp, skip, TRUE);
    printf("Next read >>%s<<\n",tok);

/* Ignore Selective Dynamics line if present */

    if (!strncmp(tok,"sel",3))
      {
        tok = tok_get( fp, skip, TRUE);
        printf("Trying to skip Selective Dynamics, read >>%s<<\n",tok);
      }

    if (!strcmp(tok,"direct"))
    {
       printf("We have fractionals\n");
       *p_is_fract = TRUE;
    }
    else
    {
       printf("We have cartesians\n");
       *p_is_fract = FALSE;
    }

/* Read in the atomic data */

      for (iloop=0; iloop < *p_num_atoms; iloop++)
       {
        good_read= read_atom_data_vasp( fp, p_molecule+iloop, p_mol_number,  
                                                  p_fix_flags );
        p_fix_flags++;
       }

      printf("Co-ordinates read in:\n");
      for (iloop=0; iloop < *p_num_atoms; iloop++)
	      printf("%10.6f %10.6f %10.6f\n", (p_molecule+iloop)->x,
			                       (p_molecule+iloop)->y,
					       (p_molecule+iloop)->z);

  return 0;
}
