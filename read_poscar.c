/***************************************************************************/
/*** Read in a VASP POSCAR and related files *******************************/
/*** Dave Willock June 03 **************************************************/
/*** Dave Willock Last updated May 2016 ************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
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
                 int *p_have_labels, int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number, 
                 double *p_latt_vec, double *p_recip_latt_vec, double *p_abc, int *p_been_before, 
		 int *p_is_fract, int *p_is_cart, atom_number *p_types, int *p_num_types, 
                 double *p_scale_factor, coord_flags *p_fix_flags, int just_count)

{
  int ichar[LINESIZ];
  int place,itsanum, ntyp;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int at_end, skip, done, lower_case;
  int start,start2,start3,start4;
  int error, good_read;

  char *p_key;
  char *tok;

  atom *p_atom;

  atom_number *p_this_type;

  labels *p_this_label;

/* assume poscar file only ever contains one molecule */

     *p_num_of_mols=0;

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
      tok = tok_get( fp, skip, FALSE);
      printf("Next read >>%s<<\n",tok);
   
      if (isalpha(*tok))
        {
          printf("This is a new POSCAR file...reading labels\n");
/*** Read labels until you run out ****/
          *p_have_labels = TRUE;
  
          skip = FALSE;
          done = FALSE;
          *p_num_types=-1; p_this_type = p_types;

          while ( !done )
            {
              printf("tok >>%s<<\n", tok);
              strcpy( (p_this_type->atom_type), tok);
              p_this_type++;
              ++*p_num_types;

/*** Need to save the labels as types ***/

              tok = tok_get( fp, skip, FALSE);
              done = tok == NULL;
            }
        }
      else
        {
          printf("This is an old POSCAR file....will have to get labels elsewhere\n");
          *p_have_labels = FALSE;
/*** Rewind file and skip back down to the number line ***/
          rewind(fp);
          for ( iloop=1; iloop < 6; iloop++) tok = tok_get( fp, skip, FALSE);
        }
	      
          error = FALSE;
          *p_num_types=-1;
          p_this_type =p_types;
          while (!error)
            {
              p_this_type->num=get_integer( fp, skip, &error );

              skip = FALSE;
              ++*p_num_types;
              if (*p_num_types < MAXTYPES)
                {
                   p_this_type++;
                }
              else
                {
                   printf("ERROR: MAXTYPES exceeded while reading VASP POSCAR file\n");
                   exit(0);
                }
            }
          --*p_num_types;

      printf("Have %d types of ion...\n",*p_num_types);
      if (!*p_have_labels) printf(".....but not know what they are yet...\n");

      *p_num_atoms=0;
      p_this_type =p_types;
      for (iloop=0; iloop <= *p_num_types; iloop++)
	{
           if (*p_have_labels)
             {
	        printf("%s : %d, ", p_this_type->atom_type, p_this_type->num);
             }
           else
             {
	        printf("%d, ", p_this_type->num);
             }
	   *p_num_atoms += p_this_type->num;
           p_this_type++;
	}

      printf("\nTotaling %d\n", *p_num_atoms);

/* Read in the direct/cartesian flag data */

    skip=TRUE;
    tok = tok_get( fp, skip, TRUE);

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

   if (!just_count)
     {
       p_this_type = p_types;
       p_atom=p_molecule; ntyp=1;

       for (iloop=0; iloop < *p_num_atoms; iloop++)
        {
         good_read= read_atom_data_vasp( fp, p_atom, p_mol_number,  
                                                   p_fix_flags );
         strcpy(p_atom->pot, "?");
         strcpy(p_atom->group, "VASP");
         strcpy(p_atom->group_no, "X");
  
/*** If we have the labels put them on now ***/
         if (*p_have_labels)
           {
             strcpy(p_atom->elem, p_this_type->atom_type);  
             sprintf(p_atom->label,"%s%d",p_this_type->atom_type,ntyp);

             if (ntyp == p_this_type->num)
               {
                 ntyp=0;
                 p_this_type++;
               }
           } 

         ntyp++; p_fix_flags++; p_atom++;
        }

       printf("Co-ordinates read in:\n");
       if (*p_have_labels)
         {
           for (iloop=0; iloop < *p_num_atoms; iloop++)
 	            printf("%s (%s) %10.6f %10.6f %10.6f\n", (p_molecule+iloop)->label,
                                                             (p_molecule+iloop)->elem,
                                                             (p_molecule+iloop)->x,
 			                                     (p_molecule+iloop)->y,
                                                             (p_molecule+iloop)->z);
         }
       else
         {
           for (iloop=0; iloop < *p_num_atoms; iloop++)
 	            printf("%10.6f %10.6f %10.6f\n", (p_molecule+iloop)->x,
 			                             (p_molecule+iloop)->y,
                                                     (p_molecule+iloop)->z);
         }
     }

  return 0;
}
