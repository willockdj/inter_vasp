/***************************************************************************/
/*** Read in a VASP OUTCAR *************************************************/
/*** Dave Willock & Ed Jeffery May 04 **************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

void cut_to_uscore( char *p_word );

void find_line(FILE *fp, char *p_key, char *p_key2, int sep, int *p_found, int max_lines);

void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                            double *p_abc);

int read_force_block( FILE *fp, e_vec *p_forces, int num_atoms);
/*---------------------------------------------------------------------------*/

/* read in a outcar file for VASP movies */

int read_siesta_vectors( FILE *fp, double *p_eigenvals, 
                         e_vec *p_eigenvecs, int *p_num_atoms, int num_free,
                         int *p_num_modes) 
{
  int iloop, jloop, skip=TRUE;
  int iatom;
  int noskip=FALSE;
  int num_to_skip;
  int itype, imode;
  int error, good_read;
  int have_all_labels;
  int sep, found, found_mode, endoffile;

  double *p_this_latt_vec;
  double *p_this_eigenval;
  double freq_sign;

  char *p_key, *p_key2, *p_key3;
  char *tok, *tok2, *p_letter;
  char *p_label;

  atom_number type_nums[10];
  atom *p_atom;

  e_vec *p_this_eigenvec;

  printf("Arrived in read_siesta_vectors.......\n");

  p_this_eigenval = p_eigenvals;
  p_this_eigenvec = p_eigenvecs;

  *p_num_modes = -1;

  endoffile=FALSE;
  while (!endoffile)
   {
     p_key= "Eigenvector";
     p_key2= "=";
     sep = 0;

     found=FALSE;
     find_line( fp, p_key, p_key2, sep, &found, (*p_num_atoms)+5 );

     if (found)
       {
         p_key= "Frequency";
         p_key2= "=";
         sep = 0;

         found=FALSE;
         find_line( fp, p_key, p_key2, sep, &found, 1);

/* Read in the frequency      */
         tok= tok_get( fp, noskip, FALSE);
         *p_this_eigenval= atof(tok);
         printf("Mode frequency %10.6f\n", *p_this_eigenval);

         (*p_num_modes)++;

         if (*p_num_modes > MAXMODES)
            {
              printf("ERROR: Number of modes in SIESTA vectors file ");
              printf("exceeds maximum of %d for inter_vasp\n", MAXMODES);
              exit(0);
            }
  
/*** Get corresponding eigenvector *****/
         p_key= "Eigenmode";
         p_key2= "(real";
         sep = 0;

         found_mode=FALSE;
         find_line( fp, p_key, p_key2, sep, &found_mode, -1 );

         for (iatom = 0; iatom <= *p_num_atoms; iatom++)
           {
             if ((num_free >= 0 && iatom < num_free) || num_free < 0)
               {
                 p_this_eigenvec->dx[iatom] = atof( tok_get( fp, skip, FALSE ));
                 p_this_eigenvec->dy[iatom] = atof( tok_get( fp, noskip, FALSE ));
                 p_this_eigenvec->dz[iatom] = atof( tok_get( fp, noskip, FALSE ));
               }
             else
               {
                 p_this_eigenvec->dx[iatom] = 0.0;
                 p_this_eigenvec->dy[iatom] = 0.0;
                 p_this_eigenvec->dz[iatom] = 0.0;
               }

                 printf("Atom %d %10.6f %10.6f %10.6f\n", iatom, 
                                                    p_this_eigenvec->dx[iatom],
                                                    p_this_eigenvec->dy[iatom],
                                                    p_this_eigenvec->dz[iatom] );
                                        
           }
        p_this_eigenval++;
        p_this_eigenvec++;
      }
     else
      {
        endoffile=TRUE;
      }
   }
  (*p_num_modes)++;
  return 0;
}
