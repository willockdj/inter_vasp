/***************************************************************************/
/*** Read in a ChemShell punch file to extract structure *******************/
/*** Dave Willock Nov. 2019 ************************************************/
/***************************************************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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

void cut_to_digit( char *p_word );

void find_line(FILE *fp, char *p_key, char *p_key2, int sep, int *p_found, int max_lines);

void find_line_with_stopper(FILE *fp, char *p_key, char *p_key2, char *p_stopper, int sep, int *p_found, int max_lines);

void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                            double *p_abc);

int read_force_block( FILE *fp, e_vec *p_forces, int num_atoms);
/*---------------------------------------------------------------------------*/

/* read in a outcar file for VASP movies */

int read_punch( FILE *fp, atom *p_molecule, atom *p_shells, int *p_num_shells, coord_flags *p_fix_flags, int just_count )
{
  int iloop, jloop, skip, iatom, noskip;

  int i, sep, have_num, found;
  int num_charges,num_atoms;

  char *p_key, *tok, *tok_nonum; 
  char *p_key_coords;
  char *p_key_charges;
  char *p_key_shells;

  atom *p_atom;

/* DEBUG */
  printf("Entered read_punch routine ");

  if (just_count) printf("just counting this time\n");
  else            printf("will read data this time\n");
  printf("Will see end of file from >>%s<<\n", END_OF_INPUT);
/* get all atom labels from POTCAR lines */
/* up to stopper                         */

  p_key= "block";
  p_key_coords= "coordinates";
  p_key_charges= "atom_charges";
  p_key_shells= "shells";

//  p_masskey = "POMASS";
//  p_masskey2 = "none";

  sep = 1;

  skip = TRUE;
  noskip = FALSE;
//  *p_num_types=-1;

  have_num = FALSE;
  tok= tok_get( fp, skip, FALSE);

/* Read structure expecting order coords, atom charges and shells */

  found=FALSE;

/*-----------------------------------------------*/
/*-- Jump to co-ordinates -----------------------*/
/*-----------------------------------------------*/
      find_line_with_stopper( fp, p_key, p_key_coords, p_key_charges, sep, &found, -1 );
          
      if (found)
        {
           for (iloop=0; iloop<3; iloop++) tok= tok_get( fp, noskip, FALSE);
           num_atoms = atoi(tok);
           printf("Found line next token is >>%s<< interpret as %d atoms\n", tok, num_atoms);

           if (!just_count)
             {

/* To get here must have just_count false */
/* so read in the atom list               */
/* convert Bohr to Angstrom units         */

                p_atom=p_molecule;
                for ( iatom=0; iatom < num_atoms; iatom++)
                  {
                    tok=tok_get(fp, skip, FALSE);
                    strcpy(p_atom->label, tok);
                    strcpy(p_atom->pot, "?"); 
                    strcpy(p_atom->group, "VASP"); 
                    strcpy(p_atom->group_no, "X"); 

                    strcpy(p_atom->elem,tok);
                    cut_to_digit(p_atom->elem);
      
                    p_atom->x = atof(tok_get(fp, noskip, FALSE))/BOHR;
                    p_atom->y = atof(tok_get(fp, noskip, FALSE))/BOHR;
                    p_atom->z = atof(tok_get(fp, noskip, FALSE))/BOHR;

                    p_fix_flags->fx =FALSE;
                    p_fix_flags->fy =FALSE;
                    p_fix_flags->fz =FALSE;

                    printf("Read Atom %d %s %10.6f  %10.6f  %10.6f\n", 
                                iatom, p_atom->label, p_atom->x, p_atom->y, p_atom->z);
                    p_atom++; p_fix_flags++;
                 }
             }
        }
      else
        {
           printf("ERROR: Failed to find number of atom coordinates in punch file\n");
           exit(0);
        }
/*-----------------------------------------------*/
/*-- Jump to charges ----------------------------*/
/*-----------------------------------------------*/
      found= FALSE;
      find_line_with_stopper( fp, p_key, p_key_charges, p_key_shells, sep, &found, -1 );
          
      if (found)
        {
           for (iloop=0; iloop<3; iloop++) tok= tok_get( fp, noskip, FALSE);
           num_charges = atoi(tok);

           if (num_charges != num_atoms)
             {
                printf("ERROR: Number of atoms (%d) and number of charges (%d) do not match in punch file\n",
                                                                                       num_atoms, num_charges);
                exit(0);
             }

           if (!just_count)
             {
                p_atom=p_molecule;
                for ( iatom=0; iatom < num_charges; iatom++)
                  {
                    p_atom->part_chge = atof(tok_get(fp, skip, FALSE));
                    p_atom++;
                  }
             }
        }
      else
        {
           printf("ERROR: Failed to find number of charges in punch file\n");
           exit(0);
        }
/*-----------------------------------------------*/
/*-- Jump to shells -----------------------------*/
/*-----------------------------------------------*/
      found= FALSE;
      find_line( fp, p_key, p_key_shells, sep, &found, -1 );

      if (found)
        {
           for (iloop=0; iloop<3; iloop++) tok= tok_get( fp, noskip, FALSE);
           *p_num_shells = atoi(tok);

           if (!just_count)
             {
                p_atom=p_shells;
                for ( iatom=0; iatom < *p_num_shells; iatom++)
                  {
                    tok=tok_get(fp, skip, FALSE);
                    strcpy(p_atom->label, tok);
                    strcpy(p_atom->pot, "?");
                    strcpy(p_atom->group, "VASP");
                    strcpy(p_atom->group_no, "X");

                    strcpy(p_atom->elem,tok);
                    cut_to_digit(p_atom->elem);

                    p_atom->x = atof(tok_get(fp, noskip, FALSE))/BOHR;
                    p_atom->y = atof(tok_get(fp, noskip, FALSE))/BOHR;
                    p_atom->z = atof(tok_get(fp, noskip, FALSE))/BOHR;
                    p_atom->part_chge = atof(tok_get(fp, noskip, FALSE));
 
/*** store final index as only Neighbour ***/
                    p_atom->num_neigh=1;
                    p_atom->neighb[0]=atoi(tok_get(fp, noskip, FALSE));

                    printf("Shell %d (%s) has index %d\n", iatom, p_atom->label, p_atom->neighb[0]);

                    p_atom++;
                  }
             }
        }
      else
        {
           printf("ERROR: Failed to find number of charges in punch file\n");
           exit(0);
        }

  printf("Returning from read_punch.c\n");
  return num_atoms;
}
