/***************************************************************************/
/*** Read in a VASP XDATCAR file for MD trajectories ***********************/
/*** Dave Willock May 06 ***************************************************/
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

int read_atom_data_vasp( FILE *fp, atom *p_atom, int *p_mol_number,
                          coord_flags *p_fix_flags );

double get_double(FILE *input_fp, int skip_lines, int *p_error);

int get_integer( FILE *input_fp, int skip_lines, int *p_error );

void put_string (FILE *fp, int *p_ichar, int length);

/*---------------------------------------------------------------------------*/

/* read in a xdatcar vasp formatted file */

int read_xdatcar( FILE *fp, int *p_num_frames, atom *p_molecule, int num_atoms,
                  coord_flags *p_fix_flags )
{
  int ichar[LINESIZ];
  int iatom,idummy;
  int skip, lower_case=FALSE;
  int error;

  char *p_key;
  char *tok;

/* assume xdatcar file is just a read after atoms have been sorted out */

   printf("Trying to read XDATCAR file format\n");

   if (*p_num_frames < 0)
     {
/* get number of frames */
        skip = TRUE;
        idummy = get_integer( fp, skip, &error ); 
        printf("%d\n",idummy);
        skip = FALSE;
        idummy = get_integer( fp, skip, &error ); 
        printf("%d\n",idummy);
        *p_num_frames = get_integer( fp, skip, &error ); 
        printf("%d\n",*p_num_frames);

        if (error)
          {
             printf("ERROR: Bad file format on first line of XDATCAR\n");
             exit(0);
          }
        else
          {
/*** Jump lines to set up for first frame read next time in here ***/
             skip=TRUE;
             tok= tok_get(fp, skip, lower_case);
             printf("TOK: %s\n",tok);
             tok= tok_get(fp, skip, lower_case);
             printf("TOK: %s\n",tok);
             tok= tok_get(fp, skip, lower_case);
             printf("TOK: %s\n",tok);
             tok= tok_get(fp, skip, lower_case); 
             printf("TOK: %s\n",tok);

             return 0;
          }
    }
  else
    {
 /**** Must already have number of frames and am now just reading one in ****/
      skip = TRUE;
      printf("Back in read_xdatcar for another frame\n");
 /**** First skip blank line ***/
/*      while (!(tok= tok_get(fp, skip, lower_case))); */
/*      printf("TOK: %s\n",tok);  */

      for (iatom=0; iatom < num_atoms; iatom++)
        {
          read_atom_data_vasp( fp, p_molecule, &idummy, p_fix_flags );

          printf("Read: %10.6f  %10.6f  %10.6f\n", p_molecule->x,
                                                   p_molecule->y,
                                                   p_molecule->z);
          p_molecule++;
        }
    }

   printf("Frame read back to main.....\n") ;

  return 0;
}
