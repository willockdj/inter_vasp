/***************************************************************************/
/*** Read in a VASP INCAR **************************************************/
/*** Dave Willock **********************************************************/
/***                                  **************************************/
/*** Get information from the INCAR file supplied. Currently just look   ***/
/*** for and read in the MAGMON values that are set.                     ***/
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

void find_line(FILE *fp, char *p_key, char *p_key2, int sep, int *p_found, int max_lines);

void find_line_with_stopper(FILE *fp, char *p_key, char *p_key2, char *p_stopper, int sep, int *p_found, int max_lines);

int read_star_format( char *p_word, char *p_after );

int read_incar( FILE *fp, double *p_magmom, int *p_num_magmom )
{
  int iloop, skip;
  int noskip=FALSE;
  int num_to_fill;
  int itype, imode;
  int have_magmom, have_star;
  int sep, found;

  char *p_key, *p_key2, *p_key3, *p_stopper;
  char *tok;
  char after[10];

  printf("Arrived in read_incar\n");

  p_key= "MAGMOM";
  p_key2= "none";
  sep = 0;
  skip = TRUE;
  *p_num_magmom=-1;

  have_magmom = FALSE;
  tok= tok_get( fp, skip, FALSE);

/* Read through the file until you find the MAGMON keyword   */

  *p_num_magmom = -1;
  printf("Starting while loop...\n");
  found=FALSE;
  printf("finding line\n");
  find_line( fp, p_key, p_key2, sep, &found, -1 );

  if (found)
    {
       have_magmom=TRUE;
       printf("Found MAGMON line in INCAR file\n");
       while ((tok= tok_get( fp, FALSE, FALSE)) != NULL) 
         {
            printf("Read: %s ", tok);
            if (strcmp(tok, "=") != 0)
              {
                printf("...copying\n");

                have_star = read_star_format( tok, &after[0] );
                
                if (have_star) 
                 {
                   num_to_fill=atoi(tok);
                   printf("star format, filling %d elements with >>%s<<\n", num_to_fill, after);
                   for (iloop=0; iloop<num_to_fill; iloop++)
                     {
                       *p_magmom=atof(after);
                       p_magmom++;
                       ++(*p_num_magmom);
                     }
                 }
                else
                 { 
                   *p_magmom=atof(tok);
                    p_magmom++;
                    ++(*p_num_magmom);
                 }
              }
            printf("\n");
         }
    }

   printf("DEBUG>> Leaving read_incar\n");
  return 0;
}
