/***************************************************************************/
/*** Read in a VASP POTCAR and related files *******************************/
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

/*---------------------------------------------------------------------------*/


/* read in a potcar file for VASP labels */

int read_potcar( FILE *fp, atom *p_molecule,
		 atom_number *p_types, int *p_num_types)

{
  int iloop, skip, iatom;
  int error, good_read;

  char *p_key;
  char *tok, *p_letter;
  char *p_label;

  atom *p_atom;

/* get title lines */

  p_key= "VRHFIN";

  skip = TRUE;

  for ( iloop=0; iloop <= *p_num_types; iloop++) 
     {
        tok= tok_get( fp, skip, FALSE);

/****************************************************/
/*** Ignore all lines except TITEL ******************/
/*** !tok catches NULL !           ******************/
/****************************************************/
        while ( !tok || strcmp(tok,p_key) != 0 ) tok= tok_get( fp, skip, FALSE);
/***************************************/
/*** Debug for VASP64 ******************/
/*** Assume element is second object ***/
/***************************************/
         tok= tok_get( fp, FALSE, FALSE); 
         p_label = tok; 

/**** remove leading equals sign if present ****/
  
         if (*p_label == '=') 
           {
             printf("Found equals sign in label moving pointer\n");
             p_label++;
             printf("label=>>%s<<\n",p_label);
           }

         p_letter=p_label;
         while ( *p_letter != '\0' ) 
            {
               p_letter++;
               if ( *p_letter == '_' || *p_letter == ' ' || *p_letter == ':') *p_letter = '\0'; 
            }

/****************************************************/
/*** Put labels on the atom co-ordinate *************/
/****************************************************/
	for (iatom=0; iatom < p_types->num; iatom++)
           {
               strcpy(p_molecule->elem, p_label); 
	       sprintf(p_molecule->label,"%s%d",p_label,iatom+1);
               strcpy(p_molecule->pot, "?"); 
               strcpy(p_molecule->group, "VASP"); 
               printf("Assigned label %s and element %s to atom %d based on POTCAR information\n", 
                                 p_molecule->label, p_molecule->elem, iatom);
	       p_molecule++;
           }
        p_types++;
    }

  return 0;
}
