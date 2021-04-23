/* routine to decipher a VASP POSCAR atom data line */

#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"
#include "reader.h"

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

/* ------------------------------------------------------------ */

double get_double(FILE *input_fp, int skip_lines, int *p_error);

/* ------------------------------------------------------------ */

int read_atom_data_vasp( FILE *fp, atom *p_atom, int *p_mol_number,
                          coord_flags *p_fix_flags )

{
int good_line, i, skip, no_flags_here;

char cfx[3], cfy[3], cfz[3];
char *tok;

/* get co-ordinates */

       fscanf(fp, "%le %le %le", &(p_atom->x),  &(p_atom->y),  &(p_atom->z));

       printf("DEBUG>> Read co-ords: %20.16f  %20.16f  %20.16f \n",
                               p_atom->x, p_atom->y, p_atom->z);

/*** Allow skip to read remainder of line into buffer within tok_get ***/
       skip= TRUE;
       no_flags_here=FALSE;

       tok = tok_get( fp, skip, TRUE);
       if (tok) strcpy(cfx,tok);
       else 
        {
          no_flags_here=TRUE; 
          strcpy(cfx,"T");
          strcpy(cfy,"T");
          strcpy(cfz,"T");
        }
       
       if (!no_flags_here)
        {
           skip= FALSE;
           tok = tok_get( fp, skip, TRUE);
//           printf("Read....y.....tok.....>>%s<<\n", tok);
           if (tok) strcpy(cfy,tok);
                        else strcpy(cfy,"T");
       
           tok = tok_get( fp, skip, TRUE);
//           printf("Read....z.....tok.....>>%s<<\n", tok);
           if (tok) strcpy(cfz,tok);
                        else strcpy(cfz,"T");

//           printf("Found flags : %s %s %s\n", cfx, cfy, cfz );
            
        }

       p_fix_flags->fx = FALSE;
       p_fix_flags->fy = FALSE;
       p_fix_flags->fz = FALSE;

       for (i = 0; cfx[i] != '\0'; i++) cfx[i] = tolower(cfx[i]);
       for (i = 0; cfy[i] != '\0'; i++) cfy[i] = tolower(cfy[i]);
       for (i = 0; cfz[i] != '\0'; i++) cfz[i] = tolower(cfz[i]);

//       printf("DEBUG>> Read flags: %s %s %s\n", cfx, cfy, cfz);

       if (strcmp(cfx,"f") == 0) p_fix_flags->fx = TRUE;
       if (strcmp(cfy,"f") == 0) p_fix_flags->fy = TRUE;
       if (strcmp(cfz,"f") == 0) p_fix_flags->fz = TRUE;

/***************** Set all atom data to dummy values ******/

   p_atom->part_chge  = 0.0;

/* Fill in all extra data required */

   *p_mol_number= 0;


   return good_line;
}
