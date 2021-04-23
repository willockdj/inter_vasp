/**********************************************************************/
/**************** WRITES OUT DATA TO ONETEP FILE **********************/
/**************** CHRISTIAN REECE DAVID WILLOCK ***********************/
/**************** SEPT 2015 *******************************************/
/**********************************************************************/

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "constants.h"
#include "structures.h"
#include "reader.h"

/* ------Prototype-list---------------------------------------- */

void put_string(FILE *fp, int *p_ichar, int length);
char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

/* ------------------------------------------------------------ */

void write_onetep( FILE *fp, atom *p_molecule, double *p_fract_coords,
                    atom_number *p_types, int num_types,
                   double *p_latt_vec, double *p_scale_factor, int num_atoms,
                   int *p_title_line, char *p_c_title_line, int pbc, int is_fract, 
                   int is_onetep, char *file_line_ptrs[], int num_lines)
{
   int iloop, jloop, mol_current, this_atom;
   int itype, havefixed; 
   int token, token2, skip, lower_case;
   char tok[BUFFER], *tok2;

   double *p_this_latt_vec;

   coord_flags *p_this_flag;

if (is_onetep)
  {
    printf("Data is from a onetep.dat file.. this will be much easier %d lines\n", num_lines);
    for (iloop = 0; iloop <= num_lines; iloop++)
      {
      printf("Taking line %d >>%s<<\n", iloop, file_line_ptrs[iloop]);

      strcpy(tok, file_line_ptrs[iloop]);
      tok2 = strtok(tok," ,.-");

      if (tok2 == NULL) 
        {
          printf("Whoops....\n");
          exit(0);
        }

      printf("current token is >>%s<< for line %d\n",tok2, iloop);
      token = find_kind(tok, ONETEP_DIRECTIVES);
      switch (token)
        {
          case ENDBLOCK   : printf("found Block\n");
                          tok2 = strtok(NULL," ,.-");
                          printf("Secondary token : %s\n",tok2);
                          token2 = find_kind(tok2, ONETEP_2ND_DIRECTIVES);
                          printf("Checks out as %d\n", token2);
  
                          switch (token2)
                           {
                              case LATTICE : printf("This block is for lattice parameters\n");
                              fprintf(fp, "%BLOCK LATTICE_CART\n");

                              p_this_latt_vec = p_latt_vec;
                              for (jloop = 0; jloop < 3; jloop++)
                                {
                                 for (this_atom = 0; this_atom < 3; this_atom++)
                                   {
                                      printf("looping in kloop\n");
                                      fprintf(fp, "   %14.10f ", *p_this_latt_vec);
                                      p_this_latt_vec++;
                                   } 
                                   fprintf(fp,"\n");
                                }
                              printf("Out again..this_atom %d..\n", this_atom);
                              break;
  
                              case ATMCOORDS : printf("This block is for atomic co-oridates\n");
                              fprintf(fp, "%BLOCK POSITION_ABS\n");
                              fprintf(fp, "ang\n");
                              for (jloop = 0; jloop < num_atoms; jloop++)
                                {
                                   fprintf(fp,"%3s     %14.10f  %14.10f  %14.10f\n", 
                                              p_molecule->label,  p_molecule->x, p_molecule->y, p_molecule->z);
                                   p_molecule++;
                                }
                              break;
                           }
 
                          fprintf(fp, "%%ENDBLOCK %s\n",tok2);
                          break;
          default : printf("token not recognised printing >>%s\n",file_line_ptrs[iloop]);
            fprintf(fp, "%s",file_line_ptrs[iloop]);
       }
     }
  }
else
  {
    /******** Need to include in later version ********/
    printf("Data is from another source, making writing to onetep a bit harder...\n");
    printf("Going to just write in atomistic and unit cell data into .dat template\n");
    fprintf(fp, "!\n");
    fprintf(fp, "! inter_Vasp generated ONETEP file\n");
    fprintf(fp, "! Please input keywords below\n!\n\n");
    fprintf(fp, "!\n");
    fprintf(fp, "! Species to be used in this calculation, line format:\n");
    fprintf(fp, "! (i) my symbol for the atomic species (ii) the standard element symbol\n");
    fprintf(fp, "! (iii) the atomic number Z (iv) the number of NGWFs per atom (v) the NGWF radius\n");
    fprintf(fp, "!\n");
    fprintf(fp, "%BLOCK SPECIES\n");
    fprintf(fp, "(i) (ii) (iii) (iv) (v)\n");
    fprintf(fp, "%%ENDBLOCK SPECIES\n");
    fprintf(fp, "!\n");
    fprintf(fp, "! Potentails used for each atom type (i) the atom label (ii) the potential file\n");
    fprintf(fp, "!\n");
    fprintf(fp, "%BLOCK SPECIES_POT\n");
    fprintf(fp, "(i) (ii)\n");
    fprintf(fp, "%%ENDBLOCK SPECIES_POT\n\n");
    fprintf(fp, "%BLOCK LATTICE_CART\n");
    fprintf(fp, "  ang\n");
    for (iloop = 0; iloop < 3; iloop++)
      {
        for (jloop = 0; jloop < 3; jloop++)
          {
            fprintf(fp, "   %14.10f ", *p_latt_vec);
            p_latt_vec++;
          }
        fprintf(fp,"\n");
      }
    fprintf(fp, "%%ENDBLOCK LATTICE_CART\n\n");
    fprintf(fp, "%BLOCK POSITIONS_ABS\n");
    fprintf(fp, "  ang\n");
    for (jloop = 0; jloop < num_atoms; jloop++)
      {
        fprintf(fp,"%3s     %14.10f  %14.10f  %14.10f\n", 
                p_molecule->label,  p_molecule->x, p_molecule->y, p_molecule->z);
        p_molecule++;
      }
    fprintf(fp, "%%ENDBLOCK POSITIONS_ABS\n");
    printf("Done printing into template\n");
  }

return;
}
