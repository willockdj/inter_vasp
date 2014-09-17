/***************************************************************************/
/*** Read in a VASP DOSCAR *************************************************/
/*** Dave Willock February 2006 ********************************************/
/*** Updated April 06 to read PDOS for requested atoms in part_dos_list ****/
/*** Modified Feb 07 to take the total number of columns from **************/
/*** count_doscar.c and read in the DOSCAR file using this information *****/
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

/*---------------------------------------------------------------------------*/

/* read in a doscar file from VASP */

int read_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                 int num_atoms_pdos, dos *p_pdos, int *p_spd, int num_columns )
{
  int iloop, jloop, iatom, iorb;
  int found, sep;
  int skip=TRUE;
  int noskip=FALSE;
  int num_to_skip;
  int error, good_read;
  int num_dos, this_num_dos;
  int pdos_atom;

  dos *p_this_pdos;

  double e_min, e_max;
  double this_e_min, this_e_max;
  double pop, pop2;

  char *p_key, *p_key2, *p_key3;
  char *tok, *tok2, *p_letter;

  printf("Entered read_dos\n");
  p_key= "none";
  p_key2= "system";
  sep = 0;
  num_dos = -1;

  tok= tok_get( fp, skip, FALSE);

  find_line( fp, p_key, p_key2, sep, &found, -1 );

  if (found)
    {
       if (strcmp(p_key,"none") != 0 )
          printf("Found >>%s<< in DOSCAR file\n", p_key);
       else
          printf("Found >>%s<< in DOSCAR file\n", p_key2);

       tok= tok_get( fp, skip, FALSE);
       e_max = atof(tok);
       printf("got e_max as %10.6f\n", e_max);

       tok= tok_get( fp, noskip, FALSE);
       e_min = atof(tok);
       printf("got e_min as %10.6f\n", e_min);
   
       tok= tok_get( fp, noskip, FALSE);
       num_dos = atoi(tok);
       printf("got num_dos as %d\n", num_dos);

       if (num_dos > MAX_DOS)
          {
             printf("ERROR : Number of points in DOSCAR file (%d) exceeds MAX_DOS (%d).\n",
                          num_dos, MAX_DOS);
             exit(0);
          }

       printf("File e_min = %10.6f, e_max = %10.6f, num_dos = %d\n",
                        e_min, e_max, num_dos);
       printf("Number of columns = %d\n", num_columns);

       if (num_dos > 0)
         {
            for(iloop=0; iloop < num_dos; iloop++)
              {
/*** If spin restricted ... ***/
                 if (num_columns == 2 || num_columns == 5 || num_columns == 11)
                   {
                     p_ndos->energy = atof(tok_get( fp, skip, FALSE));

                     p_ndos->up_dos = atof(tok_get( fp, noskip, FALSE));
                     p_ndos->up_totdos = atof(tok_get( fp, noskip, FALSE));

                     printf("Read from doscar: %10.6f  %10.6f  %10.6f\n",
                                                  p_ndos->energy,
                                                  p_ndos->up_dos,
                                                  p_ndos->up_totdos); 
                   }
/*** If spin unrestricted ... ***/
                 else if (num_columns == 4 || num_columns == 10 || num_columns == 22)
                   {
                     p_ndos->energy = atof(tok_get( fp, skip, FALSE));

                     p_ndos->up_dos = atof(tok_get( fp, noskip, FALSE));
                     p_ndos->down_dos = atof(tok_get( fp, noskip, FALSE));
                     p_ndos->up_totdos = atof(tok_get( fp, noskip, FALSE));
                     p_ndos->down_totdos = atof(tok_get( fp, noskip, FALSE));

                     printf("Read from doscar: %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n",
                                                  p_ndos->energy,
                                                  p_ndos->up_dos,
                                                  p_ndos->down_dos,
                                                  p_ndos->up_totdos,
                                                  p_ndos->down_totdos);
                   }
                        
                 p_ndos++;
              }
         }
       else
         {
           printf("ERROR: Bad num_dos read in DOSCAR file\n");
           exit(0);
         }
    }
  else
    {
       printf("Keyword >>%s<< not found when reading DOSCAR file\n", p_key);
       exit(0);
    }

  if ( need_pdos )
    {
      printf("Building pdos\n");

      iatom=1;
      pdos_atom=0;
      while (iatom <= *(p_part_dos_list+num_atoms_pdos))
        {
/**** Next line should be pdos header line ****/

          tok= tok_get( fp, skip, FALSE);
/**** Alert the user if this fails to check that PDOS was created ****/
          if (!tok || strcmp(tok, "Thatsit") == 0)
           {
             printf("ERROR: Failure while reading PDOS data, did you ask for PDOS (LORBIT = 12) in INCAR?\n");
             exit(0);
           }
          this_e_max = atof(tok);

          tok= tok_get( fp, noskip, FALSE);
          this_e_min = atof(tok);
   
          tok= tok_get( fp, noskip, FALSE);
          this_num_dos = atoi(tok);

/**** Test consistency of this pdos header line with the original information from the dos ****/
          if (fabs(e_min-this_e_min) > 0.0001)
             {
                printf("ERROR: Minimum energy given in PDOS title is inconsistent with DOSCAR header\n");
                printf("ERROR: Read %10.6f DOSCAR header %10.6f\n", this_e_min, e_min);
                exit(0);
             }
          if (fabs(e_max-this_e_max) > 0.0001)
             {
                printf("ERROR: Maximum energy given in PDOS title is inconsistent with DOSCAR header\n");
                printf("ERROR: Read %10.6f DOSCAR header %10.6f\n", this_e_max, e_max);
                exit(0);
             }
          if ( num_dos - this_num_dos != 0 )
             {
                printf("ERROR: Number of entries reported in PDOS title is inconsistent with DOSCAR header\n");
                printf("ERROR: Read %d DOSCAR header %d\n", this_num_dos, num_dos);
                exit(0);
             }

/**** See if we need this atoms data ****/
          printf("iatom = %d\n", iatom);
          if (iatom == *(p_part_dos_list+pdos_atom))
            {
              printf("Reading this pdos contribution\n");
              p_this_pdos = p_pdos;
                       
                   printf("Entered first atom in PDOS list\n");
                   for(iloop=0; iloop < num_dos; iloop++) 
                     {

                       tok=tok_get( fp, skip, FALSE); 

                       if (pdos_atom == 0)
                         {
/**** First acceptable contribution to pdos ****/
                          p_this_pdos->energy = atof(tok);
/**** Second or higher acceptable contribution to pdos ****/
                         }
                       else
                         {
                            if (fabs(atof(tok)-p_this_pdos->energy) > 0.0001)
                             {
                                printf("ERROR: Energy scale changed during PDOS read\n");
                                printf("ERROR: Read %10.6f reference is %10.6f\n", atof(tok),
                                            p_this_pdos->energy);
                                exit(0);
                             }
                         }

                       if (pdos_atom == 0) p_this_pdos->up_dos =0.0;
                       iorb = 0;

                       while ( (tok = tok_get( fp, noskip, FALSE)) )
                         {
/**** s, p, d orbitals ****/
                           if (num_columns == 5)
                             {
                               if (*(p_spd+iorb)) p_this_pdos->up_dos += atof(tok);
                               iorb++;
                             } 
/**** s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2 orbitals ****/
                           else if (num_columns == 11)
                             {

                               if (iorb == 0)
                                 {
                                   pop = atof(tok);
                                 }
                               else if (iorb == 1)
                                 {
                                   pop = atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok); 
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                 }
                               else if (iorb == 2)
                                 {
                                   pop = atof(tok); 
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok); 
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok); 
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok); 
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                 }
                               else if (iorb > 2)
                                 {
                                    printf("ERROR: Too many orbitals in PDOS entries of DOSCAR file.\n");
                                    printf("ERROR: Only expecting spd level but iorb = %d.\n", iorb);
                                    exit(0);
                                 }
                               if (*(p_spd+iorb)) p_this_pdos->up_dos += pop;
                               iorb++;   
                             }
/**** s, p, d orbitals, spin unrestricted ****/
                           else if (num_columns == 10) 
                             {
                               if (*(p_spd+iorb)) 
                                 {
                                   p_this_pdos->up_dos += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   p_this_pdos->down_dos += atof(tok);
                                 }
                               iorb++;
                             }
/**** s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2 orbitals, spin unrestricted ****/
                           else if (num_columns == 22)
                             {
                               if (iorb == 0)
                                 {
                                   pop = atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 = atof(tok);
                                 }
                               else if (iorb == 1)
                                 {
                                   pop = atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 = atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 += atof(tok);
                                 }
                               else if (iorb == 2)
                                 {
                                   pop = atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 = atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop += atof(tok);
                                   tok = tok_get( fp, noskip, FALSE);
                                   pop2 += atof(tok);
                                 }
                               else if (iorb > 2)
                                 {
                                    printf("ERROR: Too many orbitals in PDOS entries of DOSCAR file.\n");
                                    printf("ERROR: Only expecting spd level but iorb = %d.\n", iorb);
                                    exit(0);
                                 }

                               if (*(p_spd+iorb)) 
                                 {
                                   p_this_pdos->up_dos += pop;
                                   p_this_pdos->down_dos += pop2;
                                 }
                               iorb++;
                             }
                           else
                             {
                               printf("ERROR: Unrecognised number of columns in pdos section of DOSCAR file\n");
                               exit(0);
                             }    
                           }
                           p_this_pdos++;
                             
                    }
              pdos_atom++;
            }
      
          else
            {
/**** If not needed read in the lines to ensure proper skipping *****/

               for(iloop=0; iloop < num_dos; iloop++) tok=tok_get( fp, skip, FALSE);
            }

          iatom++;
         }
    }

  return num_dos;
}
