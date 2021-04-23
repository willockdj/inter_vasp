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

void read_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                  int num_atoms_pdos, dos *p_pdos, int *p_spd, int num_ndos_columns,
                  int num_pdos_columns, int num_dos )
{
  int i, iloop, jloop, iatom, iorb;
  int found, sep;
  int skip=TRUE;
  int noskip=FALSE;
  int started_pop;
  int num_to_skip;
  int error, good_read;
  int this_num_dos;
  int pdos_atom;

  dos *p_this_pdos;

  double e_min, e_max;
  double this_e_min, this_e_max;
  double pop, pop2;

  char *p_key, *p_key2, *p_key3;
  char *tok, *tok2, *p_letter;

  printf("Entered read_dos\n");

// Position for the info line

  for (iloop=0; iloop < 6; iloop++)
                   tok= tok_get( fp, skip, FALSE);

       e_max = atof(tok);
       printf("got e_max as %10.6f\n", e_max);

       tok= tok_get( fp, noskip, FALSE);
       e_min = atof(tok);
       printf("got e_min as %10.6f\n", e_min);
   
       tok= tok_get( fp, noskip, FALSE);
       this_num_dos = atoi(tok);
       printf("got this_num_dos as %d\n", this_num_dos);

       if ( this_num_dos != num_dos ) 
         {
            printf("ERROR : Number of points in DOSCAR file (%d) according to read_doscar\n", this_num_dos);
            printf("ERROR : does not agree with that from count_doscar (%d).\n", num_dos);
             exit(0);
         }

       printf("File e_min = %10.6f, e_max = %10.6f, num_dos = %d\n",
                        e_min, e_max, num_dos);
       printf("Number of ndos_columns = %d\n", num_ndos_columns);
       printf("Number of pdos_columns = %d\n", num_pdos_columns);

       if (num_dos > 0)
         {
            for(iloop=0; iloop < num_dos; iloop++)
              {
/*** If spin restricted ... note that num_ndos_columns does not include the energy column ***/
                 if (num_ndos_columns == 2 ) 
                   {
                     p_ndos->energy = atof(tok_get( fp, skip, FALSE));

                     p_ndos->up_dos = atof(tok_get( fp, noskip, FALSE));
                     p_ndos->up_totdos = atof(tok_get( fp, noskip, FALSE));

                     printf("Read from doscar: %10.6f  %10.6f  %10.6f\n",
                                                  p_ndos->energy,
                                                  p_ndos->up_dos,
                                                  p_ndos->up_totdos); 
                   }
/*** If spin unrestricted ... note that num_ndos_columns does not include the energy column ***/
                 else if (num_ndos_columns == 4 )
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
                 else
                   {
                      printf("ERROR: Unrecognised number of ndos columns in DOSCAR file: num_ndos_columns = %d\n", num_ndos_columns);
                      printf("ERROR: Check pdos section of the DOSCAR file.\n");
                      exit(0);
                   }
                 p_ndos++;
              }
          }

  if ( need_pdos )
    {
      printf("Building pdos for %d atoms\n", num_atoms_pdos);

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

          started_pop=FALSE;
          if (iatom == *(p_part_dos_list+pdos_atom))
            {
              printf("Reading this pdos contribution for atom %d, index %d num_dos= %d\n", 
                                                                          iatom, pdos_atom, num_dos);
              p_this_pdos = p_pdos;
                       
              for (iloop=0; iloop < num_dos; iloop++) 
                 {
                   tok=tok_get( fp, skip, FALSE); 

                   if (pdos_atom == 0)
                     {
/**** First acceptable contribution to pdos, record the energy column ****/
                       p_this_pdos->energy = atof(tok);
                     }
/**** Second or higher acceptable contribution to pdos test there are no changes ****/
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

                   if (iloop == 0) 
                     {
                        printf("Gathering data for pdos,from %d columns.\n", num_pdos_columns);
                        if ( *p_spd )     printf("spd[0] set T\n"); else printf("spd[0] set F\n");
                        if ( *(p_spd+1) ) printf("spd[1] set T\n"); else printf("spd[1] set F\n");
                        if ( *(p_spd+2) ) printf("spd[2] set T\n"); else printf("spd[2] set F\n");
                        if ( *(p_spd+3) ) printf("spd[3] set T\n"); else printf("spd[3] set F\n");
                     }

/**** read the population data for pdos lines ****/
                   while ( (tok = tok_get( fp, noskip, FALSE)) )
                     {
/**** s orbital only restricted ****/
                       if (num_pdos_columns == 1)
                         {
                           if (*p_spd) p_this_pdos->up_dos = atof(tok);
                         } 
/**** s, px, py, pz orbitals restricted ****/
                       else if (num_pdos_columns == 4)
                         {
                           p_this_pdos->up_dos = 0.0;
                           if (*p_spd) p_this_pdos->up_dos += atof(tok);

                           if (*(p_spd+1)) 
                             {
                               for (iorb=0; iorb < 3; iorb++)
                                  {
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->up_dos += atof(tok);
                                  }
                             }
                           else
                             {
                               for (iorb=0; iorb < 3; iorb++)
                                    tok = tok_get( fp, noskip, FALSE);
                             }
                         }
/**** s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2 orbitals restricted ****/
                       else if (num_pdos_columns == 9)
                         {
                           p_this_pdos->up_dos = 0.0;
                           if (*p_spd) p_this_pdos->up_dos += atof(tok);

                           if (*(p_spd+1)) 
                             {
                               for (iorb=0; iorb < 3; iorb++)
                                  {
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->up_dos += atof(tok);
                                  }
                             }
                           else
                             {
                               for (iorb=0; iorb < 3; iorb++)
                                           tok = tok_get( fp, noskip, FALSE);
                             }


                           if (*(p_spd+2)) 
                             {
                               for (iorb=0; iorb < 5; iorb++)
                                  {
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->up_dos += atof(tok);
                                  }
                             }
                           else
                             {
                               for (iorb=0; iorb < 5; iorb++)
                                    tok = tok_get( fp, noskip, FALSE);
                             }
                          }
/**** s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2 and f-set orbitals restricted ****/
                       else if (num_pdos_columns == 16)
                         {
                           p_this_pdos->up_dos = 0.0;
                           if (*p_spd) p_this_pdos->up_dos += atof(tok);

                           if (*(p_spd+1)) 
                             {
                               for (iorb=0; iorb < 3; iorb++)
                                  {
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->up_dos += atof(tok);
                                  }
                             }
                           else
                             {
                               for (iorb=0; iorb < 3; iorb++)
                                           tok = tok_get( fp, noskip, FALSE);
                             }


                               if (*(p_spd+2)) 
                                 {
                                   for (iorb=0; iorb < 5; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                      }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 5; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }


                               if (*(p_spd+3)) 
                                 {
                                   for (iorb=0; iorb < 7; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                      }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 7; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }
                             }
/**** s orbital only unrestricted ****/
                           else if (num_pdos_columns == 2)
                             {
                               if (*p_spd) 
                                 {
                                    p_this_pdos->up_dos = atof(tok);
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->down_dos = atof(tok);
                                 }
                             } 
/**** s, px, py, pz orbitals unrestricted ****/
                           else if (num_pdos_columns == 8)
                             {
                               p_this_pdos->up_dos   = 0.0;
                               p_this_pdos->down_dos = 0.0;

                               if (*p_spd) 
                                 {
                                    p_this_pdos->up_dos = atof(tok);
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->down_dos = atof(tok);
                                 }

                               if (*(p_spd+1)) 
                                 {
                                   for (iorb=0; iorb < 3; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->down_dos += atof(tok);
                                      }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 6; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }
                             }
/**** s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2 orbitals unrestricted ****/
                           else if (num_pdos_columns == 18)
                             {
                               if (iloop == 0) printf("reading from 18 columns\n");
                               p_this_pdos->up_dos = 0.0;
                               p_this_pdos->down_dos = 0.0;

                               if (*p_spd) 
                                 {
                                    p_this_pdos->up_dos = atof(tok);
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->down_dos = atof(tok);

                                    if ( iloop == 0 ) 
                                       {
                                 printf("Including s-orbital contributions of %10.6f (up) %10.6f (down)\n", 
                                          p_this_pdos->up_dos ,p_this_pdos->down_dos );
                                       }
                                 }
                               else
                                 {
                                    tok = tok_get( fp, noskip, FALSE);
                                 }

                               if (*(p_spd+1)) 
                                 {
                                   for (iorb=0; iorb < 3; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->down_dos += atof(tok);
                                      }
                                    if ( iloop == 0 ) 
                                       {
                                 printf("Including p-orbital contributions total now %10.6f (up) %10.6f (down)\n", 
                                          p_this_pdos->up_dos ,p_this_pdos->down_dos );
                                       }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 6; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }


                               if (*(p_spd+2)) 
                                 {
                                   for (iorb=0; iorb < 5; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->down_dos += atof(tok);
                                      }
                                    if ( iloop == 0 ) 
                                       {
                                 printf("Including d-orbital contributions total now %10.6f (up) %10.6f (down)\n", 
                                          p_this_pdos->up_dos ,p_this_pdos->down_dos );
                                       }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 10; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }
                              }
/**** s, px, py, pz, dxy, dyz, dxz, dx2-y2, dz2 and f-set orbitals unrestricted ****/
                           else if (num_pdos_columns == 32)
                             {
                               p_this_pdos->up_dos = 0.0;
                               p_this_pdos->down_dos = 0.0;

                               if (*p_spd) 
                                 {
                                    p_this_pdos->up_dos = atof(tok);
                                    tok = tok_get( fp, noskip, FALSE);
                                    p_this_pdos->down_dos = atof(tok);
                                 }
                               else
                                 {
                                    tok = tok_get( fp, noskip, FALSE);
                                 }

                               if (*(p_spd+1)) 
                                 {
                                   for (iorb=0; iorb < 3; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->down_dos += atof(tok);
                                      }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 6; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }


                               if (*(p_spd+2)) 
                                 {
                                   for (iorb=0; iorb < 5; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->down_dos += atof(tok);
                                      }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 10; iorb++)
                                               tok = tok_get( fp, noskip, FALSE);
                                 }


                               if (*(p_spd+3)) 
                                 {
                                   for (iorb=0; iorb < 7; iorb++)
                                      {
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->up_dos += atof(tok);
                                        tok = tok_get( fp, noskip, FALSE);
                                        p_this_pdos->down_dos += atof(tok);
                                      }
                                 }
                               else
                                 {
                                   for (iorb=0; iorb < 14; iorb++)
                                        tok = tok_get( fp, noskip, FALSE);
                                 }
                             }
                          else 
                             {
                                printf("ERROR: Unexpected number of columns in pdos of DOSCAR file.\n");
                                printf("ERROR: Counted: %d columns?\n", num_pdos_columns);
                                exit(0);
                             }
                           }
//                           printf("Exit for DEBUG \n");
//                           exit(0);
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

  return;
}
