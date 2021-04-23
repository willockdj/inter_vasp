/***************************************************************************/
/*** Read in a VASP DOSCAR and count the number of columns of data *********/
/*** Dave Willock February 2006 ********************************************/
/*** Updated April 06 to read PDOS for requested atoms in part_dos_list ****/
/*** Modified Feb 07 to count in the total number of dos and pdos columns **/
/*** for use by read_doscar ************************************************/
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

int count_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                 int num_atoms_pdos, dos *p_pdos, int *p_spd, int *p_num_dos_points,
                 int *p_num_ndos_columns, int just_count ) 

{
  int iloop, jloop, iorb;
  int found, sep;
  int skip=TRUE;
  int noskip=FALSE;
  int num_to_skip;
  int error, good_read;
  int this_num_dos;

  dos *p_this_pdos;

  double e_min, e_max;
  double this_e_min, this_e_max;
  double pop, dummy;

  char *p_key, *p_key2, *p_key3;
  char *tok, *tok2, *p_letter;

  printf("Entered count_dos\n");
//  p_key= "none";
//  p_key2= "system";
//  Assume the sixth line of the DOSCAR is the info on energy range and number of points
//
  *p_num_dos_points = -1;

  for (iloop=0; iloop < 6; iloop++)
                   tok= tok_get( fp, skip, FALSE);

//  sep = 0;


//  find_line( fp, p_key, p_key2, sep, &found, -1 );

//  if (found)
//    {
//     if (strcmp(p_key,"none") != 0 )
//        printf("Found >>%s<< in DOSCAR file\n", p_key);
//     else
//        printf("Found >>%s<< in DOSCAR file\n", p_key2);
//
       e_max = atof(tok);
       printf("got e_max as %10.6f\n", e_max);

       tok= tok_get( fp, noskip, FALSE);
       e_min = atof(tok);
       printf("got e_min as %10.6f\n", e_min);
   
       tok= tok_get( fp, noskip, FALSE);
       *p_num_dos_points = atoi(tok);
       printf("got num_dos_points as %d in count_doscar\n", *p_num_dos_points);

       printf("File e_min = %10.6f, e_max = %10.6f, num_dos_points = %d\n",
                        e_min, e_max, *p_num_dos_points);

       *p_num_ndos_columns = 0;

       if (*p_num_dos_points > 0)
         {
            if (!just_count)
               p_ndos->energy = atof(tok_get( fp, skip, FALSE));
            else
               dummy = atof(tok_get( fp, skip, FALSE));
               

/*** Count number of columns ***/

            while ( tok = tok_get( fp, noskip, FALSE) ) (*p_num_ndos_columns)++;

            printf("Number of dos columns = %d\n", *p_num_ndos_columns);

            if (!just_count)
              {
                for(iloop=0; iloop < (*p_num_dos_points)-1; iloop++)
                  {
                    p_ndos->energy = atof(tok_get( fp, skip, FALSE));
                  }   
                printf("Read from doscar\n"); 
              }
            else
              {
                for(iloop=0; iloop < (*p_num_dos_points)-1; iloop++)
                  {
                    dummy = atof(tok_get( fp, skip, FALSE));
                  }   
                printf("Skipped ndos lines for just_count in doscar\n"); 
              }
         }
//       else
//         {
//           printf("ERROR: Bad num_dos_points read in DOSCAR file in count_doscar routine.\n");
//           exit(0);
//         }
//    }
//  else
//    {
//       printf("Keyword >>%s<< not found when reading DOSCAR file\n", p_key);
//       exit(0);
//    }

  if ( need_pdos )
    {
      printf("Building pdos\n");

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

          printf("DEBUG: min_e %10.6f num_dos %d\n", this_e_min, this_num_dos);

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
          if ( *p_num_dos_points - this_num_dos != 0 )
             {
                printf("ERROR: Number of entries reported in PDOS title is inconsistent with DOSCAR header\n");
                printf("ERROR: Read %d DOSCAR header %d\n", this_num_dos, *p_num_dos_points);
                exit(0);
             }

                       
          tok=tok_get( fp, skip, FALSE);
                
/*** Count number of pdos columns and add to total ***/
 
          printf("Counting pdos orbitals\n");

          iorb=0;
          while ( tok = tok_get( fp, noskip, FALSE) )
            {
              printf("tok: %s\n" ,tok);
              iorb++;
             }
                     
          printf("Total number of columns = %d\n", iorb);
         
     }
  return iorb;
}
