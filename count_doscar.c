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
                 int num_atoms_pdos, dos *p_pdos, int *p_spd ) 

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
  double pop;

  char *p_key, *p_key2, *p_key3;
  char *tok, *tok2, *p_letter;

  printf("Entered count_dos\n");
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

       printf("File e_min = %10.6f, e_max = %10.6f, num_dos = %d\n",
                        e_min, e_max, num_dos);

       iorb = 0;

       if (num_dos > 0)
         {
            p_ndos->energy = atof(tok_get( fp, skip, FALSE));

/*** Count number of columns ***/

            while ( tok = tok_get( fp, noskip, FALSE) )
                   {
                      iorb++;
                   }
            printf("Number of dos columns = %d\n", iorb);

            for(iloop=0; iloop < num_dos-1; iloop++)
              {
                 p_ndos->energy = atof(tok_get( fp, skip, FALSE));
              }   
            printf("Read from doscar\n"); 
  
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

          printf("Reading this pdos contribution\n");
          p_this_pdos = p_pdos;
         
                       
          tok=tok_get( fp, skip, FALSE);
                
/**** First acceptable contribution to pdos ****/
          p_this_pdos->energy = atof(tok);

          p_this_pdos->up_dos =0.0;
         
/*** Count number of pdos columns and add to total ***/
 
          printf("Counting pdos orbitals\n");

          while ( tok = tok_get( fp, noskip, FALSE) )
            {
              iorb++;
             }
                     
          printf("Total number of columns = %d\n", iorb);
            
         
     }
  return iorb;
}
