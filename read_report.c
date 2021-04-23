/***************************************************************************/
/*** Read in a VASP REPORT file from and MD run ****************************/
/*** Dave Willock Nov. 2017           **************************************/
/***                                  **************************************/
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

void cut_to_uscore( char *p_word );

void find_line(FILE *fp, char *p_key, char *p_key2, int sep, int *p_found, int max_lines);

void find_line_with_stopper(FILE *fp, char *p_key, char *p_key2, char *p_stopper, int sep, int *p_found, int max_lines);

void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                            double *p_abc);

int read_force_block( FILE *fp, e_vec *p_forces, int num_atoms);
/*---------------------------------------------------------------------------*/

/* read in a REPORT file for VASP MD analysis */

/** fp is the pointer to the file that the calling programme has set ***/

int read_report( FILE *fp, int *p_num_frames, int just_count,
                 double *p_etot, double *p_epot, double *p_ekin, double *p_tsim, double *p_tinst)

{
/***********************************************************************/
/** All variables must be declared in c-codes. *************************/
/***********************************************************************/
  int iloop, jloop, skip, iatom, found, sep, Tsep;

  char *p_key, *p_key2, *p_key3, *p_stopper;
  char *p_Tkey, *p_Tkey2;
  char *tok, cdum[10], cdum2[15];
  char *p_label;

  double epot, ekin, econst, eps, es, tsim, tinst; 

  atom *p_atom;

/*** Capture the md algorithm flag in this variable ***/
  int mdalgo;


/* DEBUG */
  printf("Arrived in read_report...Will see end of file from >>%s<<\n", END_OF_INPUT);
/* get all atom labels from POTCAR lines */
/* up to stopper                         */

/*** These keys are used to position the file pointer in the file to allow data to be found ***/

/*** This code is taken from another reading routine and requires re-structuring ***/

      p_key= "MDALGO";
      p_key2= "none";
      p_stopper= "original";
      sep = 0;

      find_line_with_stopper( fp, p_key, p_key2, p_stopper, sep, &found, -1 );

      if (found)
        {
          printf("Found the line\n");
          tok= tok_get( fp, FALSE, FALSE); 
          printf("Debug --- tok >>%s<< \n", tok);
          tok= tok_get( fp, FALSE, FALSE); 

          mdalgo=atoi(tok);
          printf("Debug --- tok >>%s<< that is mdalgo=%d \n", tok, mdalgo);

          
/*** convert string to value ****/
 
        }
      else
        {
          printf("NOT found the line\n");
        } 

      found=TRUE;
      p_key= "MD";
      p_key2= "No.";
      sep = 1;


    if (just_count)
      {
        *p_num_frames=0;       

        	while (found) {
            find_line( fp, p_key, p_key2, sep, &found, -1 );
            if (found) ++*p_num_frames;       
          	}

        printf("Found %d frames\n", *p_num_frames);
        return 0;
     }
/*** Now read the data for real ***/
     
   found= TRUE;
   p_key="E_tot";
   p_key2= "E_const";
   sep = 2;	
/*** Keys for temperature reading ***/
   p_Tkey="T_sim";
   p_Tkey2="T_inst";
   Tsep=0;

/*!feof(fp) **/	
   while (found) {
		find_line(fp, p_key,p_key2,sep, &found, -1);
		if (found)
			{
//
//		 Note that fscanf is expecting pointers to the variables it reads
//
              fscanf(fp, "%s %le %le %le %le %le %le", &cdum[0], p_etot, p_epot, p_ekin, &econst, &eps, &es); 
              printf( "%s %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",cdum, *p_etot,*p_epot,*p_ekin, econst, eps, es); 

              p_etot++;
	      p_epot++;
              p_ekin++;
/*** Look for temperature data ***/
              find_line(fp, p_Tkey, p_Tkey2, Tsep, &found, -1);
              if (found)
                {
                 fscanf(fp, "%s %le %le",&cdum[0], p_tsim, p_tinst);
                 printf( "%s %10.6f %10.6f\n",cdum, *p_tsim, *p_tinst);
		p_tsim++;
		p_tinst++;
                }
            }	
       	 }
   return 0;
}

