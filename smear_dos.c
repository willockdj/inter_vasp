/***************************************************************************/
/*** Routine to add smearing to the contents of a dos array ****************/
/*** Dave Willock February 2006 ********************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "own_maths.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */

/*---------------------------------------------------------------------------*/

void smear_dos( dos *p_dos, int num_dos, double smear )
{
  int iloop, jloop, skip, iatom;
  int noskip=FALSE;
  int num_to_skip;

  double pre_exp, centre;
  double arg;
  double two_smear2, disp, disp2;
  double updos_max, uptotdos_max, upscale, downdos_max, downtotdos_max, downscale;
  double *p_temp, *p_temp2;
  double *p_this_temp, *p_this_temp2;

  double de;

  dos *p_this_dos, *p_now_dos;

  printf("Entered smear_dos will apply smearing width %10.6f to %d dos data points\n",
                               smear, num_dos);

// Malloc the temporary arrays

  p_temp =(double*)malloc(num_dos*sizeof(double));
  p_temp2=(double*)malloc(num_dos*sizeof(double));

  pre_exp = 1/ ( smear * sqrt( 2.0 * pi));
  two_smear2 = 2.0 * smear * smear;
  updos_max = 0.0;
  downdos_max = 0.0;

  printf("pre_exponential = %10.6f\n", pre_exp);

  de = (p_dos+1)->energy - p_dos->energy;

  printf("Delta E between raw data points = %10.6f\n", de);

  p_this_dos= p_dos; p_this_temp= p_temp; p_this_temp2= p_temp2;

  for(iloop=0; iloop < num_dos; iloop++)
    {
       *p_this_temp=0.0;
       *p_this_temp2=0.0;
       centre = p_this_dos->energy;

       p_now_dos= p_dos;
       for (jloop=0; jloop < num_dos; jloop++)
         {
           disp = p_now_dos->energy - centre;
           disp2 = disp * disp;
           arg = disp2 / two_smear2;
       
/***************************************************/
/** only consider cases for which exp of gaussian **/
/** will be 1E-6 or greater                       **/
/***************************************************/
           if ( arg < 14 )
             {
               *p_this_temp += p_now_dos->up_dos * pre_exp * exp( -arg);
               *p_this_temp2 += p_now_dos->down_dos * pre_exp * exp( -arg);
             }
           p_now_dos++;
         }

/***************************************************/
/** Record maximum in smeared spectrum to allow ****/
/** scaling of the tot_dos curve                ****/
/***************************************************/
       if (*p_this_temp  > updos_max   ) updos_max   = *p_this_temp;
       if (*p_this_temp2 > downdos_max ) downdos_max = *p_this_temp2;

       p_this_dos++; p_this_temp++; p_this_temp2++;
    }

  printf("temporary smeared version created\n");
/****************************************************/
/** Copy back the smeared DOS ***********************/
/** and include smearing in integration of dos ******/
/****************************************************/
  p_this_temp= p_temp; p_this_temp2= p_temp2;
 
  p_this_dos= p_dos;
  p_this_dos->up_dos = *p_this_temp;
  p_this_dos->down_dos = *p_this_temp2;
  p_this_dos->up_totdos = *p_this_temp * de;
  p_this_dos->down_totdos = *p_this_temp2 * de;

  p_this_dos++; p_this_temp++; p_this_temp2++;

/****************************************************/
/** tot_dos is an integration so must have maximum **/
/** value at the end of the array. Use this to     **/
/** work out a scaling factor which allows dos and **/
/** tdos to be shown on same plot nicely           **/
/****************************************************/
 
  for(iloop=1; iloop < num_dos; iloop++)
    {
       p_this_dos->up_dos = *p_this_temp;
       p_this_dos->down_dos = *p_this_temp2;
       p_this_dos->up_totdos = ( (p_this_dos-1)->up_totdos + de * *p_this_temp);
       p_this_dos->down_totdos = ( (p_this_dos-1)->down_totdos + de * *p_this_temp);

       p_this_dos++; p_this_temp++; p_this_temp2++;
    }
       
  uptotdos_max = (p_dos + num_dos -1 )->up_totdos;
  upscale = updos_max / uptotdos_max;

  downtotdos_max = (p_dos + num_dos -1 )->down_totdos;
  downscale = downdos_max / downtotdos_max;

  printf("updos_max = %10.6f max up_totdos = %10.6f upscale = %10.6f downdos_max = %10.6f max down_totdos = %10.6f downscale = %10.6f\n",
                    updos_max, uptotdos_max, upscale, downdos_max, downtotdos_max, downscale);

  p_this_dos= p_dos;
  for(iloop=0; iloop < num_dos; iloop++)
    {
      p_this_dos->up_totdos = p_this_dos->up_totdos * upscale;
      p_this_dos->down_totdos = p_this_dos->down_totdos * downscale;

      printf("Calculated : %10.6f %10.6f %10.6f %10.6f\n", p_this_dos->up_dos,
                                               p_this_dos->down_dos,
                                               p_this_dos->up_totdos,
                                               p_this_dos->down_totdos);
      p_this_dos++;
    }

  free(p_temp); free(p_temp2);

  return;
}
