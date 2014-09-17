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
  double temp[MAX_DOS], temp2[MAX_DOS], arg;
  double two_smear2, disp, disp2;
  double updos_max, uptotdos_max, upscale, downdos_max, downtotdos_max, downscale;

  double de;

  dos *p_this_dos, *p_now_dos;

  printf("Entered smear_dos will apply smearing width %10.6f\n",
                               smear);

  pre_exp = 1/ ( smear * sqrt( 2.0 * pi));
  two_smear2 = 2.0 * smear * smear;
  updos_max = 0.0;
  downdos_max = 0.0;

  printf("pre_exponential = %10.6f\n", pre_exp);

  de = (p_dos+1)->energy - p_dos->energy;

  printf("Delta E between raw data points = %10.6f\n", de);

  p_this_dos= p_dos;
  for(iloop=0; iloop < num_dos; iloop++)
    {
       temp[iloop]=0.0;
       temp2[iloop]=0.0;
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
               temp[iloop] += p_now_dos->up_dos * pre_exp * exp( -arg);
               temp2[iloop] += p_now_dos->down_dos * pre_exp * exp( -arg);
             }
           p_now_dos++;
         }

/***************************************************/
/** Record maximum in smeared spectrum to allow ****/
/** scaling of the tot_dos curve                ****/
/***************************************************/
       if (temp[iloop] > updos_max ) updos_max = temp[iloop];
       if (temp2[iloop] > downdos_max ) downdos_max = temp2[iloop];

       p_this_dos++;
    }

/****************************************************/
/** Copy back the smeared DOS ***********************/
/** and include smearing in integration of dos ******/
/****************************************************/
  p_this_dos= p_dos;
  p_this_dos->up_dos = temp[0];
  p_this_dos->down_dos = temp2[0];
  p_this_dos->up_totdos = temp[0] * de;
  p_this_dos->down_totdos = temp2[0] * de;
  p_this_dos++;

/****************************************************/
/** tot_dos is an integration so must have maximum **/
/** value at the end of the array. Use this to     **/
/** work out a scaling factor which allows dos and **/
/** tdos to be shown on same plot nicely           **/
/****************************************************/
 
  for(iloop=1; iloop < num_dos; iloop++)
    {
       p_this_dos->up_dos = temp[iloop];
       p_this_dos->down_dos = temp2[iloop];
       p_this_dos->up_totdos = ( (p_this_dos-1)->up_totdos + de * temp[iloop]);
       p_this_dos->down_totdos = ( (p_this_dos-1)->down_totdos + de * temp[iloop]);
       p_this_dos++;
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

  return;
}
