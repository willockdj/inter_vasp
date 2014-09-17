#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include "maxima.h"
#include "data.h"

/* routine to assign atomic mass according to element type label */

double atomic_mass_list( char *element )
{

int iloop;

for (iloop = 0; iloop < NUM_ELEMENTS; iloop++)
   {
      if (( *element ==  period_table[iloop].elem[0]
        && *(element+1) == period_table[iloop].elem[1] )
        || ( *element ==  period_table[iloop].elem[0]
        && *(element+1) == ' ' && period_table[iloop].elem[1] == '\0' ))
       
         {
         return period_table[iloop].mass ;
         }
   }

return -1.0;

}









