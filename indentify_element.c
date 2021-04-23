#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <limits.h>
#include "maxima.h"
#include "data.h"

/* routine to pick out the element type from the list */

int identify_element( char *element )
{

int iloop;

for (iloop = 0; iloop < NUM_ELEMENTS; iloop++)
   {
      if (( *element ==  period_table[iloop].elem[0]
        && *(element+1) == period_table[iloop].elem[1] )
        || ( *element ==  period_table[iloop].elem[0]
        && *(element+1) == ' ' && period_table[iloop].elem[1] == '\0' ))
       
         {
         return iloop ;
         }
   }

return -1.0;

}









