/* routine to copy one string array segment to another */
#include <stdio.h>
#include <limits.h>

 void copy_int( int *p_ichar1, int *p_ichar2, 
                                  int min_position, int max_position )
{
 int iloop;

   if (min_position < max_position)
       for (iloop=min_position; iloop <= max_position; iloop++)
          {
             *p_ichar2= *(p_ichar1+iloop);
             p_ichar2++;
           }
   else 
       printf(" Warning: attempt to copy string with range minimum greater than range maximum. ");
}
