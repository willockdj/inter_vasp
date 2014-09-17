/* routine to convert an integer array to a string array */
#include <stdio.h>

 void int_to_string(int *p_ichar1, char *p_ichar2, int max_position )
{
 int iloop;

       for (iloop=0; iloop < max_position; iloop++)
          {
             *p_ichar2= *(p_ichar1+iloop);
             p_ichar2++;
           }
       *p_ichar2= '\0';

  return;
}
