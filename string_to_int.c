/* routine to convert an string array to a int array */
#include <stdio.h>

 void string_to_int(char *p_ichar2, int *p_ichar1, int max_position )
{
 int iloop;

       while ( *p_ichar2 != '\0' && iloop < max_position )
          {
             *p_ichar1= *p_ichar2;
             p_ichar1++;
             p_ichar2++;
             iloop++;
           }

  return;
}
