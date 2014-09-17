#include <stdio.h>
#include <limits.h>
#include "global_values.h"

/* subroutine for comparing two strings over a given length */

int compare_strings( char *p_ichar1, char *p_ichar2 )

{

int iloop;

   iloop=0;
   while (*(p_ichar1+iloop) != '\0' && *p_ichar2+iloop != '\0')
      {
      if ( *(p_ichar1+iloop) != *(p_ichar2+iloop) ) return FALSE;
      iloop++;
      }

  return TRUE;
}




