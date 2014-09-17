/* routine to move an indexing variable to the next space in a string
 * or the end of string by \n */
#include <stdio.h>
#include <limits.h>

 int next_space( int *p_ichar, int start, int num_of_chars )
{
 int iloop;
 
 for (iloop= start; iloop <= num_of_chars+1; iloop++)
   {
      if ( *(p_ichar+iloop) == ' ' || *(p_ichar+iloop) == '\n') return iloop;
   }
 return -1;
}
