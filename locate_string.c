#include <stdio.h>
#include <limits.h>

/* subroutine for searhing the read in arrays for a string of defined length */

int locate_string( char *p_key, int *p_char, int num_of_chars )

{

int iloop;
char *p_start;

  p_start= p_key;

  for (iloop=0; iloop <= num_of_chars; ++iloop)
    {
      if ( *(p_char+iloop) == *p_key ) p_key++; 
                                     else p_key = p_start;
      if (*(p_key) == '\0') return 1; 
    }

  return 0;
}




