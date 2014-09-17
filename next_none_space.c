/* routine to move an indexing variable to the next none space in a string */
#include <stdio.h>
#include <limits.h>

 int next_none_space( int *p_ichar, int start, int num_of_chars )
{
 int iloop;
 
 for (iloop= start; iloop <= num_of_chars; iloop++)
                         if ( *(p_ichar+iloop) != ' ') return iloop;
    
 return -1;
}
