/*** Routine to print out an integer array in defined blocks ***/
#include <stdio.h>

void block_print_int(int *p_array, int num, int num_per_row)
  {
    int iii, iret;

    iret=0;
      for ( iii=0; iii<num; iii++)
       {
          printf("%6d ", *p_array);
          iret++; p_array++;
          if ( iret == num_per_row ) { printf("\n"); iret=0; }
       }
    if (iret != 0 ) printf("\n");

    return;
  }

