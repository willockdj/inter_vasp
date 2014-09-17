#include <stdio.h>
#include <limits.h>
#include "global_values.h"

/* routine to read in a line into the character array ichar */
/* now altered to read from the file pointed to by fp
 * if this is the standard input send in stdin */

int read_line(FILE *fp, int *p_ichar)

{
 int idave,num_of_chars;

/* read in one line upto 81 characters 
  note that c arrays have an element zero and this is used here */

      num_of_chars=-1;
      while ((idave= getc(fp)) != '\n' && num_of_chars<=LINESIZ )
      {
          if (idave == EOF)
            {
               printf("End of file detected in read_line\n");
               return -10;
            }
          ++num_of_chars;
          *(p_ichar+num_of_chars) = idave;   
      }

/* mark the end of the string read in with \0 */

  *(p_ichar+num_of_chars+1)= '\n';
  *(p_ichar+num_of_chars+2)= '\0';
  return num_of_chars;
}

