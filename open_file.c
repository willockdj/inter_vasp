/******************************************************************************/
/* open_file.c opens file from file pointer                                   */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

void open_file(FILE **p_file, char *p_filename, char *p_status)
{
if ((*p_file = fopen(p_filename,p_status)) == NULL)
          {
            fprintf(stderr, "Unable to open file >>%s<<.\tAborting\n",
                        p_filename);
            exit(EXIT_FAILURE);
          }
return;
}

