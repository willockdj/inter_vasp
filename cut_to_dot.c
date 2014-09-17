#include <stdio.h>

  void cut_to_dot( char *p_word )
    {
   char *p_letter;

         p_letter=p_word;
         while ( *p_letter != '\0' ) 
            {
               p_letter++;
               if ( *p_letter == '.') *p_letter = '\0'; 
            }
  return;
    }

