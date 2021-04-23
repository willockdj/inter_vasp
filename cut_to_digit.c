#include <stdio.h>
#include <ctype.h>

  void cut_to_digit( char *p_word )
    {
   char *p_letter;

         p_letter=p_word;
         while ( *p_letter != '\0' ) 
            {
               p_letter++;
               if ( isdigit(*p_letter) ) *p_letter = '\0'; 
            }
  return;
    }

