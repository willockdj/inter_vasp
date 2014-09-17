#include <stdio.h>
#include "global_values.h"

int read_star_format( char *p_word, char *p_after )
 {
   char *p_letter, *p_after_letter;
   int have_star;

   have_star=FALSE;

   p_letter=p_word;
   p_after_letter=p_after;

   while ( *p_letter != '\0' ) 
     {
        if (have_star)
          {
            *p_after_letter= *p_letter;
            p_after_letter++;
          }
        else if ( *p_letter == '*') 
          {
            have_star=TRUE;
            *p_letter = '\0'; 
          }
        p_letter++;
     }

  if (have_star) *p_after_letter = '\0'; 
  return have_star;
 }

