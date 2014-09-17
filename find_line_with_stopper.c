/*******************************************************/
/*** Generic routine to find a line containing up to ***/
/*** two key words at a specified spacing            ***/
/*** found logical added to allow multiple reading   ***/
/*** of force lines from OUTCAR file.                ***/
/*** max_lines added to allow limited search forward ***/
/*** if max_lines is positive search will continue   ***/
/*** for that number of lines, if negative search    ***/
/*** continues to end of file                        ***/
/*** Update to allow first key to be "none". This    ***/
/*** allows searching for a keyword that is not the  ***/
/*** first on the line!                              ***/
/***                                                 ***/
/*** With stopper added Feb. 09                      ***/
/*** The stopper is a line that stops the search     ***/
/*** and returns false.                              ***/
/*** Dave Willock May 04 Last update Feb. 09         ***/
/*******************************************************/
  #include <stdio.h>
  #include <string.h>
  #include "maxima.h"
  #include "global_values.h"
  #include "structures.h"
  #include "debug.h"
  #include "reader.h"

  char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

  void find_line_with_stopper(FILE *fp, char *p_key, char *p_key2, char *p_stopper, int sep, int *p_found, int max_lines)
    {
      int skip = TRUE;
      int noskip = FALSE;
      int num_read;
      int iloop, done;
      int hit_stopper;

      char *tok, *tok2;

      printf("In find_line_with_stopper looking for %s and %s %d apart, will watch for stopper: %s\n", p_key, p_key2, sep, p_stopper);
      done = FALSE;

/****************************************************/
/*** If p_key is not "none" then read a new tok for */
/*** testing. Otherwise we must be only testing on  */
/*** p_key2, so tok can be anything and read in a   */
/*** new tok2.                                      */
/****************************************************/
      if (strcmp(p_key, "none") != 0) 
        {
          tok = tok_get( fp, skip, FALSE);
        }
      else
        {
           tok2 = tok_get( fp, skip, FALSE);
           tok = "none";
        }
          
      num_read=1;

/****************************************************/
/*** If we are testing on p_key2 need to get the ****/
/*** right position (sep) on the line            ****/
/****************************************************/
      if ( strcmp(p_key2, "none") != 0) 
         {
           for (iloop=0; iloop <= sep; iloop++) tok2= tok_get( fp, noskip, FALSE);
         }
      else
         {
           tok2 = "none";
         }
/****************************************************/
/*** See if first line has hit the spot *************/
/****************************************************/
     if (strcmp(p_key, "none") != 0 && strcmp(p_key2, "none") != 0)   
       {
          printf("Testing first line >>%s<< >>%s<<\n", tok, tok2);
          if (tok && strcmp(p_key, tok) == 0 && tok2 && strcmp(p_key2, tok2) == 0) done=TRUE;

          if (tok && strcmp(p_stopper, tok) == 0)
            {
               done=TRUE;
               hit_stopper=TRUE;
            }
       }
/****************************************************/
/** Case where p_key is not "none" but p_key2 is  ***/
/****************************************************/
     else if (strcmp(p_key, "none") != 0)
       {
          if (tok && strcmp(p_key, tok) == 0) done=TRUE;

          if (tok && strcmp(p_stopper, tok) == 0)
            {
               done=TRUE;
               hit_stopper=TRUE;
            }
       }
/****************************************************/
/** Case where p_key2 is not "none" but p_key is  ***/
/****************************************************/
     else if (strcmp(p_key2, "none") != 0)
       {
          if (tok2 && strcmp(p_key2, tok2) == 0) done=TRUE;

          if (tok2 && strcmp(p_stopper, tok2) == 0)
            {
               done=TRUE;
               hit_stopper=TRUE;
            }
       }


      if (done) printf("Returning after first line\n");
/****************************************************/
/*** Loop to test remaining lines *******************/
/****************************************************/

      while (!done)
         {
            if (num_read > 1 ) 
              {
                if (strcmp(p_key, "none") != 0) 
                   {
                        tok= tok_get( fp, skip, FALSE);
                   }
                else
                   {
                        tok2= tok_get( fp, skip, FALSE);
                   }
              }

            num_read++;
/****************************************************/
/*** If requested only go forward max_lines *********/
/****************************************************/
            if ( max_lines > 0 && num_read > max_lines ) 
              {
                done = TRUE;
                printf("find_line checked %d lines returning\n", num_read);
                *p_found= FALSE;
                return;
              }
            if (debug) printf(".....................Testing >>%s<< >>%s<<\n", tok, tok2);

            if ( tok && strcmp(tok,"none") != 0 && strcmp(tok,p_key) == 0 && strcmp(p_key2, "none") != 0) 
              {
                 if (num_read > 0) 
                            for (iloop=0; iloop <= sep; iloop++) tok2= tok_get( fp, noskip, FALSE);

                 done = done || strcmp(p_key2, tok2) == 0;      
                 printf("Read >>%s<< and >>%s<<\n", tok, tok2);
                 if ( done ) printf("Happy with that!\n");

              }
            else if ( tok2 && strcmp(p_key,"none") == 0) 
              {

                 if (num_read > 0) 
                            for (iloop=0; iloop <= sep; iloop++) if (tok2) tok2= tok_get( fp, noskip, FALSE);

                 done = done || (tok2 && strcmp(p_key2, tok2) == 0);      
                 if (tok2) printf("Read >>%s<<\n", tok2);
                 if ( done ) 
                   {
                      printf("Happy with that!\n");
                      *p_found = TRUE;
                      return;
                   }

/*** Check for stopper ***/
                 if (tok2 && strcmp(p_stopper, tok) == 0)
                   {
                      done=TRUE;
                      hit_stopper=TRUE;
                   }

              }
/*** Check for stopper ***/
            else if (tok && strcmp(p_stopper, tok) == 0)
              {
                done=TRUE;
                hit_stopper=TRUE;
              }
            else if ( tok && strcmp(tok,"none") != 0 && strcmp(tok,p_key) == 0) 
              {
                 done = TRUE;
                 *p_found = TRUE;
                 return;
              }

            if ( (tok && strcmp(tok,END_OF_INPUT) == 0) || (tok2 && strcmp(tok2,END_OF_INPUT) == 0) ) 
              {
		 if (!*p_found)
	           {
                     if ( strcmp(p_key2, "none") == 0) 
                       {
                         printf("ERROR the key >>%s<< not found in file\n",
                                    p_key);
                       }
                     else if ( strcmp(p_key, "none") == 0) 
                       {
                         printf("ERROR the key >>%s<< not found in position %d in file\n",
                                    p_key2, sep);
                       }
                     else
                       {
                         if (max_lines > 0)
                           {
                              printf("WARNING the keys >>%s<< and >>%s<< not found %d words apart in file\n",
                                                 p_key, p_key2, sep); 
                              return;
                           }
                         else
                           {
                              printf("ERROR the keys >>%s<< and >>%s<< not found %d words apart in file\n",
                                                 p_key, p_key2, sep); 
                           }
                       }
		   }
		 else
	           {
                      *p_found=FALSE;
		      return;
	           }
                 exit(0);
              }
         }
     if (done && !hit_stopper)
       {
         *p_found=TRUE;
       }
     else
       {
         printf("Hit stopper in find_line\n");
         *p_found=FALSE;
       }
     return;
    }
