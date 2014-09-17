#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "maxima.h"
#include "structures.h"
#include "header.h"
#include "global_values.h"
#include "reader.h"

/******************************************************************************
  tok_get.c : parses command lines in a semi intelligent fashion
IMportant change 22/6/95 DWL Now ignores EVERYTHING after a #

Added lower_case logical to allow atom names to remain in upper case after
reading (DJW June 98)
******************************************************************************/

int find_kind(char *token, int level);


char * tok_get(FILE *input_fp, int skip_lines, int lower_case)
{
  int          i, force_read;
  char         *tok;
  char         *comment;	

  i = 0;
  read_new_line = FALSE;

/**************************************************************************/
/*** If the last token read was a NULL character last_tok will act ********/
/*** like a FALSE logical and force us to read the next line   ************/
/*** If we call for data from the next line get it anyway      ************/
/*** AND in this if block altered to OR to do this June 98 DJW ************/
/**************************************************************************/
/*  printf("Arrived in tok_get with last_tok = >>%s<<\n", last_tok); */

/*  if (skip_lines==TRUE) printf("tok_get can skip lines\n");     */
/*                    else printf("tok_get cannot skip lines\n"); */

  if (!last_tok || skip_lines)
    {

/**************************************************************************/
/**** Read in the next line from input_fp, if it is NULL return NULL ******/
/**************************************************************************/

      line = fgets(target, BUFFER, input_fp);

/*    printf("tok_get reads line >>%s<<\n", line); */

      if (line == NULL ) return(END_OF_INPUT);
 
/**************************************************************************/
/**** Remove commented lines or parts of lines that are comments **********/
/**************************************************************************/

      if ((comment = strchr(line,'#')) != NULL) 
        {

/**************************************************************************/
/****** chop off the bits after the hash **********************************/
/**************************************************************************/

           *(comment) = '\0';

/**************************************************************************/
/****** if this is a simple comment line it will now have length 0 ********/
/****** so simply return a single # character                      ********/
/**************************************************************************/

           if (strlen(line) < 1) 
             {
               *line= '#';
               *(line+1)= '\0';
               return(line);
             }
         }

/**************************************************************************/
/***** read_new_line says we have just read a new line in !! **************/
/**************************************************************************/

       read_new_line = TRUE;

      /****** dont bother counting lines ****/
       /**** line_no++; ****/
    }
  else if (!last_tok && !skip_lines)
    {
       return(NULL);
    }

/**************************************************************************/
/****** If we have just read a new line copy it to the buffer *************/
/**************************************************************************/

  if (read_new_line) strncpy(buf, line, BUFFER);

/**************************************************************************/
/***** Now process the latest string of information ***********************/
/**************************************************************************/

  for(;;)
    {

/**************************************************************************/
/***** tok is read from the last line if some tokens are left *************/
/**************************************************************************/

    tok = (read_new_line) ? strtok(buf, " :;,\n\r\t") : strtok(NULL, " :;,\n\r\t");

      read_new_line = (!tok) ? TRUE : FALSE; 
      last_tok = tok;
      if (!tok || !isalpha(*tok)) break;

      if (lower_case) for (i = 0; tok[i] != '\0'; i++) tok[i] = tolower(tok[i]);

      if (find_kind(tok, REMAINING_DIRECTIVE) != UNIT) break;
      continue;
    }
  return(tok);
}
