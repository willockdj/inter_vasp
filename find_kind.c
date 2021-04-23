#include <string.h>

#include "maxima.h"
#include "reader.h"
  
/******************************************************************************
  Find_kind takes a parsed token, and matches it against the directive list
  contained in reader.h, returns a BLANK_DIRECT value if token is unrecognised,
  and a UNIT value if the token is extraneous clutter..
******************************************************************************/

list null_namelist[]   = {NULL_DIRECTIVE_LIST};
list first_namelist[]  = {PRIME_DIRECTIVE_LIST};
list second_namelist[] = {SECOND_DIRECTIVE_LIST};
list third_namelist[] = {THIRD_DIRECTIVE_LIST}; 
list siesta_namelist[] = {SIESTA_DIRECTIVE_LIST}; 
list siesta_2nd_namelist[] = {SIESTA_2ND_DIRECTIVE_LIST}; 
list onetep_namelist[] = {ONETEP_DIRECTIVE_LIST}; 
list onetep_2nd_namelist[] = {ONETEP_2ND_DIRECTIVE_LIST}; 

/*****************************************************************************/
/*** find_kind altered to be context aware by knowing the level from *********/
/*** which the token was read. Dave Willock March 1997 ***********************/
/*****************************************************************************/

int find_kind(char *token, int level)
{
  int i,j,k,l,len;
  
  i = 0;
  j = 0;
  k = 0;
  l = 0;

  if (!token) return(-1);

if (level == PRIME_DIRECTIVE || level == REMAINING_DIRECTIVE)
  {
    do
      {
        if (!strncmp(token, first_namelist[j].directive,4))
          {
	     return(first_namelist[j].token_index);
	  }
      }
  while(first_namelist[j++].token_index != BLANK_DIRECT);
  }
else if (level == SECONDARY_DIRECTIVE)
  {
    do
      {
        if (!strncmp(token,second_namelist[k].directive,4))
  	  {
	    return(second_namelist[k].token_index);
	  }
      }
    while(second_namelist[k++].token_index != BLANK_DIRECT); 
  }
else if (level == TERTIARY_DIRECTIVE)
  {
    do 
     { 
/**********************************************************/
/**** For 3rd directive compare only to 2nd character *****/
/**********************************************************/
       if (!strncmp(token,third_namelist[l].directive,2))
	 { 
	  return(third_namelist[l].token_index); 
	 } 
     } 
    while(third_namelist[l++].token_index != BLANK_DIRECT);
  }

/*********************************************************/
/** For SIESTA force exact keyword match (except case) ***/
/*********************************************************/
else if (level == SIESTA_DIRECTIVES)
  {
    do
      {
        len = strlen(siesta_namelist[j].directive);
        if (!strncmp(token, siesta_namelist[j].directive,len))
          {
	     return(siesta_namelist[j].token_index);
	  }
      }
  while(siesta_namelist[j++].token_index != BLANK_DIRECT);
  }
else if (level == SIESTA_2ND_DIRECTIVES)
  {
    do
      {
        len = strlen(siesta_2nd_namelist[j].directive);
        if (!strncmp(token, siesta_2nd_namelist[j].directive,len))
          {
	     return(siesta_2nd_namelist[j].token_index);
	  }
      }
  while(siesta_2nd_namelist[j++].token_index != BLANK_DIRECT);
  }

else if (level == ONETEP_DIRECTIVES)
  {
    do
      {
        len = strlen(onetep_namelist[j].directive);
        if (!strncmp(token, onetep_namelist[j].directive,len))
          {
	     return(onetep_namelist[j].token_index);
	  }
      }
  while(onetep_namelist[j++].token_index != BLANK_DIRECT);
  }
else if (level == ONETEP_2ND_DIRECTIVES)
  {
    do
      {
        len = strlen(onetep_2nd_namelist[j].directive);
        if (!strncmp(token, onetep_2nd_namelist[j].directive,len))
          {
	     return(onetep_2nd_namelist[j].token_index);
	  }
      }
  while(onetep_2nd_namelist[j++].token_index != BLANK_DIRECT);
  }

    do
      {
        if (!strncmp(token, null_namelist[i].directive,2))
	  {
	    return(UNIT);
	  }
      }
    while(null_namelist[i++].token_index != BLANK_DIRECT);

  return(BLANK_DIRECT);
}

