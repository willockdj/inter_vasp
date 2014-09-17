#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "global_values.h"
int get_int( int *p_ichar, int *point_j,
                               int *itsanum, int *ndigi,int i, int *sign);

/* routine to pick out signed double from garbage */

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum )
{

int idave,sign,num_mmbrs;
int itsanum,iafter_dec,ndigi,iloop,istart,sign_predec;
double fdave,fafter_dec;

*p_itsanum=0;

/* step through line looking for negative signs and numbers */

while ( *p_place <= num_of_chars && !*p_itsanum )
    {
         idave= get_int( p_ichar, p_place, p_itsanum, &ndigi, num_of_chars, 
                                                                    &sign);

/* if next is not a '.' we have an integer or a pure decimal with
   no leading zero eg. .5670  */

         if (*p_itsanum && *(p_ichar+*p_place) != '.') 
          {
          istart= *p_place -ndigi;
          if (*(p_ichar+istart-1) == '.')
	    {
             fafter_dec= fabs(idave);
             for (iloop=0; iloop < ndigi; ++iloop) fafter_dec= 0.1*fafter_dec;
             fdave= fafter_dec;
            }
          else fdave=sign*idave;
     	  }

/* otherwise its a double  */
         else if (*p_itsanum)
         {
         sign_predec= sign;
        
         iafter_dec = get_int(p_ichar , p_place, p_itsanum, &ndigi, 
                                                  num_of_chars, &sign);

         fafter_dec= fabs(iafter_dec);

         for (iloop=0; iloop < ndigi; ++iloop) fafter_dec= 0.1*fafter_dec;
         fdave= sign_predec * (fabs(idave)+fafter_dec);
         }

    }
if ( !*p_itsanum) fdave=0.0;
return fdave;
}









