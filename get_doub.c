#include <stdio.h>
#include <math.h>
#include <limits.h>
#include "global_values.h"
long int get_long_int( int *p_ichar, int *point_j,
                               int *itsanum, int *ndigi,int i, int *sign);

/* routine to pick out signed double from garbage */

double get_doub( int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum )
{

long int idave;
long int iafter_dec,sign_predec;
int itsanum,ndigi,iloop,istart, sign;
double fdave,fafter_dec;

printf("In get_doub num_of_chars %d, place %d and itsanum %d\n", num_of_chars, *p_place, *p_itsanum);

*p_itsanum=0;

/* step through line looking for negative signs and numbers */

while ( *p_place <= num_of_chars && !*p_itsanum )
    {
         idave= get_long_int( p_ichar, p_place, p_itsanum, &ndigi, num_of_chars, 
                                                                    &sign);

         printf("get_doub read int %ld with sign %d\n", idave, sign);
/* if next is not a '.' we have an integer or a pure decimal with
   no leading zero eg. .5670  */

         if (*p_itsanum && *(p_ichar+*p_place) != '.') 
          {
          istart= *p_place -ndigi;
          if (*(p_ichar+istart-1) == '.')
	    {
             fafter_dec= labs(idave);
             for (iloop=0; iloop < ndigi; ++iloop) fafter_dec= 0.1*fafter_dec;
             fdave= fafter_dec;
            }
          else fdave=sign*idave;
     	  }

/* otherwise its a double  */
         else if (*p_itsanum)
         {
         sign_predec= sign;
        
         iafter_dec = get_long_int(p_ichar, p_place, p_itsanum, &ndigi, 
                                                  num_of_chars, &sign);

         printf("get_doub read int after decimal %ld with %d digits\n", iafter_dec, ndigi);
         fafter_dec= labs(iafter_dec);
         printf("Converts to %10.6f\n",fafter_dec);

         for (iloop=0; iloop < ndigi; ++iloop) fafter_dec= 0.1*fafter_dec;
         printf("after scaling %10.6f\n",fafter_dec);
         fdave= sign_predec * (labs(idave)+fafter_dec);
         }

    }
if ( !*p_itsanum) fdave=0.0;
return fdave;
}









