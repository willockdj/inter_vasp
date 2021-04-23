#include <stdio.h>

void set_binary(int *p_distrib, int state, int num_sites)
   {
     int iloop, signif;
     int num_state_ints, icount;
     int istate;
     int *p_place;
     int int_bits;

/************************************************************************/
/*** The accuracy with which intergers are held determines how many *****/
/*** intergers are needed to index this site                        *****/
/************************************************************************/

     int_bits= 32;
     int_bits--;

     signif=1;
     
     p_place= p_distrib;
     num_state_ints = 0;
     for (iloop=0; iloop <= num_sites; iloop++)
       {
         *p_place=0;
         if ( state & signif ) *p_place=1 ;

         p_place++;
         signif *=2;
       }

   return;
   }

