/*************************************************************************/
/*** This routine gives the integer/set of integers corresponding to the */
/*** state in the distrib array. Dave Willock, July 2013.                */
/*************************************************************************/
#include <stdio.h>

void detail_state(int *p_distrib, int *p_state, int num_sites)
   {
     int iloop, signif;
     int num_state_ints, icount;
     int istate;
     int *p_place;
     int *p_state_digi;
     int int_bits;

/************************************************************************/
/*** The accuracy with which intergers are held determines how many *****/
/*** intergers are needed to index this state                       *****/
/************************************************************************/

     int_bits= 32;
     int_bits--;

     signif=1;
     icount=0;
     
     for (iloop=0; iloop < num_sites/int_bits + 2; iloop++ ) *(p_state+iloop)= 0;

     p_place= p_distrib;
     p_state_digi= p_state;
     num_state_ints = 0;
/************************************************************/
/** Bug fix August 2000 allow highest sig. bit set **********/
/** num_sites zero reference! DJW                  **********/
/************************************************************/
     for (iloop=0; iloop <= num_sites; iloop++)
       {
         icount++;
         if (icount== int_bits)
           {
             icount=0;
             p_state_digi++;
             signif=1;
             num_state_ints++;
            
           }

         if (*p_place==1) *p_state_digi += signif;
         p_place++;
         signif *=2;
       }

   return;
   }

