/***************************************************************************/
/*** Based on the cartesian components of the lattice vectors  *************/
/*** calculate the reciprocal space vectors and a,b,c,alpha,beta,gamma *****/
/*** Dave Willock May 04                                               *****/
/***************************************************************************/

#include <stdio.h>
#include <math.h>
#include "constants.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

  void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                            double *p_abc )
    {
    double a_cross_b[3], b_cross_c[3], c_cross_a[3];
    double cell_volume, dot;

    printf("\nLattice vectors as read in:\n");
    printf("%10.6f  %10.6f  %10.6f  \n",
                   *p_latt_vec,     *(p_latt_vec+1), *(p_latt_vec+2));
    printf("%10.6f  %10.6f  %10.6f  \n",
                   *(p_latt_vec+3), *(p_latt_vec+4), *(p_latt_vec+5));
    printf("%10.6f  %10.6f  %10.6f  \n",
                   *(p_latt_vec+6), *(p_latt_vec+7), *(p_latt_vec+8));
     
/****** generate reciprocal lattice vectors ********************************/

    vec_cross( p_latt_vec,   p_latt_vec+3, &a_cross_b[0]);
    vec_cross( p_latt_vec+3, p_latt_vec+6, &b_cross_c[0]);
    vec_cross( p_latt_vec+6, p_latt_vec  , &c_cross_a[0]);

    cell_volume= vec_dot(p_latt_vec, &b_cross_c[0]);


    printf("\nCell Volume = %10.6f\n", cell_volume);

    *p_recip_latt_vec    = b_cross_c[0]/cell_volume;
    *(p_recip_latt_vec+1)= b_cross_c[1]/cell_volume;
    *(p_recip_latt_vec+2)= b_cross_c[2]/cell_volume;

    *(p_recip_latt_vec+3)= c_cross_a[0]/cell_volume;
    *(p_recip_latt_vec+4)= c_cross_a[1]/cell_volume;
    *(p_recip_latt_vec+5)= c_cross_a[2]/cell_volume;

    *(p_recip_latt_vec+6)= a_cross_b[0]/cell_volume;
    *(p_recip_latt_vec+7)= a_cross_b[1]/cell_volume;
    *(p_recip_latt_vec+8)= a_cross_b[2]/cell_volume;

    printf("\nCalculated reciprocal space vectors:\n");
    printf("%10.6f %10.6f %10.6f \n",
          *p_recip_latt_vec, *(p_recip_latt_vec+1), *(p_recip_latt_vec+2));
    printf("%10.6f %10.6f %10.6f \n",
          *(p_recip_latt_vec+3), *(p_recip_latt_vec+4), *(p_recip_latt_vec+5));
    printf("%10.6f %10.6f %10.6f \n", 
          *(p_recip_latt_vec+6), *(p_recip_latt_vec+7), *(p_recip_latt_vec+8));
/*******************************************************************/
/*** Cell vector lengths and angles ********************************/
/*******************************************************************/

     *p_abc =    *p_latt_vec     * *p_latt_vec
              +  *(p_latt_vec+1) * *(p_latt_vec+1)
              +  *(p_latt_vec+2) * *(p_latt_vec+2);
     *p_abc = sqrt(*p_abc);
            
     *(p_abc+1) =    *(p_latt_vec+3) * *(p_latt_vec+3)
                  +  *(p_latt_vec+4) * *(p_latt_vec+4)
                  +  *(p_latt_vec+5) * *(p_latt_vec+5);
     *(p_abc+1) = sqrt(*(p_abc+1));
            
     *(p_abc+2) =    *(p_latt_vec+6) * *(p_latt_vec+6)
                  +  *(p_latt_vec+7) * *(p_latt_vec+7)
                  +  *(p_latt_vec+8) * *(p_latt_vec+8);
     *(p_abc+2) = sqrt(*(p_abc+2));

/*** alpha ***/

     dot=  *(p_latt_vec+3) * *(p_latt_vec+6)
         + *(p_latt_vec+4) * *(p_latt_vec+7)
         + *(p_latt_vec+5) * *(p_latt_vec+8);

     dot= dot / (*(p_abc+1) * *(p_abc+2));

     *(p_abc+3) = RAD_TO_DEG*acos(dot);

/*** beta  ***/

     dot=  *p_latt_vec     * *(p_latt_vec+6)
         + *(p_latt_vec+1) * *(p_latt_vec+7)
         + *(p_latt_vec+2) * *(p_latt_vec+8);

     dot= dot / (*p_abc * *(p_abc+2));

     *(p_abc+4) = RAD_TO_DEG*acos(dot);

/*** gamma ***/

     dot=  *p_latt_vec     * *(p_latt_vec+3)
         + *(p_latt_vec+1) * *(p_latt_vec+4)
         + *(p_latt_vec+2) * *(p_latt_vec+5);

     dot= dot / (*p_abc * *(p_abc+1));

     *(p_abc+5) = RAD_TO_DEG*acos(dot);

      printf("\nabc= %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n", 
                    *p_abc, *(p_abc+1), *(p_abc+2), *(p_abc+3), 
                                                    *(p_abc+4), 
                                                    *(p_abc+5));
    return;
   } 
