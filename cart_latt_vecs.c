/*********************************************************************************/
/********  Routine to work out cartesian real and reciprocal  ********************/
/********  lattice vectors from the a,b,c,alpha,beta,gamma    ********************/
/********  line of a biosym .car vector, assuming that the XYZ********************/
/********  convention is followed in the .mdf file.           ********************/
/********                                                     ********************/
/********  p_latt_vec and p_recip_latt_vec point to 9 comp.   ********************/
/********  arrays to be filled                                ********************/
/********  ax ay az bx by bz cx cy cz                         ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

#include <math.h>
#include <stdio.h>
#include "constants.h"

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

void cart_latt_vecs( double *p_abc, double *p_latt_vec, double *p_recip_latt_vec)
  {
    double cell_volume, bx_cx, mag_a, mag_b, mag_c;
    double a_cross_b[3], b_cross_c[3], c_cross_a[3]; 
    double alpha, beta, gamma;
    double dot_set[9];

    int iloop;

/********* make alpha beta gamma in radians ****************************/

    mag_a = *p_abc;
    mag_b = *(p_abc+1);
    mag_c = *(p_abc+2);
    alpha = *(p_abc+3)/RAD_TO_DEG;
    beta  = *(p_abc+4)/RAD_TO_DEG;
    gamma = *(p_abc+5)/RAD_TO_DEG;

/********* X indicates that a is to be aligned with the x axis *********/

    *p_latt_vec    = mag_a;
    *(p_latt_vec+1)= 0.0;
    *(p_latt_vec+2)= 0.0;

/********* XY indicates that the b axis is to lie in the plane formed ******/
/********* by the y-axis and the x-axis i.e. perp. to Z               ******/

    *(p_latt_vec+3)= mag_b * cos( gamma);
    *(p_latt_vec+4)= mag_b * sin( gamma);
    *(p_latt_vec+5)= 0.0;

/******** XYZ indicates that the c vector position is now defined **********/

    *(p_latt_vec+6)= mag_c * cos( beta); 

    bx_cx = ( *(p_latt_vec+3)) * ( *(p_latt_vec+6));

    *(p_latt_vec+7)= ( mag_b*mag_c*cos( alpha)  - bx_cx )/ ( *(p_latt_vec+4));

    *(p_latt_vec+8)= sqrt ( mag_c*mag_c - *(p_latt_vec+6) * *(p_latt_vec+6) 
                                        - *(p_latt_vec+7) * *(p_latt_vec+7)); 

/****** Print out lattice vectors ******************************************/
  
    printf ("Lattice vectors : \n");
    printf ("a: %10.6f %10.6f %10.6f\n", *p_latt_vec    , *(p_latt_vec+1), *(p_latt_vec+2));
    printf ("b: %10.6f %10.6f %10.6f\n", *(p_latt_vec+3), *(p_latt_vec+4), *(p_latt_vec+5));
    printf ("c: %10.6f %10.6f %10.6f\n", *(p_latt_vec+6), *(p_latt_vec+7), *(p_latt_vec+8));

/****** generate reciprocal lattice vectors ********************************/

    vec_cross( p_latt_vec,   p_latt_vec+3, &a_cross_b[0]); 
    vec_cross( p_latt_vec+3, p_latt_vec+6, &b_cross_c[0]); 
    vec_cross( p_latt_vec+6, p_latt_vec  , &c_cross_a[0]); 

    cell_volume= vec_dot(p_latt_vec, &b_cross_c[0]); 

    printf(" a cross b = %10.6f %10.6f %10.6f \n", a_cross_b[0],a_cross_b[1],a_cross_b[2]);
    printf(" b cross c = %10.6f %10.6f %10.6f \n", b_cross_c[0],b_cross_c[1],b_cross_c[2]);
    printf(" c cross a = %10.6f %10.6f %10.6f \n", c_cross_a[0],c_cross_a[1],c_cross_a[2]);

    printf("Cell Volume = %10.6f\n", cell_volume);

    *p_recip_latt_vec    = b_cross_c[0]/cell_volume;
    *(p_recip_latt_vec+1)= b_cross_c[1]/cell_volume;
    *(p_recip_latt_vec+2)= b_cross_c[2]/cell_volume;

    *(p_recip_latt_vec+3)= c_cross_a[0]/cell_volume;
    *(p_recip_latt_vec+4)= c_cross_a[1]/cell_volume;
    *(p_recip_latt_vec+5)= c_cross_a[2]/cell_volume;

    *(p_recip_latt_vec+6)= a_cross_b[0]/cell_volume;
    *(p_recip_latt_vec+7)= a_cross_b[1]/cell_volume;
    *(p_recip_latt_vec+8)= a_cross_b[2]/cell_volume;

/**************** Check by forming all dot products this should *******/
/**************** give the identity matrix ****************************/

    dot_set[0] = vec_dot(p_recip_latt_vec, p_latt_vec);    
    dot_set[1] = vec_dot(p_recip_latt_vec+3, p_latt_vec);    
    dot_set[2] = vec_dot(p_recip_latt_vec+6, p_latt_vec);    

    dot_set[3] = vec_dot(p_recip_latt_vec, p_latt_vec+3);    
    dot_set[4] = vec_dot(p_recip_latt_vec+3, p_latt_vec+3);    
    dot_set[5] = vec_dot(p_recip_latt_vec+6, p_latt_vec+3);    

    dot_set[6] = vec_dot(p_recip_latt_vec, p_latt_vec+6);    
    dot_set[7] = vec_dot(p_recip_latt_vec+3, p_latt_vec+6);    
    dot_set[8] = vec_dot(p_recip_latt_vec+6, p_latt_vec+6);    


       printf ("Check on recip/real lat vecs\n\n");

       for (iloop = 0; iloop < 9; iloop++)
         {
            printf("%10.6f  ", dot_set[iloop]);
            if (iloop == 2 || iloop == 5 || iloop == 8) printf("\n");
         }


    return;
  }
