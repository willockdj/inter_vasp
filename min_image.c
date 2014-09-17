#include <stdio.h>
#include <math.h>
#include "global_values.h"
#include "maxima.h"
#include "structures.h"

/*********************************************************************************/
/********  Routine to work out the minimum image vector for   ********************/
/********  periodic systems                                   ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *frac_a, double *frac_b, double *frac_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

void min_image( double *x, double *y, double *z, double *p_recip_latt_vec, double *p_latt_vec)
  {
double fract_a, fract_b, fract_c;
int done;

/******* find the fractional co-ords of this vector *************/
 
    cart_to_fract( *x, *y, *z, &fract_a, &fract_b, &fract_c, p_recip_latt_vec);

/******* adjust to min image ************************************/
    if (fract_a > 0.5) fract_a = fract_a - floor(fract_a);
    if (fract_b > 0.5) fract_b = fract_b - floor(fract_b);
    if (fract_c > 0.5) fract_c = fract_c - floor(fract_c);

    if (fract_a < -0.5) fract_a = fract_a - ceil(fract_a);
    if (fract_b < -0.5) fract_b = fract_b - ceil(fract_b);
    if (fract_c < -0.5) fract_c = fract_c - ceil(fract_c);

    if (fract_a >  0.5)  fract_a= fract_a - 1.0; 
    if (fract_a < -0.5)  fract_a= fract_a+ 1.0;

    if (fract_b >  0.5)  fract_b= fract_b- 1.0;
    if (fract_b < -0.5)  fract_b= fract_b+ 1.0;

    if (fract_c >  0.5)  fract_c= fract_c- 1.0;
    if (fract_c < -0.5)  fract_c= fract_c+ 1.0;

/****** convert back to cartessian ******************************/

    fract_to_cart( x, y, z, fract_a, fract_b, fract_c, p_latt_vec );

    return;   

  }
