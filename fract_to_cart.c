#include <math.h>

/*********************************************************************************/
/********  Routine go from fractional to cartessian co-ords   ********************/
/********  for a single 3-component vector                    ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z, 
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec )
  {

     *cart_x =   *p_latt_vec     * frac_a
               + *(p_latt_vec+3) * frac_b
               + *(p_latt_vec+6) * frac_c;
 
     *cart_y =   *(p_latt_vec+1) * frac_a
               + *(p_latt_vec+4) * frac_b
               + *(p_latt_vec+7) * frac_c;
 
     *cart_z =   *(p_latt_vec+2) * frac_a
               + *(p_latt_vec+5) * frac_b
               + *(p_latt_vec+8) * frac_c;

     return;
   }
