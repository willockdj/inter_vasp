#include <math.h>

/*********************************************************************************/
/********  Routine go from cartessian to fractional co-ords   ********************/
/********  for a single 3-component vector                    ********************/
/********                                                     ********************/
/********  Dave Willock May 1995                              ********************/
/*********************************************************************************/

void cart_to_fract( double cart_x,  double cart_y,  double cart_z, 
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec )
  {
 
     *fract_a =   *p_recip_latt_vec    * cart_x
               + *(p_recip_latt_vec+1) * cart_y
               + *(p_recip_latt_vec+2) * cart_z;
 
     *fract_b =  *(p_recip_latt_vec+3) * cart_x
               + *(p_recip_latt_vec+4) * cart_y
               + *(p_recip_latt_vec+5) * cart_z;
 
     *fract_c =  *(p_recip_latt_vec+6) * cart_x
               + *(p_recip_latt_vec+7) * cart_y
               + *(p_recip_latt_vec+8) * cart_z;

     return;
   }
