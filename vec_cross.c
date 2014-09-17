#include <stdio.h>

/* the dot product between two vectors */

void vec_cross(double *p_A, double *p_B, double *p_cross)
{

*p_cross=     ( *(p_A+1) )*( *(p_B+2) )-( *(p_A+2) )*( *(p_B+1) );
*(p_cross+1)= ( *(p_A+2) )*( *p_B     )-(  *p_A    )*( *(p_B+2) );
*(p_cross+2)= ( *p_A     )*( *(p_B+1) )-( *(p_A+1) )*( *p_B     );

  return;

}

