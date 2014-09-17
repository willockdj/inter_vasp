#include <stdio.h>
#include <math.h>

/* routine to give the unit vector in the direction of the supplied vector */

/* Functions required for this program -------------------------*/

double size_vector(double *p_vector);

/*--------------------------------------------------------------*/

void unit_vector(double *p_vector)
{

  int iloop;
  double size, *p_start;

  size = size_vector( p_vector );

  if (size < 0.000001) printf("Warning a vector with magnitude less than 1E-6 is being made unit!!\n");

   for (iloop=0; iloop < 3; iloop++)
    {
         *p_vector= (*p_vector)/size;
         p_vector++;
    }
}

