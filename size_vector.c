#include <stdio.h>
#include <math.h>

/* routine to give the size of a vector */

double size_vector(double *p_vector)
{

  int iloop;
  double size;

  size = 0.0;
  for (iloop=0; iloop < 3; iloop++)
         size= size + (*(p_vector+iloop)) * (*(p_vector+iloop));

  return sqrt(size);

}

