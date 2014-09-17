#include <stdio.h>

/* the dot product between two vectors */

double vec_dot(double *p_A, double *p_B)
{

  int iloop;
  double dot;

  dot=0.0;
  for (iloop=0; iloop < 3; iloop++)
                        dot+= (*(p_B+iloop))*(*(p_A+iloop));

  return dot;

}

