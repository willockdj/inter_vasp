/***************************************************/
/***** Random number generator to give a      ******/
/***** random double 0.0 -> 1.0 using the     ******/
/***** standard library routines random and   ******/
/***** srandom with optional seeding on the   ******/
/***** current calander time                  ******/
/*****                                        ******/
/***** Dave Willock May 95                    ******/
/***************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double real_random(int done)
{
long irand;
double real_rand;
time_t time_now;

/**** seed if done = 0 **************************************/

if (!done)
   {
      time_now= time(&time_now);

      if (time_now == -1)
        {
          printf("Warning: Calander time not available seeding with 1\n");
          time_now = 1;
        }
         
     srand((unsigned int) time_now ); 
   }

    irand = rand();
    real_rand = irand;

return real_rand/RAND_MAX;
}
