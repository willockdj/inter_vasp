#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"
#include "math.h"

void write_csv(FILE *fp, char *p_title_x, char *p_title_y, 
               char *p_title_z,
               double *p_x, double *p_y, double *p_z,
               int have_y, int have_z, int num)
  {
int i;

printf("write_csv expecting %d as highest index of a list with %d members.\n", num, num+1);

fprintf(fp,"Index, %s", p_title_x); 

if (have_y)
        fprintf(fp,",%s", p_title_y); 

if ( have_z )
        fprintf(fp,",%s", p_title_z); 

fprintf(fp,"\n"); 


for (i=0; i<=num; i++)
  {
    fprintf(fp,"%d",i+1); 
    
    if (fabs(*p_x) < 1.0e-4 || fabs(*p_x) >1.0e4 )
      {
        fprintf(fp,",%e",*p_x); 
      }
    else
      {
        fprintf(fp,",%f",*p_x); 
      }

    if (have_y)
      {
        if ( fabs(*p_y) < 1.0e-4 || fabs(*p_y) > 1.0e4 )
          {
            fprintf(fp,",%e", *p_y); 
          }
        else
          {
            fprintf(fp,",%f", *p_y); 
          }
        p_y++;
      }

    if ( have_z )
      {
        if (fabs(*p_z) < 1.0e-4 || fabs(*p_z) > 1.0e4 )
          {
            fprintf(fp,",%e", *p_z); 
          }
        else
          {
            fprintf(fp,",%f", *p_z); 
          }
        p_z++;
      }
    fprintf(fp,"\n"); 
    p_x++;
  }

return;
  }

