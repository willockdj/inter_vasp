#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"
#include "data.h"

void write_doscsv(FILE *fp, char *p_title_x, char *p_title_y, 
                  char *p_title_z,
                  dos *p_dos,
                  int have_tot, int num, int is_restricted, double fermi)
  {
int i;

double val;

printf("Fermi Energy = %f\n Correct?\n",fermi);

fprintf(fp,"Fermi Energy (eV), %f\n",fermi);

fprintf(fp,"\n");

if (is_restricted)
  {
    fprintf(fp,"%s corrected by Fermi Energy, %s",p_title_x,p_title_y); 

    if ( have_tot )
        fprintf(fp,", %s\n",p_title_z); 

    fprintf(fp,"\n"); 
  }
else 
  {
    fprintf(fp,"%s corrected by Fermi Energy, up_%s, down_%s",p_title_x,p_title_y,p_title_y);

    if ( have_tot )
        fprintf(fp,", up_%s, down_%s\n",p_title_z, p_title_z);

    fprintf(fp,"\n");
  }


for (i=0; i<=num; i++)
  {

    val= p_dos->energy;
    if (fabs(val) < 1.0e-4 || fabs(val) >1.0e4 )
      {
        fprintf(fp,"%e",val-fermi); 
      }
    else
      {
        fprintf(fp,"%f",val-fermi); 
      }
    fprintf(fp,","); 

    val= p_dos->up_dos;
    if ( fabs(val) < 1.0e-4 || fabs(val) > 1.0e4 )
      {
        fprintf(fp,"%e", val); 
      }
    else
      {
        fprintf(fp,"%f", val); 
      }
    
    if (!is_restricted)
      {
        fprintf(fp,",");
        val= -p_dos->down_dos;
        if ( fabs(val) < 1.0e-4 || fabs(val) > 1.0e4 )
          {
            fprintf(fp,"%e", val);
          }
        else
          {
            fprintf(fp,"%f", val);
          }
      }

    if ( have_tot )
      {
        fprintf(fp,","); 
        val= p_dos->up_totdos;

        if (fabs(val) < 1.0e-4 || fabs(val) > 1.0e4 )
          {
            fprintf(fp,"%e", val); 
          }
        else
          {
            fprintf(fp,"%f", val); 
          }

        if (!is_restricted)
          {  
            fprintf(fp,",");
            val= -p_dos->down_totdos;
  
            if (fabs(val) < 1.0e-4 || fabs(val) > 1.0e4 )
              {
                fprintf(fp,"%e", val);
              }
            else
              {
                fprintf(fp,"%f", val);
              }
          }
      }
    fprintf(fp,"\n"); 
    p_dos++;
  }

return;
  }

