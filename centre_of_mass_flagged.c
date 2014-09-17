#include <stdio.h>
#include "maxima.h"
#include "structures.h"

/* protype list for this routine */

/*---------------------------------------------------------------------------*/


/* work out the centre of mass for a molecule in .car in insight formatted file */

void centre_of_mass_flagged(double *p_c_of_m, double *p_total_mass, 
                            atom *p_molecule, 
                            int *p_chosen_indices, int num_chosen )
{
  int icomp, this_atom;
  atom *p_atom;

  double atomic_mass;
/*----------------------------------------------------------------------------*/

   *p_c_of_m= 0.0;
   *(p_c_of_m+1)= 0.0; 
   *(p_c_of_m+2)= 0.0;


   *p_total_mass= 0.0;

   for (this_atom=0; this_atom <= num_chosen; this_atom++)
   {
       p_atom= p_molecule+ *(p_chosen_indices+this_atom);
       atomic_mass = p_atom->mass;

/*       printf("Atom %s has mass %10.6f\n", p_atom->label, atomic_mass); */

       *p_c_of_m     += atomic_mass * (p_atom->x);
       *(p_c_of_m+1) += atomic_mass * (p_atom->y);
       *(p_c_of_m+2) += atomic_mass * (p_atom->z);

       *p_total_mass += atomic_mass;

   }

/* normalise to total mass */

    if (*p_total_mass > 0.0)
      {
        for (icomp=0; icomp < 3; icomp++)
                        *(p_c_of_m+icomp)= *(p_c_of_m+icomp) / *p_total_mass;
      }
    else
      {
        printf("ERROR: Total mass of molecule is zero when attempting to find centre of mass\n");
        printf("ERROR: Do not use moment of inertia option unless atom masses are set\n");
        exit(0);
      }

}
