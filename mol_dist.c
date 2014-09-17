/***************************************************************************/
/*** Work out total distance between structures.                      ******/
/*** Dave Willock December 05 **********************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */

/*---------------------------------------------------------------------------*/

double mol_dist(atom *p_mol1, atom *p_mol2, int num_atoms)
{
  int iatom;

  double dist, dx,dy,dz;

  dist= 0.0;
  for (iatom=0; iatom < num_atoms; iatom++)
    {

       dx = p_mol1->x - p_mol2->x;
       dy = p_mol1->y - p_mol2->y;
       dz = p_mol1->z - p_mol2->z;

       dist+=  dx*dx + dy*dy + dz*dz;

       p_mol1++;
       p_mol2++;
     }

  return sqrt(dist);
}
