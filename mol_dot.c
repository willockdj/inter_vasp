/***************************************************************************/
/*** Take molecule dot product, i.e. the products of the set of x,y,z ******/
/*** coords summed.                                                   ******/
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

/******************************************************************/
/** Note that num_images is just the number input on images *******/
/** directive so is the proper number of images.            *******/
/******************************************************************/

double mol_dot(atom *p_mol1, atom *p_mol2, int num_atoms)
{
  int iatom;

  double dot;

  dot= 0.0;
  for (iatom=0; iatom < num_atoms; iatom++)
    {

       dot+= p_mol1->x * p_mol2->x
            +p_mol1->y * p_mol2->y
            +p_mol1->z * p_mol2->z;

       p_mol1++;
       p_mol2++;
     }

  return dot;
}
