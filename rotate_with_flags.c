/****************************************************************************/
/******* routine to rotate a set of cartesian co-ordinates around     *******/
/******* a given vector, with the set distributed in an array so that *******/
/******* those with rot_flag = 1 (TRUE) get rotated those with        *******/ 
/******* rot_flag = 0 (FALSE) do not                                  *******/
/*******                                                              *******/
/******* Dave Willock and Dewi Lewis May 1995                         *******/
/****************************************************************************/

#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"


/* required functions ------------------------------------*/

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

/*--------------------------------------------------------*/

void rotate_with_flags(atom *p_molecule, double *p_axis, double *p_origin,
                       double theta, int *p_flag_list, int num_atoms)
{

  int iloop,index,icomp,icoloum,irow;
  int *p_flag;

  double *p_origin_comp;
  double theta_on_2, e[4], quarternion[9], product[10], hold_coord[3];

  atom *p_atom;

   p_origin_comp= p_origin;
  *p_origin_comp= (*p_origin_comp) * -1.0;
   p_origin_comp++;
  *p_origin_comp= (*p_origin_comp) * -1.0;
   p_origin_comp++;
  *p_origin_comp= (*p_origin_comp) * -1.0;

  move_molecule( p_molecule, num_atoms, p_origin); 

/* Build up the quarternion matrix */

  theta_on_2 = theta/2.0;

  e[0] = cos(theta_on_2); 
  e[1] = (*p_axis)* sin(theta_on_2);
  p_axis++;
  e[2] = (*p_axis)* sin(theta_on_2);
  p_axis++;
  e[3] = (*p_axis)* sin(theta_on_2);

/* arrange quarternion elements :
 *
 *     0  1  2
 *     3  4  5
 *     6  7  8
 *
 * products set as
 *
 *    0 = e0e0,  1 = e1e1,  2 = e2e2,  3 = e3e3,
 *    4 = e0e1,  5 = e0e2,  6 = e0e3,
 *    7 = e1e2,  8 = e1e3,
 *    9 = e2e3
 */

     for (iloop=0; iloop <= 3; iloop++) product[iloop]= e[iloop]*e[iloop];
     for (iloop=4; iloop <= 6; iloop++) product[iloop]= e[0]*e[iloop-3];
     for (iloop=7; iloop <= 8; iloop++) product[iloop]= e[1]*e[iloop-5];
     product[9]= e[2]*e[3];

/*********** diagonals *********************************************/

     quarternion[0]= product[0] + product[1] - product[2] - product[3];
     quarternion[4]= product[0] - product[1] + product[2] - product[3];
     quarternion[8]= product[0] - product[1] - product[2] + product[3];

/************ upper right off diagonals ****************************/

     quarternion[1]= 2.0*(product[7] - product[6]);
     quarternion[2]= 2.0*(product[8] + product[5]);
     quarternion[5]= 2.0*(product[9] - product[4]);

/************ lower left off diagonals *****************************/

     quarternion[3]= 2.0*(product[7] + product[6]);
     quarternion[6]= 2.0*(product[8] - product[5]);
     quarternion[7]= 2.0*(product[9] + product[4]);

/************ perform the rotation on flagged atoms only ***********/

  p_flag= p_flag_list-1;
  for (iloop= 0; iloop <= num_atoms; iloop++)
  {
    p_flag++;
    if ( *p_flag )
      {
         p_atom= p_molecule+iloop; 
         for (icoloum=0; icoloum< 3; icoloum++)
           {
              hold_coord[icoloum]= 0.0;
              hold_coord[icoloum]+= quarternion[3*icoloum]*(p_atom->x);
              hold_coord[icoloum]+= quarternion[3*icoloum+1]*(p_atom->y);
              hold_coord[icoloum]+= quarternion[3*icoloum+2]*(p_atom->z);
           } 
         p_atom->x = hold_coord[0];
         p_atom->y = hold_coord[1];
         p_atom->z = hold_coord[2];
      }
  }

/*********** move atom set back to origin's position **************/

   p_origin_comp= p_origin;
  *p_origin_comp= (*p_origin_comp) * -1.0;
   p_origin_comp++;
  *p_origin_comp= (*p_origin_comp) * -1.0;
   p_origin_comp++;
  *p_origin_comp= (*p_origin_comp) * -1.0;

   move_molecule( p_molecule, num_atoms, p_origin);

 return;
}

