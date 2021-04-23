/***************************************************************************/
/*** Read in a block of forces from VASP OUTCAR file ***********************/
/*** Dave Willock May 04 ***************************************************/
/*** updated for mallocing: num_atoms is proper counter not highest index **/
/*** Dave Willock Nov 17 ***************************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

/*-------------------------------*/

void read_force_block( FILE *fp, e_vec *p_forces, int num_atoms)
{
  int iloop, jloop, iatom;
  int noskip=FALSE, skip=TRUE;

  char *tok;
/******** skip dotted lines *****************************/

  tok= tok_get( fp, skip, FALSE);

  printf("Reading forces for %d atoms\n", num_atoms);
  printf("throwing >>%s<<\n",tok);
  for (iatom = 0; iatom < num_atoms; iatom++)
    {
      tok_get( fp, skip, FALSE);
      for ( jloop = 0; jloop <= 1; jloop++ ) tok= tok_get( fp, noskip, FALSE);

      p_forces->dx[iatom] = atof( tok_get( fp, noskip, FALSE ));
      p_forces->dy[iatom] = atof( tok_get( fp, noskip, FALSE ));
      p_forces->dz[iatom] = atof( tok_get( fp, noskip, FALSE ));

      printf("Just read %10.6f %10.6f %10.6f\n", p_forces->dx[iatom]
                                               , p_forces->dy[iatom]
                                               , p_forces->dz[iatom]);
    }
  return;
}
