#include <stdio.h>
#include <math.h>
#include "maxima.h"
#include "constants.h"
#include "structures.h"

void write_pdb(FILE *pdb_fp, atom *p_molecule, double *p_abc, int num_atoms,
               int *p_super, double *p_latt_vec)
{
int iloop, iatom, ineigh, neigh_index;
int ia, ib, ic, image;
char elem[3];
atom *p_atom;
double dx, dy, dz, dist;
double image_vec[3];
               
fprintf(pdb_fp,"CRYST1  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f P1\n", *p_abc,
                     *(p_abc+1),*(p_abc+2),*(p_abc+3),
                                           *(p_abc+4),
                                           *(p_abc+5));

/***********************************************************/
/*** p_super points to an integer array defining the *******/
/*** extent of the required super cell               *******/
/***********************************************************/

for ( ia = 0; ia <= *p_super-1; ia++ )
 {

   for ( ib = 0; ib <= *(p_super+1)-1 ; ib++ )
    {

      for ( ic = 0; ic <= *(p_super+2)-1; ic++ )
        {
          image_vec[0] =  ia * *p_latt_vec 
                         +ib * *(p_latt_vec+3)
                         +ic * *(p_latt_vec+6);

          image_vec[1] =  ia * *(p_latt_vec+1)
                         +ib * *(p_latt_vec+4)
                         +ic * *(p_latt_vec+7);

          image_vec[2] =  ia * *(p_latt_vec+2)
                         +ib * *(p_latt_vec+5)
                         +ic * *(p_latt_vec+8);

          p_atom = p_molecule;


          for ( iloop=0; iloop < num_atoms; iloop++)
            {
    printf("%10.6f %10.6f\n", p_atom->z , image_vec[2]);

    fprintf(pdb_fp,"HETATM %4d  %2s        %4d  %7.3f %7.3f %7.3f\n",
                                        iloop+1, p_atom->elem, iloop+1,
                                        p_atom->x + image_vec[0],
                                        p_atom->y + image_vec[1], 
                                        p_atom->z + image_vec[2]);

              p_atom++;
            }
        }
    }
 }

image=0;
for ( ia = 0; ia <= *p_super-1; ia++ )
 {

   for ( ib = 0; ib <= *(p_super+1)-1 ; ib++ )
    {

      for ( ic = 0; ic <= *(p_super+2)-1; ic++ )
        {
           p_atom = p_molecule;
           for ( iloop=0; iloop < num_atoms; iloop++)
             {
               fprintf(pdb_fp,"CONECT %4d", iloop+1+image*num_atoms);
     
               for ( ineigh=0; ineigh < p_atom->num_neigh; ineigh++)
                 {
                   neigh_index = p_atom->neighb[ineigh]; 

                   dx = p_atom->x - (p_molecule+neigh_index)->x;
                   dy = p_atom->y - (p_molecule+neigh_index)->y;
                   dz = p_atom->z - (p_molecule+neigh_index)->z;

                   dist = sqrt( dx*dx + dy*dy + dz * dz);

                   if (dist < 3.0) fprintf(pdb_fp," %4d", 
                                   neigh_index+1+image*num_atoms);
                 }

               fprintf(pdb_fp,"\n");
               p_atom++;
            }
          image++;
        }
     }
  }

fprintf(pdb_fp,"END\n");

return;
}

