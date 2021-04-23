/***************************************************************************/
/*** Define simple hemi-spherical cluster for ChemShell applications       */
/*** Dave Willock March/April 2018 *****************************************/
/***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include "maxima.h"
#include "structures.h"
#include "constants.h"
#include "global_values.h"
#include "reader.h"

/* protype list for this routine */
void unit_vector(double *p_vector);

double size_vector(double *p_vector);

void vec_cross(double *p_A, double *p_B, double *p_cross);

void latt_vecs_from_cart( double *p_latt_vec, 
                          double *p_recip_latt_vec,
                          double *p_abc );

void rotate_vector(double *p_vec, double *p_axis, double theta);

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

void rotate_with_flags(atom *p_molecule, double *p_axis, double *p_origin,
                       double theta, int *p_flag_list, int num_atoms);

void min_image( double *x, double *y, double *z,
                           double *p_recip_latt_vec, double *p_latt_vec);

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

int pbc_interactions(double *p_dx, double *p_dy, double *p_dz, double cutoff,
                                   double cutoff_2, double *p_pos_separations,  
                                   double *p_pos_vectors, int count_flag,
                                   int is_self, double *p_recip_latt_vec, 
                                   double *p_latt_vec, int need_hemi );

void define_hemi(atom *p_molecule, int *p_num_atoms, double *p_latt, 
                 double *p_recip_latt, double *p_abc, double radius,
                 double *p_origin, atom *p_cluster_mol, int *p_num_cluster_atoms)
{
  int iatom, jatom, image, end_list;
  int num_new,need_rot, flag[MAXATOMS];
  int min_vecs[3], max_vecs[3], ivec, ivec2, ivec3;
  int ilayers, map_top[MAXATOMS], map_bot[MAXATOMS];
  int num_top, num_bot, num_interactions;
  int found, doit, got_layers;
  int imiller, h_miller, k_miller, l_miller;
  int this_miller[3], super[3], just_count;
  int need_hemi, use;

  double dot, theta;
  double x,y,z;
  double vec[3], this_vec[3], frac[3];
  double u,v,w,dhkl, sum_vec[3], q[3], qi[3];
  double trans_a[3], trans_b[3], trans_c[3]; 
  double zmax, zmin, d, radius2;
  double pos_separations[MAX_PERI_IMAGES], pos_vectors[MAX_PERI_IMAGES3]; 

  atom *p_atom, *p_end_atom, *p_top_atom, *p_bot_atom;
  atom *p_this_atom, *p_next_atom, *p_this_clus_atom;

  *p_num_cluster_atoms = -1;
  radius2=radius*radius;
  need_hemi=TRUE;
  just_count=0; // for pbc_interactions this setting means fill the vectors straight away //

/**** define supercell dimensions that will be required to get all of the **/
/**** hemi-sphere atoms.                                                  **/

      p_atom=p_molecule;
      p_this_clus_atom = p_cluster_mol;
      for (iatom=0; iatom <= *p_num_atoms; iatom++)
        {
/*** shift atom to position relative to origin ***/

          x= p_atom->x - *p_origin;
          y= p_atom->y - *(p_origin+1);
          z= p_atom->z - *(p_origin+2);

/*** Get pbc within the cut-off **/

          num_interactions = pbc_interactions(&x, &y, &z, radius, radius2,
                                              &pos_separations[0], &pos_vectors[0], just_count,
                                              FALSE, p_recip_latt, p_latt, need_hemi);

          if (num_interactions < MAX_PERI_IMAGES)
            {
          printf("Atom %d (%s) %10.6f %10.6f %10.6f has %d images in hemi-sphere test djw\n", iatom, p_atom->label,
                                                                                     x,y,z, num_interactions);

// Check this atom will be unique
//

              ivec=0;
              for (image=0; image <= num_interactions; image++)
                 {
                    use = TRUE; 
                    if ( *p_num_cluster_atoms >= 0 )
                      {
                        p_this_atom= p_cluster_mol;
                        for (jatom=0; jatom <= *p_num_cluster_atoms; jatom++)
                          {
                            x=pos_vectors[ivec] - p_this_atom->x;
                            y=pos_vectors[ivec+1] - p_this_atom->y;
                            z=pos_vectors[ivec+2] - p_this_atom->z;

                            if ( x*x + y*y + z*z < 1E-3 ) 
                                {
                                   printf("REJECTING A MATCHED IMAGE\n");
                                   use = FALSE;
                                }
                      
                            p_this_atom++;
              
                          }
                      }

                    if (use)
                      {
                        *p_this_clus_atom   = *p_atom; 
                        p_this_clus_atom->x = pos_vectors[ivec]; 
                        p_this_clus_atom->y = pos_vectors[ivec+1];
                        p_this_clus_atom->z = pos_vectors[ivec+2];
    
                        printf("cluster atom: %s %10.6f  %10.6f  %10.6f \n", p_this_clus_atom->label, 
                                                                             p_this_clus_atom->x,
                                                                             p_this_clus_atom->y,
                                                                             p_this_clus_atom->z);
 
                        ++*p_num_cluster_atoms;
                        p_this_clus_atom++;
                    }
                  ivec+=3;
                 }
            }
          else
            {
              printf("ERROR: Too many periodic images needed when building cluster, decrease radius or increase\n");
              printf("ERROR: MAX_PERI_IMAGES (%d)  and MAX_PERI_IMAGES3 (%d) in maxima.h.\n", 
                                                                                 MAX_PERI_IMAGES, MAX_PERI_IMAGES3);
              exit(0);
            }

          p_atom++;
        }

  return; 
   }

