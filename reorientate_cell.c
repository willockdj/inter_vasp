/***************************************************************************/
/*** Reorientate cell to make c perpendicular to the chosen Miller surface */
/*** Dave Willock March/April 2007 *****************************************/
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

double vec_dot(double *p_A, double *p_B);

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

void reorientate_cell(atom *p_molecule, int *p_num_atoms, double *p_latt, 
                      double *p_recip_latt, double *p_abc, int *p_miller, int num_miller,
                      double *p_new_latt, double *p_new_recip_latt, 
                      double *p_new_abc, atom *p_slab_mol, int *p_num_slab_atoms)
{
  int iatom, jatom, iloop, end_list;
  int num_new, need_rot1, need_rot, flag[MAXATOMS];
  int min_vecs[3], max_vecs[3], ivec1, ivec2, ivec3;
  int ilayers, map_top[MAXATOMS], map_bot[MAXATOMS];
  int num_top, num_bot, imiller, h_miller, k_miller, l_miller;
  int found, doit, got_layers, need_faces_same, max_layers;

  double dot, theta, theta1, siz;
  double surf_norm[3], norm[3];
  double old_latt[9], old_recip_latt[9], old_abc[6];
  double axis[3], axis1[3], origin[3], cross_check[3];
  double vec[3], this_vec[3], frac[3];
  double u,v,w,dhkl, sum_vec[3];
  double trans_a[3], trans_b[3], trans_c[3]; 
  double zmax, zmin, new_miller[MAX_MILLER3];
  double hmil, kmil, lmil;

  atom *p_atom, *p_end_atom, *p_top_atom, *p_bot_atom;
  atom *p_this_atom, *p_next_atom;

/**** Assume that the Miller indices define a real space vector which is ***/
/**** perpendicular to the plane                                         ***/

  surf_norm[0] = *p_miller *(*p_recip_latt)    + *(p_miller+1)* *(p_recip_latt+3) + *(p_miller+2)* *(p_recip_latt+6);
  surf_norm[1] = *p_miller * *(p_recip_latt+1) + *(p_miller+1)* *(p_recip_latt+4) + *(p_miller+2)* *(p_recip_latt+7);
  surf_norm[2] = *p_miller * *(p_recip_latt+2) + *(p_miller+1)* *(p_recip_latt+5) + *(p_miller+2)* *(p_recip_latt+8);

  siz=size_vector(&surf_norm[0]);

  unit_vector(&surf_norm[0]);

  printf("Unit normal to surface of interest: %10.6f, %10.6f %10.6f\n", surf_norm[0], 
                                                                        surf_norm[1], 
                                                                        surf_norm[2]);

  printf("Original size: %10.6f\n", siz);

  hmil = siz*vec_dot(&surf_norm[0], p_latt);
  kmil = siz*vec_dot(&surf_norm[0], p_latt+3);
  lmil = siz*vec_dot(&surf_norm[0], p_latt+6);

  printf("Recovered Miller indices          : %10.6f, %10.6f %10.6f\n", hmil, kmil, lmil);

/** Work out angle between normal and c-vector ***/

  dot = surf_norm[0]* *(p_latt+6) + surf_norm[1]* *(p_latt+7) + surf_norm[2]* *(p_latt+8); 

  dot = dot / size_vector((p_latt+6));

  theta1 = acos(dot);

  printf("Angle between Miller plane normal and c-vector is: %10.6f degrees\n", theta1*RAD_TO_DEG);

  need_rot1 = TRUE;
  if (theta1 < 1.0E-4)
    {
      printf("This cell already has c-vector perpendicular to surface, no rotation required\n");
      need_rot1 = FALSE;
      axis1[0] = 0.0;
      axis1[1] = 0.0;
      axis1[2] = 1.0;
    }
  else
    {
      vec_cross(&surf_norm[0], (p_latt+6), &axis1[0]);

      unit_vector(&axis1[0]);
    }

/*** Calculate surface vectors as difference between intercepts from Miller ***/
/*** Indices                                                                ***/
/*** Catch special cases first                                              ***/

      printf("Old Lattice vectors: \n");
      printf("\n %10.6f %10.6f %10.6f Len: %10.6f\n", *p_latt    , *(p_latt+1), *(p_latt+2), size_vector(p_latt));
      printf(" %10.6f %10.6f %10.6f Len: %10.6f\n",   *(p_latt+3), *(p_latt+4), *(p_latt+5), size_vector(p_latt+3));
      printf(" %10.6f %10.6f %10.6f Len: %10.6f\n",   *(p_latt+6), *(p_latt+7), *(p_latt+8), size_vector(p_latt+6));


  if (*p_miller == 0 && *(p_miller+1) == 0 && *(p_miller+2) == 0)
    {
       printf("ERROR: All three Miller indices zero in reorientate.\n");
       exit(0);
    }
  else if (*p_miller == 0 && *(p_miller+1) == 0)
    {
/*** This plane is the a-b plane ****/
       *p_new_latt = *p_latt;
       *(p_new_latt+1) = *(p_latt+1);
       *(p_new_latt+2) = *(p_latt+2);

       *(p_new_latt+3) = *(p_latt+3);
       *(p_new_latt+4) = *(p_latt+4);
       *(p_new_latt+5) = *(p_latt+5);

    }
  else if (*(p_miller+1) == 0 && *(p_miller+2) == 0)
    {
/*** This plane is the b-c plane        ***/
/*** a replaces old b, b replaces old c ***/
       *p_new_latt = *(p_latt+3);
       *(p_new_latt+1) = *(p_latt+4);
       *(p_new_latt+2) = *(p_latt+5);

       *(p_new_latt+3) = *(p_latt+6);
       *(p_new_latt+4) = *(p_latt+7);
       *(p_new_latt+5) = *(p_latt+8);

    }
  else if (*p_miller == 0 && *(p_miller+2) == 0)
    { 
/*** This plane is the a-c plane        ***/
/*** a replaces old c, b replaces old a ***/
       *p_new_latt     = *(p_latt+6);
       *(p_new_latt+1) = *(p_latt+7);
       *(p_new_latt+2) = *(p_latt+8);

       *(p_new_latt+3) = *p_latt;
       *(p_new_latt+4) = *(p_latt+1);
       *(p_new_latt+5) = *(p_latt+2);
    }
  else if (*p_miller == 0)
    {
/*** This plane is parallel to a-vector ***/
/*** because h is zero                  ***/
/*** So only the old a-vector is in the ***/
/*** plane                              ***/

       *p_new_latt = *p_latt;
       *(p_new_latt+1) = *(p_latt+1);
       *(p_new_latt+2) = *(p_latt+2);

/*** for new b-vector use difference between b and c intercepts ***/
       v = 1.0/(double) *(p_miller+1);
       w = 1.0/(double) *(p_miller+2);

       *(p_new_latt+3) = v* *(p_latt+3) - w*  *(p_latt+6);
       *(p_new_latt+4) = v* *(p_latt+4) - w*  *(p_latt+7);
       *(p_new_latt+5) = v* *(p_latt+5) - w*  *(p_latt+8);

    }
  else if (*(p_miller+1) == 0)
    {
/*** This plane is parallel to b-vector ***/
/*** because k is zero                  ***/
/*** So only the old b-vector is in the ***/
/*** plane                              ***/

       *(p_new_latt+3) = *(p_latt+3);
       *(p_new_latt+4) = *(p_latt+4);
       *(p_new_latt+5) = *(p_latt+5);

/*** for new a-vector use difference between a and c intercepts ***/
       u = 1.0/(double) *p_miller;
       w = 1.0/(double) *(p_miller+2);

       *p_new_latt     = u* *(p_latt)   - w*  *(p_latt+6);
       *(p_new_latt+1) = u* *(p_latt+1) - w*  *(p_latt+7);
       *(p_new_latt+2) = u* *(p_latt+2) - w*  *(p_latt+8);

    }
  else if (*(p_miller+2) == 0)
    {
/*** This plane is parallel to c-vector ***/
/*** because l is zero                  ***/
/*** So only the old c-vector is in the ***/
/*** plane                              ***/

       printf("This plane is parallel to the c-vector\n");

/*** for new a-vector use the c-vector ****************************/
       *p_new_latt = *(p_latt+6);
       *(p_new_latt+1) = *(p_latt+7);
       *(p_new_latt+2) = *(p_latt+8);

/*** for new b-vector use vector between a and b intercepts ***/
       u = 1.0/(double) *p_miller;
       v = 1.0/(double) *(p_miller+1);

       *(p_new_latt+3) =  u*  *p_latt     - v* *(p_latt+3);
       *(p_new_latt+4) =  u*  *(p_latt+1) - v* *(p_latt+4);
       *(p_new_latt+5) =  u*  *(p_latt+2) - v* *(p_latt+5);

    }
  else 
    {
       printf("The general case %d %d %d\n",*p_miller, *(p_miller+1), *(p_miller+2));
/*** This plane is parallel to none of the current axes ***/
       u = 1.0/(double) *p_miller;
       v = 1.0/(double) *(p_miller+1);
       w = 1.0/(double) *(p_miller+2);

       printf("Intercepts %10.6f %10.6f %10.6f\n",u,v,w);

/*** for new a-vector use vector between c and a intercept ****/
       *p_new_latt     =  u * *p_latt    - w * *(p_latt+6);
       *(p_new_latt+1) =  u * *(p_latt+1)- w * *(p_latt+7);
       *(p_new_latt+2) =  u * *(p_latt+2)- w * *(p_latt+8);

/*** for new b-vector use vector between c and b intercepts ***/
       *(p_new_latt+3) =  v * *(p_latt+3)- w * *(p_latt+6);
       *(p_new_latt+4) =  v * *(p_latt+4)- w * *(p_latt+7);
       *(p_new_latt+5) =  v * *(p_latt+5)- w * *(p_latt+8);
    }

/*** Set c-vector to be perpendicular to the surface and with length dhkl ***/
/*** Work out dhkl as 1/ | h a* + k b* + l c* |                           ***/

      sum_vec[0] = *p_miller * *p_recip_latt     + *(p_miller+1) * *(p_recip_latt+3) + *(p_miller+2) * *(p_recip_latt+6);
      sum_vec[1] = *p_miller * *(p_recip_latt+1) + *(p_miller+1) * *(p_recip_latt+4) + *(p_miller+2) * *(p_recip_latt+7);
      sum_vec[2] = *p_miller * *(p_recip_latt+2) + *(p_miller+1) * *(p_recip_latt+5) + *(p_miller+2) * *(p_recip_latt+8);

      dhkl = 1.0/size_vector(&sum_vec[0]);

      printf("dhkl= %10.6f\n", dhkl);

      *(p_new_latt+6) = dhkl*surf_norm[0];
      *(p_new_latt+7) = dhkl*surf_norm[1];
      *(p_new_latt+8) = dhkl*surf_norm[2];

      printf("New Lattice vectors: \n");
      printf("\n %10.6f %10.6f %10.6f Len: %10.6f\n", *p_new_latt, *(p_new_latt+1), 
                                                      *(p_new_latt+2), size_vector(p_new_latt));

      printf(" %10.6f %10.6f %10.6f Len: %10.6f\n",   *(p_new_latt+3), *(p_new_latt+4), 
                                                      *(p_new_latt+5), size_vector(p_new_latt+3));

      printf(" %10.6f %10.6f %10.6f Len: %10.6f\n",   *(p_new_latt+6), *(p_new_latt+7),
                                                      *(p_new_latt+8), size_vector(p_new_latt+6));

/*** Cross check a cross b is same direction as c ***/
      vec_cross(p_new_latt, (p_new_latt+3), &cross_check[0]);

      printf("a x b = %10.6f %10.6f %10.6f\n", cross_check[0],  cross_check[1], cross_check[2]);

      for (iloop=0; iloop < 9; iloop++) old_latt[iloop] = *(p_latt+iloop);

/*** Make sure all co-ordinates are min_images in original cell ***/
      for (iatom=0; iatom < *p_num_atoms; iatom++)
        {
          *(p_slab_mol+iatom)= *(p_molecule+iatom);

          vec[0] = (p_slab_mol+iatom)->x;
          vec[1] = (p_slab_mol+iatom)->y;
          vec[2] = (p_slab_mol+iatom)->z;

          min_image( &vec[0], &vec[1], &vec[2], p_recip_latt, p_latt);

          (p_slab_mol+iatom)->x = vec[0];
          (p_slab_mol+iatom)->y = vec[1];
          (p_slab_mol+iatom)->z = vec[2];

          flag[iatom]=TRUE;
        }


/*** Rotate the new lattice vectors and the old to the new orientation                   ***/

      if (need_rot1)
        {
          printf("Using rotation axis: %10.6f, %10.6f %10.6f angle %10.6f\n", 
                                                 axis1[0],  axis1[1],  axis1[2], theta1*RAD_TO_DEG);

          rotate_vector(&old_latt[0], &axis1[0], theta1);
          rotate_vector(&old_latt[3], &axis1[0], theta1);
          rotate_vector(&old_latt[6], &axis1[0], theta1);

          rotate_vector(p_new_latt,   &axis1[0], theta1);
          rotate_vector(p_new_latt+3, &axis1[0], theta1);
          rotate_vector(p_new_latt+6, &axis1[0], theta1);

          latt_vecs_from_cart( p_new_latt, p_new_recip_latt, p_new_abc ); 

          printf("Old Lattice vectors after rotation: \n");
          printf("\n %10.6f %10.6f %10.6f\n", old_latt[0], old_latt[1], old_latt[2]);
          printf(" %10.6f %10.6f %10.6f\n",   old_latt[3], old_latt[4], old_latt[5]);
          printf(" %10.6f %10.6f %10.6f\n",   old_latt[6], old_latt[7], old_latt[8]);

          printf("\nNew Lattice vectors after rotation: \n");
          printf("\n %10.6f %10.6f %10.6f\n", *p_new_latt    , *(p_new_latt+1), *(p_new_latt+2));
          printf(" %10.6f %10.6f %10.6f\n",   *(p_new_latt+3), *(p_new_latt+4), *(p_new_latt+5));
          printf(" %10.6f %10.6f %10.6f\n",   *(p_new_latt+6), *(p_new_latt+7), *(p_new_latt+8));
       }
     else
       {
          printf("New axes already have c-vector out of surface\n");
       }

/**** Rotate atoms by same amount ***/
     if (need_rot1)
       {
         origin[0] = 0.0;
         origin[1] = 0.0;
         origin[2] = 0.0;

         rotate_with_flags(p_slab_mol, &axis1[0], &origin[0],
                           theta1, &flag[0], *p_num_atoms);   
       }


/*** To fit with XYZ convention now align the new a with the cartessian X direction ***/
/*** and do the same rotation on the old lattice vectors                            ***/

/** Work out angle between X-axis and a-vector ***/

  dot = *p_new_latt / size_vector(p_new_latt);

  theta = acos(dot);

  printf("Angle between new a vector and cartessian X is: %10.6f\n", theta*RAD_TO_DEG);

  need_rot = TRUE;
  if (theta < 1.0E-4)
    {
      printf("This cell already has a-vector aligned with X\n");
      need_rot = FALSE;
    }

  if (need_rot)
    {
       axis[0] = 0.0;
       axis[1] = 0.0;

/*** Cross product of a and x must lie along Z and so to find out which way to rotate **/
/*** use the ay component ( as the only survivor in the cross product for the axis    **/

       if (*(p_new_latt+1) < 0.0)
         {
           axis[2] = 1.0;
         }
       else
         {
           axis[2] = -1.0;
         }

/*** Rotate the new lattice vectors and the old to the new orientation                   ***/
       printf("Old Lattice vectors before rotation: \n");
       printf("\n %10.6f %10.6f %10.6f len: %10.6f\n", old_latt[0], old_latt[1], old_latt[2], 
                                                                            size_vector(&old_latt[0]));
       printf(" %10.6f %10.6f %10.6f len: %10.6f\n",   old_latt[3], old_latt[4], old_latt[5], 
                                                                            size_vector(&old_latt[3]));
       printf(" %10.6f %10.6f %10.6f len: %10.6f\n",   old_latt[6], old_latt[7], old_latt[8], 
                                                                            size_vector(&old_latt[6]));

       printf("Using rotation axis: %10.6f, %10.6f %10.6f\n", axis[0],  axis[1],  axis[2]);

       rotate_vector(&old_latt[0], &axis[0], theta);
       rotate_vector(&old_latt[3], &axis[0], theta);
       rotate_vector(&old_latt[6], &axis[0], theta);

       rotate_vector(p_new_latt,   &axis[0], theta);
       rotate_vector(p_new_latt+3, &axis[0], theta);
       rotate_vector(p_new_latt+6, &axis[0], theta);

       printf("New Lattice vectors after rotation: \n");
       printf("\n %10.6f %10.6f %10.6f len: %10.6f\n", *p_new_latt, *(p_new_latt+1), *(p_new_latt+2), 
                                                                            size_vector(p_new_latt));

       printf(" %10.6f %10.6f %10.6f len: %10.6f\n",   *(p_new_latt+3), *(p_new_latt+4), *(p_new_latt+5), 
                                                                            size_vector(p_new_latt+3));

       printf(" %10.6f %10.6f %10.6f len: %10.6f\n",   *(p_new_latt+6), *(p_new_latt+7), *(p_new_latt+8), 
                                                                            size_vector(p_new_latt+6));

/*** If there are a larger set of miller planes rotate the normals and reform ****/

       if (num_miller > 0)
         {

          for (imiller=1; imiller< num_miller; imiller++)
            {
               h_miller= 3*imiller; k_miller= h_miller+1; l_miller= k_miller+1;

               norm[0]=  *p_recip_latt     * *(p_miller+h_miller) 
                       + *(p_recip_latt+3) * *(p_miller+k_miller) 
                       + *(p_recip_latt+6) * *(p_miller+l_miller);

               norm[1]=  *(p_recip_latt+1) * *(p_miller+h_miller) 
                       + *(p_recip_latt+4) * *(p_miller+k_miller) 
                       + *(p_recip_latt+7) * *(p_miller+l_miller);

               norm[2]=  *(p_recip_latt+2) * *(p_miller+h_miller) 
                       + *(p_recip_latt+5) * *(p_miller+k_miller) 
                       + *(p_recip_latt+8) * *(p_miller+l_miller);

               printf("Miller set %d in old orientation of %d %d %d\n", imiller, *(p_miller+h_miller), 
                                                                                 *(p_miller+k_miller), 
                                                                                 *(p_miller+l_miller));

               printf("Miller set %d gives normal vector in old orientation of %10.6f %10.6f %10.6f\n",
                                                                                      imiller, norm[0], norm[1], norm[2]);
               if (need_rot1) rotate_vector(&norm[0], &axis1[0], theta1);
               if (need_rot)  rotate_vector(&norm[0], &axis[0],  theta );

/*** norm is now a vector in the new cartessian space ****/

               new_miller[h_miller-3]= *p_new_latt     * norm[0] + *(p_new_latt+1) * norm[1] + *(p_new_latt+2) * norm[2];
               new_miller[k_miller-3]= *(p_new_latt+3) * norm[0] + *(p_new_latt+4) * norm[1] + *(p_new_latt+5) * norm[2];
               new_miller[l_miller-3]= *(p_new_latt+6) * norm[0] + *(p_new_latt+7) * norm[1] + *(p_new_latt+8) * norm[2];

               printf("Miller set %d gives normal vector in new orientation of %10.6f %10.6f %10.6f\n",
                                                                                      imiller, norm[0], norm[1], norm[2]);

               printf("Miller indices set %d reorientated to : %3.1f %3.1f %3.1f\n",
                                imiller, new_miller[h_miller-3], new_miller[k_miller-3], new_miller[l_miller-3]);

            }
        }

/*** Rotate atoms to final orientation ************************/
       origin[0] = 0.0; origin[1] = 0.0; origin[2] = 0.0;

       rotate_with_flags(p_slab_mol, &axis[0], &origin[0],
                         theta, &flag[0], *p_num_atoms); 

       printf("Old Lattice vectors after rotation: \n");
       printf("\n %10.6f %10.6f %10.6f len: %10.6f\n", old_latt[0], old_latt[1], old_latt[2], 
                                                                            size_vector(&old_latt[0]));
       printf(" %10.6f %10.6f %10.6f len: %10.6f\n",   old_latt[3], old_latt[4], old_latt[5], 
                                                                            size_vector(&old_latt[3]));
       printf(" %10.6f %10.6f %10.6f len: %10.6f\n",   old_latt[6], old_latt[7], old_latt[8], 
                                                                            size_vector(&old_latt[6]));

       printf("\nNew Lattice vectors after rotation: \n");
       printf("\n %10.6f %10.6f %10.6f\n", *p_new_latt    , *(p_new_latt+1), *(p_new_latt+2));
       printf(" %10.6f %10.6f %10.6f\n",   *(p_new_latt+3), *(p_new_latt+4), *(p_new_latt+5));
       printf(" %10.6f %10.6f %10.6f\n",   *(p_new_latt+6), *(p_new_latt+7), *(p_new_latt+8));
    }
  else
    {
       printf("New axes already have a-vector along X\n");
    }

/*** Get reciprocal lattice vectors for old system reorientated ***/

     latt_vecs_from_cart( &old_latt[0], &old_recip_latt[0], &old_abc[0] ); 

/**** loop to search for a slab with equivalent faces                     ***/
/**** The top and the bottom of the dhkl slab need not be equivalent this ***/
/**** loop allows additional layers to be added until they are.           ***/

  need_faces_same=FALSE;

  if (need_faces_same)
   {
     max_layers= 10;
   }
  else
   {
     max_layers= 1;
   }

  *(p_new_latt+8) = 0.0;
  for (ilayers=0; ilayers<max_layers; ilayers++)
   {

     *(p_new_latt+8) = *(p_new_latt+8) + dhkl;

/*** Work out new recip lattice vectors and appropriate abc for the new axes ****/

     latt_vecs_from_cart( p_new_latt, p_new_recip_latt, p_new_abc ); 

     printf("\na=%10.6f b=%10.6f c=%10.6f alpha=%10.6f beta=%10.6f gamma=%10.6f\n",
                    *p_new_abc    , *(p_new_abc+1), *(p_new_abc+2),
                    *(p_new_abc+3), *(p_new_abc+4), *(p_new_abc+5));

/*** Fill the new unit cell checking the required number of periodic images of each atom ***/
/*** First need to find out how many old vectors are needed to visit all parts of new cell */

/*** Find the fractional co-ordinates of the new vectors in the old system ***/

     for (ivec1=0; ivec1 < 3; ivec1++)
       {
         cart_to_fract(  *(p_new_latt+ivec1*3), *(p_new_latt+ivec1*3+1), *(p_new_latt+ivec1*3+2),
                        &frac[0], &frac[1], &frac[2], &old_recip_latt[0] );

         printf("New vector %d has fractional co-ords in old system of: %10.6f  %10.6f  %10.6f\n",
                                        ivec1, frac[0], frac[1], frac[2]);

         if (ivec1==0)
           {
             min_vecs[0] = floor(frac[0]);
             min_vecs[1] = floor(frac[1]);
             min_vecs[2] = floor(frac[2]);

             max_vecs[0] = ceil(frac[0]);
             max_vecs[1] = ceil(frac[1]);
             max_vecs[2] = ceil(frac[2]);
          }

         if (frac[0] < min_vecs[0]) min_vecs[0] = floor(frac[0]);
         if (frac[0] > max_vecs[0]) max_vecs[0] = ceil(frac[0]);

         if (frac[1] < min_vecs[1]) min_vecs[1] = floor(frac[1]);
         if (frac[1] > max_vecs[1]) max_vecs[1] = ceil(frac[1]);

         if (frac[2] < min_vecs[2]) min_vecs[2] = floor(frac[2]);
         if (frac[2] > max_vecs[2]) max_vecs[2] = ceil(frac[2]);
       }

     if (min_vecs[0] >= max_vecs[0])
       {
          printf("ERROR: problem defining extent of new cell in reorientate_cell\n");
          printf("       difficulty is in old a-vector direction.\n");
          exit(0);
       }
     if (min_vecs[1] >= max_vecs[1])
       {
          printf("ERROR: problem defining extent of new cell in reorientate_cell\n");
          printf("       difficulty is in old a-vector direction.\n");
          exit(0);
       }
     if (min_vecs[2] >= max_vecs[2])
       {
          printf("ERROR: problem defining extent of new cell in reorientate_cell\n");
          printf("       difficulty is in old a-vector direction.\n");
          exit(0);
       }

     if (fabs(min_vecs[0]) > max_vecs[0]) max_vecs[0] = fabs(min_vecs[0]);
     if (fabs(min_vecs[1]) > max_vecs[1]) max_vecs[1] = fabs(min_vecs[1]);
     if (fabs(min_vecs[2]) > max_vecs[2]) max_vecs[2] = fabs(min_vecs[2]);

     max_vecs[0]++;
     max_vecs[1]++;
     max_vecs[2]++;

     min_vecs[0] = -max_vecs[0];
     min_vecs[1] = -max_vecs[1];
     min_vecs[2] = -max_vecs[2];

     printf("Will fill new cell by scanning old a from %d to %d\n", min_vecs[0], max_vecs[0]);
     printf("                               old b from %d to %d\n", min_vecs[1], max_vecs[1]);
     printf("                           and old c from %d to %d\n", min_vecs[2], max_vecs[2]);

/*** Loop over the required translations and decide if each atom is in new cell ****/

     p_atom=p_slab_mol;
     end_list=*p_num_atoms;

     for (iatom=0; iatom < *p_num_atoms; iatom++)
       {
         flag[iatom]=FALSE;
         printf("\nProcessing atom %d %s coords: %10.6f %10.6f %10.6f\n",iatom, p_atom->label,
                                                                       p_atom->x,
                                                                       p_atom->y,
                                                                       p_atom->z);

         for (ivec1= min_vecs[0]; ivec1 <=max_vecs[0] ; ivec1++)
           {
             trans_a[0] = ivec1*old_latt[0]; 
             trans_a[1] = ivec1*old_latt[1]; 
             trans_a[2] = ivec1*old_latt[2]; 

             for (ivec2= min_vecs[1]; ivec2 <=max_vecs[1] ; ivec2++)
               {
                 trans_b[0] = ivec2*old_latt[3]; 
                 trans_b[1] = ivec2*old_latt[4]; 
                 trans_b[2] = ivec2*old_latt[5]; 

                 for (ivec3= min_vecs[2]; ivec3 <=max_vecs[2] ; ivec3++)
                   {
                     trans_c[0] = ivec3*old_latt[6]; 
                     trans_c[1] = ivec3*old_latt[7]; 
                     trans_c[2] = ivec3*old_latt[8]; 

/*** Work out position of atom with this translation ****/                  

                     vec[0] = p_atom->x +  trans_a[0] + trans_b[0] + trans_c[0];
                     vec[1] = p_atom->y +  trans_a[1] + trans_b[1] + trans_c[1];
                     vec[2] = p_atom->z +  trans_a[2] + trans_b[2] + trans_c[2];

                     printf("%d %d %d Translated vector %10.6f %10.6f %10.6f\n", ivec1, ivec2, ivec3,
                                                                                 vec[0], vec[1], vec[2]);

/*** Find out if this atom is in the new cell ****/
                     cart_to_fract( vec[0], vec[1], vec[2],
                                    &frac[0], &frac[1], &frac[2], p_new_recip_latt );

                     printf("Fractional in new system : %10.6f %10.6f %10.6f\n",  frac[0], frac[1], frac[2]);

                     if (   frac[0] <= 1.00001 && frac[0] >= -0.00001
                         && frac[1] <= 1.00001 && frac[1] >= -0.00001
                         && frac[2] <= 1.00001 && frac[2] >= -0.00001)
                       {
                          printf("Accepted\n");
/**** if this is the atom itself flag to keep it ***/
                          if (ivec1 == 0 && ivec2 == 0 && ivec3 == 0)
                            {
                               flag[iatom]=TRUE;
                            }
/**** if this is an image add it to the end of the list */
                          else
                            {
                               if (end_list < MAXATOMS)
                                 {
                                   p_end_atom = p_slab_mol+end_list;
                                   *p_end_atom = *p_atom;

/**** shift atoms that have shown up on 1.0 boundaries to 0.0 boundary **/
              
                                   doit=FALSE;
                                   if (frac[0] == 1.0)
                                     {
                                       frac[0] = 0.0;
                                       doit=TRUE;
                                     }
                                   if (frac[1] == 1.0)
                                     {
                                       frac[1] = 0.0;
                                       doit=TRUE;
                                     }
                                   if (frac[2] == 1.0)
                                     {
                                       frac[2] = 0.0;
                                       doit=TRUE;
                                     }

                                   doit=FALSE;
                                   if (doit)
                                     {
                                       fract_to_cart( &vec[0], &vec[1], &vec[2], 
                                                      frac[0], frac[1], frac[2],
                                                      p_new_latt );

                                     }

                                   flag[end_list]=TRUE;
                                   p_end_atom->x = vec[0];
                                   p_end_atom->y = vec[1];
                                   p_end_atom->z = vec[2];

                                   printf("New atom coords %10.6f %10.6f %10.6f\n", p_end_atom->x, p_end_atom->y, p_end_atom->z);
                  
                                   end_list++;
                                 }
                               else
                                 {
                                   printf("ERROR: Dimension of atoms list is too short for this slab\n");
                                   exit(0);
                                 }
                               
                            }
                       }
                     printf("\n");
                   }
               }
           }
          p_atom++;
        }

     *p_num_slab_atoms=end_list;

/*** Flag any duplicate atoms for removal ***/
     p_atom = p_slab_mol;
     for (iatom=0; iatom < *p_num_slab_atoms; iatom++)
       {
          if (flag[iatom])
            {
              vec[0]=p_atom->x;
              vec[1]=p_atom->y;
              vec[2]=p_atom->z;

              p_this_atom=p_slab_mol+iatom+1;
              for(jatom=iatom+1; jatom < *p_num_slab_atoms; jatom++)
                {
                  if (flag[jatom])
                   {
                    this_vec[0]=p_this_atom->x;
                    this_vec[1]=p_this_atom->y;
                    this_vec[2]=p_this_atom->z;

                    this_vec[0] -=vec[0];
                    this_vec[1] -=vec[1];
                    this_vec[2] -=vec[2];

                    min_image( &this_vec[0], &this_vec[1], &this_vec[2], 
                                           p_new_recip_latt, p_new_latt); 

                    if (size_vector(&this_vec[0]) < 1E-4) flag[jatom]=FALSE;
                   }

                  p_this_atom++;
                }
             }
           p_atom++;
         }

     printf("After translations have %d atoms in slab list\n\n", *p_num_slab_atoms);
     p_atom = p_slab_mol;
     for (iatom=0; iatom < *p_num_slab_atoms; iatom++)
       {     
          vec[0]=p_atom->x;
          vec[1]=p_atom->y;
          vec[2]=p_atom->z;

          cart_to_fract( vec[0], vec[1], vec[2],
                                    &frac[0], &frac[1], &frac[2], p_new_recip_latt );

          if (flag[iatom] ==TRUE)
             printf("%s frac: %10.6f %10.6f %10.6f flagged TRUE\n", p_atom->label,
                                                              frac[0],  
                                                              frac[1],  
                                                              frac[2]); 
          else
             printf("%s frac: %10.6f %10.6f %10.6f flagged FALSE\n", p_atom->label,
                                                              frac[0],  
                                                              frac[1],  
                                                              frac[2]); 
          p_atom++;
       }

/*** See if the two faces in the c-direction are the same ***/
     p_atom = p_slab_mol;

     zmax=0.0;
     zmin=*(p_new_latt+8);
     for (iatom=0; iatom < *p_num_slab_atoms; iatom++)
       {
         if (iatom >= *p_num_atoms || flag[iatom])
           {
             if (p_atom->z > zmax) zmax = p_atom->z; 
             if (p_atom->z < zmin) zmin = p_atom->z; 
           }

         p_atom++;
       }

     printf("Initial top/bot location finds:\n");
     printf("Top z = %10.6f\n", zmax);
     printf("Bot z = %10.6f\n", zmin);

     if (zmin >= zmax)
       {
         printf("ERROR: cannot define top and bottom of slab.\n");
         exit(0);
       }
     else if (zmin > 0.0001 || zmin < -0.0001)
       {
/*** If no atom is on boundary shift system ***/

         vec[0]=0.0;
         vec[1]=0.0;
         vec[2]= -zmin;

         move_molecule(p_slab_mol, *p_num_slab_atoms, &vec[0]);
   
         zmax -= zmin;
         zmin=0.0;
       }

     printf("After shifting slab in Z:\n");
     printf("Top z = %10.6f\n", zmax);
     printf("Bot z = %10.6f\n", zmin);

/*** identify all atoms that are top and bottom ***/
     num_top=0;
     num_bot=0;
     p_atom = p_slab_mol;
     for (iatom=0; iatom < *p_num_slab_atoms; iatom++)
       {
         if (iatom >= *p_num_atoms || flag[iatom])
           {
             printf("Looking at atom %s z:%10.6f\n",p_atom->label, p_atom->z);
             if (p_atom->z > zmax-0.0001 && p_atom->z < zmax+0.0001)
              {
                 printf("top atom: %s z=%10.6f\n",p_atom->label, p_atom->z);
                 map_top[num_top]=iatom;
                 num_top++;
              }
            if (p_atom->z > zmin-0.0001 && p_atom->z < zmin+0.0001)
              {
                 printf("bot atom: %s z=%10.6f\n",p_atom->label, p_atom->z);
                 map_bot[num_bot]=iatom;
                 num_bot++;
              }
          }
         p_atom++;
       }

printf("\nMap top:\n");
for (iloop=0; iloop<num_top; iloop++)
  {
    printf("%d  :  %d\n",iloop,map_top[iloop]);
  }
printf("\nMap bot:\n");
for (iloop=0; iloop<num_bot; iloop++)
  {
    printf("%d  :  %d\n",iloop,map_bot[iloop]);
  }
printf("\n");

/*** Check if top and bottom faces are identical ****/
     found=-1;
     got_layers=FALSE;
     if (num_top == num_bot)
       {
         found=0;
         for (iatom=0; iatom<num_top; iatom++)
           {
             p_top_atom = p_slab_mol+map_top[iatom];

             printf("\nComparing top atom %s %10.6f  %10.6f  %10.6f\n\n", p_top_atom->label,
                                                                          p_top_atom->x,
                                                                          p_top_atom->y,
                                                                          p_top_atom->z);

             for (jatom=0; jatom<num_bot; jatom++)
               {
                 p_bot_atom = p_slab_mol+map_bot[jatom];

             printf("with   bottom atom %s %10.6f  %10.6f  %10.6f\n\n", p_bot_atom->label,
                                                                          p_bot_atom->x,
                                                                          p_bot_atom->y,
                                                                          p_bot_atom->z);

/*** If these match the top should be a displaced version of the bottom by only a c-vector **/

                 this_vec[0] = p_top_atom->x - p_bot_atom->x; 
                 this_vec[1] = p_top_atom->y - p_bot_atom->y; 
                 this_vec[2] = p_top_atom->z - p_bot_atom->z - *(p_new_latt+8); 

                 printf("this_vec= %10.6f %10.6f %10.6f\n", this_vec[0], this_vec[1], this_vec[2]);

                 if (size_vector(&this_vec[0]) < 1E-4)
                   {
                     found++; 
                     printf("Matched!!! found = %d\n",found);
                   }
               }
           }
       }
     else
       {
         printf("Not good slab, c-faces differ in number of atoms\n");
       }
    if (found == num_top) got_layers=TRUE; 
    if (got_layers) break;
  }


if (!got_layers)
  {
    printf("ERROR: not found slab for this miller plane with equivalent faces\n");
    printf("DEBUG: Continuing anyway.............\n\n");
//    exit(0);
  }

/*** remove atoms flagged FALSE, i.e. original atoms that are not in new cell ***/

     num_new=0;
     for (iatom=0; iatom < *p_num_slab_atoms; iatom++)
       {     
          p_atom = p_slab_mol+num_new;
          printf("Testing %s\n",p_atom->label);
          if (flag[num_new] == FALSE)
            { 
              printf("Removing this one\n");

/*** If this is the last atom in the list just remove it ***/
              if (iatom != (*p_num_slab_atoms)-1)  
                {             
/*** Collapse list on top of this atom ***/
                   p_this_atom = p_atom;  
                   p_next_atom = p_atom+1; 
                   for (jatom=num_new; jatom < (*p_num_slab_atoms)-1; jatom++)
                    {    
                       *p_this_atom = *p_next_atom;
                       p_this_atom++;     
                       p_next_atom++;      
               
                       flag[jatom] = flag[jatom+1];
                    }
                }
            }
          else
            {
              num_new++;
            }
       }

     *p_num_slab_atoms=num_new;

     printf("After flag removal have %d atoms in slab list\n", *p_num_slab_atoms);
     p_atom = p_slab_mol;
     for (iatom=0; iatom < *p_num_slab_atoms; iatom++)
       {     
             printf("%s %10.6f %10.6f %10.6f\n", p_atom->label,
                                                              p_atom->x,
                                                              p_atom->y,
                                                              p_atom->z);
          p_atom++;
       }

  return;
}
