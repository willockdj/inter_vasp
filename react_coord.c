/***************************************************************************/
/*** Analyse a set of files for reaction co-ordinate plotting **************/
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

void open_file(FILE **p_file, char *p_filename, char *p_status);

int read_car( FILE *fp, int *p_header_line, int *p_title_line, 
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms, 
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number, 
              double *p_abc, int *p_been_before, group_lists *p_dummy_group, 
              int have_grp, int *p_num_in_group,
              int find_fixed, coord_flags *p_fix_flags);

double atom_separation_squared(atom *p_A, atom *p_B, int pbc,  
                           double *p_recip_latt_vec, double *p_latt_vec);

void min_image( double *x, double *y, double *z, 
                           double *p_recip_latt_vec, double *p_latt_vec);

double mol_dot(atom *p_mol1, atom *p_mol2, int num_atoms);

double mol_dist(atom *p_mol1, atom *p_mol2, int num_atoms);
/*---------------------------------------------------------------------------*/

/******************************************************************/
/** Note that num_images is just the number input on images *******/
/** directive so is the proper number of images.            *******/
/******************************************************************/

void react_coords(atom *p_start_molecule, atom *p_end_molecule, int num_atoms,
                  char_list *p_image_files, int num_images, 
                  double *p_recip_latt_vec, double *p_latt_vec, int pbc)
{
  int good_read, iloop, jloop, skip, iatom;
  int header_line[LINESIZ], title_line[LINESIZ], date_line[LINESIZ];
  int im_num_atoms, im_num_of_mols, im_num_mol_members[MAXMOL];
  int im_mol_number[MAXMOL], been, im_group_indices[MAXATOMS];
  int im_num_in_group, im_have_group, im_pbc;

  double d, d2, dx, dy, dz, dot, mag1, mag2;
  double im_abc[6], tot_d, theta, ab, ap, bp, a2, b2, bot_start, bot_end; 
  double amount_start, amount_end;

  atom image[MAXATOMS], linear_comb_mol[MAXATOMS];
  atom *p_start_atom, *p_end_atom;

  coord_flags dummy_flags;

  group_lists dummy_group;

  FILE *fp_image;

  printf("Arrived in react_coords with %d atoms\n", num_atoms);

  printf("First atom in ref point %s %10.6f %10.6f %10.6f\n",
                  p_start_molecule->label, p_start_molecule->x,
                  p_start_molecule->y    , p_start_molecule->z);

  printf("First atom in end point %s %10.6f %10.6f %10.6f\n",
                  p_end_molecule->label, p_end_molecule->x,
                  p_end_molecule->y    , p_end_molecule->z);

/**********************************************************************/
/*** Check out distance between start and end point *******************/
/**********************************************************************/

d = mol_dist(p_start_molecule, p_end_molecule, num_atoms); 

dot = mol_dot(p_start_molecule, p_end_molecule, num_atoms);
mag1 = sqrt(mol_dot(p_start_molecule, p_start_molecule, num_atoms));
mag2 = sqrt(mol_dot(p_end_molecule, p_end_molecule, num_atoms));

theta = RAD_TO_DEG * acos(dot / (mag1*mag2));

printf("Distance between start and end structures = %10.6f A", d );
printf(" angle = %10.6f degrees\n", theta);

/**********************************************************************/
/*** Read in images in turn and analyse *******************************/
/**********************************************************************/

  im_have_group = FALSE;
  for (iloop=0; iloop < num_images; iloop++)
    {
       printf("image %d from file %s\n", iloop+1, p_image_files->name);

/****************************/

       open_file( &fp_image, p_image_files->name, "r");

       been = FALSE;

/**** No need to worry about coord flags here because we are only ***/
/**** working out the position of this image, not converting to   ***/
/**** POSCAR files. Hence use of dummy variable and FALSE         ***/
/**** Dave Willock Sept 07.                                       ***/

       good_read= read_car( fp_image, &header_line[0], &title_line[0],
                            &image[0], &date_line[0], &im_pbc, &im_num_atoms,
                            &im_num_of_mols, &im_num_mol_members[0], &im_mol_number[0],
                            &im_abc[0], &been, &dummy_group,
                            im_have_group, &im_num_in_group, FALSE, &dummy_flags );

       fclose(fp_image);
       
       if (good_read != 0)
         {
            printf("ERROR: Problem reading this image's car file\n");
            exit(0);
         }

       tot_d;
       p_start_atom= p_start_molecule;
       for (iatom=0; iatom < num_atoms; iatom++)
         {
/***********************************************************/
/** move end point atoms to min_image with start point *****/
/***********************************************************/
           if (pbc)
             {
               dx = image[iatom].x - p_start_atom->x;
               dy = image[iatom].y - p_start_atom->y;
               dz = image[iatom].z - p_start_atom->z;

               min_image( &dx, &dy, &dz, 
                               p_recip_latt_vec, p_latt_vec);

               image[iatom].x = p_start_atom->x + dx;
               image[iatom].y = p_start_atom->y + dy;
               image[iatom].z = p_start_atom->z + dz;
             }

       d = sqrt(atom_separation_squared(p_start_atom, &image[iatom], FALSE,  
                                        p_recip_latt_vec, p_latt_vec));

       printf("Atoms %d : %s (Start) and %s (image) are %10.6f apart\n",
                          iatom, p_start_atom->label, image[iatom].label, d);

       p_start_atom++;
     }

      d = mol_dist(p_start_molecule, &image[0], num_atoms); 

      dot = mol_dot(p_start_molecule, &image[0], num_atoms);
      mag1 = sqrt(mol_dot(p_start_molecule, p_start_molecule, num_atoms));
      mag2 = sqrt(mol_dot( &image[0], &image[0], num_atoms));

      theta = RAD_TO_DEG * acos(dot / (mag1*mag2));

      printf("Distance between start and image %d structures = %10.6f A",  iloop+1,d );
      printf(" angle = %10.6f degrees\n", theta);

/**********************************************************************/
/** Get closest point to image formed by linear combination of ********/
/** start and end                                              ********/
/** in this part of the code:                                  ********/
/** a is used for start point                                  ********/
/** b is used for end point                                    ********/
/** p is used for image                                        ********/
/** so that:                                                   ********/
/** ab is the dot product of start and end etc.                ********/
/**********************************************************************/

      ab = mol_dot(p_start_molecule, p_end_molecule, num_atoms);
      ap = mol_dot(p_start_molecule, &image[0], num_atoms);
      bp = mol_dot(p_end_molecule, &image[0], num_atoms);
      a2 = mol_dot(p_start_molecule, p_start_molecule, num_atoms);
      b2 = mol_dot(p_end_molecule, p_end_molecule, num_atoms);

      bot_end   = ab * ab - a2 * b2; 

      if ((a2        > 0.001 || a2        < -0.001) 
       && (bot_end   > 0.001 || bot_end   < -0.001))
        {
          amount_end   = ( ab * ap - a2 * bp ) / bot_end;
          amount_start = ( ap - amount_end * ab ) / a2;

          p_end_atom= p_end_molecule;
          p_start_atom= p_start_molecule;
          for (iatom=0; iatom < num_atoms; iatom++)
            {
               linear_comb_mol[iatom]= *p_start_atom; 

               linear_comb_mol[iatom].x = amount_start * p_start_atom->x + amount_end * p_end_atom->x;
               linear_comb_mol[iatom].y = amount_start * p_start_atom->y + amount_end * p_end_atom->y;
               linear_comb_mol[iatom].z = amount_start * p_start_atom->z + amount_end * p_end_atom->z;

               p_end_atom++;
               p_start_atom++;
            }

          printf("Nearest linear combination: %10.6f start, %10.6f end  ",
                  amount_start, amount_end);
          printf("Distance between image and linear combination = %10.6f Angstrom\n",
                               mol_dist(&image[0], &linear_comb_mol[0], num_atoms)); 
        }
      else
        {
          printf("ERROR: This image leads to a divide by zero\n");
          printf("       bot_start= %10.6f, bot_end = %10.6f\n", bot_start, bot_end);
          exit(0);
        }

      p_image_files++;
    }
  
  exit(0);
  return;
}
