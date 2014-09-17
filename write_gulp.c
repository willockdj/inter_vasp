
#include <stdio.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

/* ------Prototype-list---------------------------------------- */

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void write_atom_data_gulp(FILE *fp, atom *p_atom, int is_core);

void put_string(FILE *fp, int *p_ichar, int length);

void min_image( double *x, double *y, double *z, double *p_recip_latt_vec, double *p_latt_vec);

/* ------------------------------------------------------------ */

void write_gulp(FILE *fp, int *p_title_line, atom *p_molecule, int num_atoms,
                atom *p_shells, int num_shells, int fract_or_cart, double *p_abc,
                int space_group, double *p_recip_latt_vec, double *p_latt_vec,
                charge_list *p_spec_charges, int num_species, int *p_super, 
                int need_shells, labels *p_shell_species, int num_shell_species)
  {
   int iloop, mol_current, this_atom;
   int ishell, put_shell, limits[3];
   int iavec, ibvec, icvec;
   double ta[3], tb[3], t[3];
   double x,y,z, fract_a, fract_b, fract_c;
   int is_core;

   atom *p_atom;
   atom current_atom;

   labels *p_this_lab;

   printf("Trying to write gulp file super cell %d %d %d\n", 
                             *p_super, *(p_super+1), *(p_super+2));

   fprintf(fp,"REPLACE WITH JOB OPTIONS\ntitle \n inter_vasp generated file ");
   fprintf(fp,"Super cell of dimensions %d %d %d based on original master file",
                             *p_super, *(p_super+1), *(p_super+2));
   fprintf(fp,"\nend\n");

   fprintf(fp,"cell\n%10.6f  %10.6f  %10.6f  %10.6f  %10.6f  %10.6f\n",
           *p_abc * *p_super, 
           *(p_abc+1) * *(p_super+1), 
           *(p_abc+2) * *(p_super+2), 
           *(p_abc+3), *(p_abc+4), *(p_abc+5));

   if (fract_or_cart == GULP_FRACT)
     {
       fprintf(fp,"fractional %d\n", space_group);
     }
   else
     {
       fprintf(fp,"cartesian\n");
     }

   limits[0]= *p_super;
   limits[1]= *(p_super+1);
   limits[2]= *(p_super+2);

   limits[0]--;
   limits[1]--;
   limits[2]--;

/*************************************************************/
/*** Move all atoms to their minimum image co-ordinates ******/
/*** with respect to atom 1                             ******/
/*************************************************************/

   p_atom=p_molecule;

   min_image( &(p_atom->x), &(p_atom->y), &(p_atom->z), p_recip_latt_vec, p_latt_vec);

   cart_to_fract( p_atom->x,  p_atom->y,  p_atom->z,
                  &fract_a, &fract_b, &fract_c,
                  p_recip_latt_vec );

   printf("fract co-ords: %10.6f %10.6f %10.6f\n", fract_a, fract_b, fract_c);

/***********************************************/
/*** Bring atoms back into unit cell ***********/
/***********************************************/
   if (fract_a < 0 && fract_a > -0.01)
     {
       p_atom->x += 0.01* *p_latt_vec;
       p_atom->y += 0.01* *(p_latt_vec+1);
       p_atom->z += 0.01* *(p_latt_vec+2);
     }
   else if (fract_a < 0)
     {
       p_atom->x += *p_latt_vec;
       p_atom->y += *(p_latt_vec+1);
       p_atom->z += *(p_latt_vec+2);
     }
   if (fract_b < 0 && fract_b > -0.01)
     {
       p_atom->x += 0.01* *(p_latt_vec+3);
       p_atom->y += 0.01* *(p_latt_vec+4);
       p_atom->z += 0.01* *(p_latt_vec+5);
     }
   else if (fract_b < 0)
     {
       p_atom->x +=  *(p_latt_vec+3);
       p_atom->y +=  *(p_latt_vec+4);
       p_atom->z +=  *(p_latt_vec+5);
     }
   if (fract_c < 0 && fract_c > -0.01)
     {
       p_atom->x += 0.01* *(p_latt_vec+6);
       p_atom->y += 0.01* *(p_latt_vec+7);
       p_atom->z += 0.01* *(p_latt_vec+8);
     }
   else if (fract_c < 0)
     {
       p_atom->x +=  *(p_latt_vec+6);
       p_atom->y +=  *(p_latt_vec+7);
       p_atom->z +=  *(p_latt_vec+8);
     }

   cart_to_fract( p_atom->x,  p_atom->y,  p_atom->z,
                  &fract_a, &fract_b, &fract_c,
                  p_recip_latt_vec );

   printf("fract co-ords: %10.6f %10.6f %10.6f\n\n", fract_a, fract_b, fract_c);

   for (iloop = 1; iloop < num_atoms; iloop++)
     {
       p_atom++;

       x= p_atom->x - p_molecule->x;
       y= p_atom->y - p_molecule->y;
       z= p_atom->z - p_molecule->z;

       min_image( &x, &y, &z, p_recip_latt_vec, p_latt_vec);

       p_atom->x = p_molecule->x + x;
       p_atom->y = p_molecule->y + y;
       p_atom->z = p_molecule->z + z;

   cart_to_fract( p_atom->x,  p_atom->y,  p_atom->z,
                  &fract_a, &fract_b, &fract_c,
                  p_recip_latt_vec );

   printf("fract co-ords: %10.6f %10.6f %10.6f\n", fract_a, fract_b, fract_c);
/***********************************************/
/*** Bring atoms back into unit cell ***********/
/***********************************************/
   if (fract_a < 0 && fract_a > -0.01)
     {
       p_atom->x += 0.01* *p_latt_vec;
       p_atom->y += 0.01* *(p_latt_vec+1);
       p_atom->z += 0.01* *(p_latt_vec+2);
     }
   else if (fract_a < 0)
     {
       p_atom->x += *p_latt_vec;
       p_atom->y += *(p_latt_vec+1);
       p_atom->z += *(p_latt_vec+2);
     }
   if (fract_b < 0 && fract_b > -0.01)
     {
       p_atom->x += 0.01* *(p_latt_vec+3);
       p_atom->y += 0.01* *(p_latt_vec+4);
       p_atom->z += 0.01* *(p_latt_vec+5);
     }
   else if (fract_b < 0)
     {
       p_atom->x +=  *(p_latt_vec+3);
       p_atom->y +=  *(p_latt_vec+4);
       p_atom->z +=  *(p_latt_vec+5);
     }
   if (fract_c < 0 && fract_c > -0.01)
     {
       p_atom->x += 0.01* *(p_latt_vec+6);
       p_atom->y += 0.01* *(p_latt_vec+7);
       p_atom->z += 0.01* *(p_latt_vec+8);
     }
   else if (fract_c < 0)
     {
       p_atom->x +=  *(p_latt_vec+6);
       p_atom->y +=  *(p_latt_vec+7);
       p_atom->z +=  *(p_latt_vec+8);
     }

   cart_to_fract( p_atom->x,  p_atom->y,  p_atom->z,
                  &fract_a, &fract_b, &fract_c,
                  p_recip_latt_vec );

   printf("fract co-ords: %10.6f %10.6f %10.6f\n\n", fract_a, fract_b, fract_c);

     }

   
   for (iavec=0; iavec <= limits[0]; iavec++)
     {
       ta[0]=iavec * *p_latt_vec;
       ta[1]=iavec * *(p_latt_vec+1);
       ta[2]=iavec * *(p_latt_vec+2);

       for (ibvec=0; ibvec <= limits[1]; ibvec++)
         {
           tb[0]=ibvec * *(p_latt_vec+3);
           tb[1]=ibvec * *(p_latt_vec+4);
           tb[2]=ibvec * *(p_latt_vec+5);

           for (icvec=0; icvec <= limits[2]; icvec++)
             {
               t[0]= ta[0] + tb[0] + icvec * *(p_latt_vec+6);
               t[1]= ta[1] + tb[1] + icvec * *(p_latt_vec+7);
               t[2]= ta[2] + tb[2] + icvec * *(p_latt_vec+8);

                p_atom= p_molecule;
                is_core = TRUE;

                printf("Writing atoms for translation %10.6f %10.6f %10.6f\n", t[0], t[1], t[2]);
                for (iloop = 0; iloop < num_atoms; iloop++)
                  {
                    current_atom = *p_atom;
                    current_atom.x += t[0];
                    current_atom.y += t[1];
                    current_atom.z += t[2];

                    if (fract_or_cart == GULP_FRACT)
                       {
                          x= current_atom.x;
                          y= current_atom.y;
                          z= current_atom.z;

                          cart_to_fract( x,  y,  z,  
                                         &fract_a, &fract_b, &fract_c, 
                                         p_recip_latt_vec );

                          p_atom->x= fract_a;
                          p_atom->y= fract_b;
                          p_atom->z= fract_c;

                          write_atom_data_gulp(fp, &current_atom, is_core);
 
                          put_shell=FALSE;
                          p_this_lab = p_shell_species; 
                          for (ishell=0; ishell<=num_shell_species; ishell++)
                            {
                             put_shell = put_shell 
                                         || (need_shells && strcmp(current_atom.label, &(p_this_lab->label[0])));
                             p_this_lab++;
                            }

                          if (put_shell)
                                write_atom_data_gulp(fp, &current_atom, FALSE);
                              
                       }
                    else
                       {
                          write_atom_data_gulp(fp, &current_atom, is_core);
 
                          put_shell=FALSE;
                          p_this_lab = p_shell_species; 
                          for (ishell=0; ishell<=num_shell_species; ishell++)
                            {
                             put_shell = put_shell 
                                         || (need_shells && strcmp(current_atom.label, &(p_this_lab->label[0])));
                             p_this_lab++;
                            }

                          if (put_shell)
                                write_atom_data_gulp(fp, &current_atom, FALSE);
                       }
     
                     p_atom++;
                  }
              }
          }
     }

/*******************************************************************************/
/*** Output species info *******************************************************/
/*******************************************************************************/

/*
   fprintf(fp,"species %d\n", num_species+1);
   for (iloop = 0; iloop <= num_species; iloop++)
     {
       if (p_spec_charges->is_core)
         {
           fprintf(fp,"%s core %10.6f\n",p_spec_charges->label, p_spec_charges->part_chge);

           if (need_shells && strcmp(p_spec_charges->label,"O") == 0)
                fprintf(fp,"%s shel %10.6f\n",
                   p_spec_charges->label, p_spec_charges->part_chge);
         }
       p_spec_charges++;
     }
*/
return;
}
