/*********************************************************************/
/** Routine sets up the required group member lists for group ********/
/** interpolation options.                                    ********/
/** Broke out of main code, April 2019, Dave Willock          ********/
/*********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "maxima.h"
#include "constants.h"
#include "structures.h"
#include "reader.h"
#include "global_values.h"

void inter_atom_vector(atom *p_atom, atom *p_atom2, double *p_vec);

double vec_dot(double *p_A, double *p_B);

void unit_vector(double *p_vector);

double size_vector(double *p_vector);


void set_up_groups( atom *p_molecule, atom *p_end_molecule, int num_atoms,
                    group_lists *p_groups, int num_groups,
                    double *p_orig_bond1, double *p_orig_bond2,
                    double *p_delta_bond1, double *p_delta_bond2,
                    have_list *p_have, need_list *p_need)
 {
 int iloop, check, check1, check2;
 int grp_mem_ind, num_in_group, igrp;

 double origin[3], axis[3], dot;
 double vec[3], vec1[3], vec2[3];

 atom *p_atom, *p_atom2, *p_end_atom, *p_end_atom2;

 printf("In set_up_groups have %d groups:\n", num_groups);

 printf("GRP1 defined as containing %d atoms:\n", p_groups->num_grp1+1);

 if (num_groups > 0)
   {
     printf("GRP2 defined as containing %d atoms:\n", p_groups->num_grp2+1);
   }

 check  = FALSE;
 check1 = FALSE;
 check2 = FALSE;

 if (p_need->angle)
  {
    printf("Need angle set\n");
    for (iloop=0; iloop<=p_groups->num_grp1; iloop++)
      {
        grp_mem_ind= p_groups->group1[iloop];
        p_atom= p_molecule + grp_mem_ind;
        printf("%d %d >>%s<<",iloop, grp_mem_ind, p_atom->label );

        if (strcmp(p_atom->label, p_groups->axis1_lab) == 0)
           {
             check1=TRUE;
             p_groups->axis1_ind[0] = grp_mem_ind;
             printf("   the first axis definition atom\n");
           }
        else if (strcmp(p_atom->label, p_groups->axis2_lab) == 0)
           {
             check2=TRUE;
             p_groups->axis2_ind[0] = grp_mem_ind;
             printf("   the second axis definition atom\n");
           }
        else
           {
             printf("\n");
           }
      }
    if (!check1 && !check2)
      {
        printf("ERROR: The axis definition atoms given in the input do not occur in GRUP\n");
        exit(0);
      }
    else if (!check1)
      {
        printf("ERROR: The first axis definition atom given in the input does not occur in GRUP\n");
        exit(0);
      }
    else if (!check2)
      {
        printf("ERROR: The second axis definition atom given in the input does not occur in GRUP\n");
        exit(0);
      }
  }
 else
  {
    for (iloop=0; iloop<=p_groups->num_grp1; iloop++)
      {

/** Identify the central atoms in the group lists ***/

        grp_mem_ind= p_groups->group1[iloop];
        p_atom= p_molecule + grp_mem_ind;
        printf("%d %d >>%s<<",iloop, grp_mem_ind, p_atom->label );

        if (strcmp(p_atom->label, p_groups->cnt_lab1) == 0)
/**** start from a zero index the space required is **/
           {
             check=TRUE;
             p_groups->centre[0] = grp_mem_ind;
             printf("   the central atom\n");
           }
        else
           {
             printf("\n");
           }
       }
     if (!check)
       {
         printf("ERROR: The central atom given in the input >>%s<< does not occur in GRUP or GRP1\n", p_groups->cnt_lab1);
         exit(0);
       }

/** Deal with second group if present **/
     if (num_groups > 0 && p_groups->num_grp2 > -1 )
       {
          printf("Have a second group with %d members\n", p_groups->num_grp2+1);
          for (iloop=0; iloop<=p_groups->num_grp2; iloop++)
            {

/** Identify the central atoms in the group lists ***/

              grp_mem_ind= p_groups->group2[iloop];
              p_atom= p_molecule + grp_mem_ind;
              printf("%d %d >>%s<<",iloop, grp_mem_ind, p_atom->label );

              if (strcmp(p_atom->label, p_groups->cnt_lab2) == 0)
/**** start from a zero index the space required is **/
                 {
                   check1=TRUE;
                   p_groups->centre[1] = grp_mem_ind;
                   printf("   the central atom for group 2\n");
                 }
              else
                 {
                   printf("\n");
                 }
            }
          if (!check1)
            {
              printf("ERROR: The central atom for group 2 given in the input >>%s<< does not occur in GRP2\n", p_groups->cnt_lab2);
              exit(0);
            }
       }
   }

/**** Now work out vectors within the group for the interpolation ****/

if (p_need->angle)
  {
/************************************************************/
/*** Use rotation of GRUP atoms around a vector defined *****/
/*** by the two atom labels supplied.                   *****/
/************************************************************/

     p_atom =p_molecule+p_groups->axis1_ind[0];
     p_atom2=p_molecule+p_groups->axis2_ind[0];
     printf("Interpolating by angle about defined axis\n");
     printf("Axis defined by atoms %d %s and %d %s\n",
                  p_groups->axis1_ind[0], p_atom->label, 
                  p_groups->axis2_ind[0], p_atom2->label );

/************************************************************/
/* Define axis vector ***************************************/
/************************************************************/

     origin[0] = p_atom->x;
     origin[1] = p_atom->y;
     origin[2] = p_atom->z;

     inter_atom_vector(p_atom, p_atom2, &axis[0]);

     unit_vector(&axis[0]); 

/*************************************************************/
/* Work out vector from axis to each of the grouped atoms ****/
/* Check first which group has this type of iteration     ****/
/*************************************************************/
    if (p_groups->group_type1 == ANGLE_TYPE) num_in_group = p_groups->num_grp1;
    if (p_groups->group_type2 == ANGLE_TYPE) num_in_group = p_groups->num_grp2;
                     
    for (igrp=0; igrp<=num_in_group; igrp++)
      {
         if (p_groups->group_type1 == ANGLE_TYPE) grp_mem_ind = p_groups->group1[igrp]; 
         if (p_groups->group_type2 == ANGLE_TYPE) grp_mem_ind = p_groups->group2[igrp]; 

         if (grp_mem_ind != p_groups->axis1_ind[0] && grp_mem_ind != p_groups->axis2_ind[0])   
            {
/*************************************************************/
/** vector from atom to one of axis defining atoms which *****/
/** must be on the rotation axis.                        *****/
/*************************************************************/
               inter_atom_vector(p_atom, p_molecule+grp_mem_ind, &vec[0]); 

               dot = vec_dot(&vec[0], &axis[0]);

/*************************************************************/
/** Take component which is along axis and remove to get *****/
/** vector to atom perpendicular to the axis.            *****/
/*************************************************************/
               vec1[0]= vec[0] - dot * axis[0];
               vec1[1]= vec[1] - dot * axis[1];
               vec1[2]= vec[2] - dot * axis[2];

               unit_vector(&vec1[0]); 

/*************************************************************/
/** Repeat for the end_point molecule                    *****/
/*************************************************************/
               inter_atom_vector(p_end_molecule+p_groups->axis1_ind[0], p_end_molecule+grp_mem_ind, &vec[0]); 

               dot = vec_dot(&vec[0], &axis[0]);

               vec2[0]= vec[0] - dot * axis[0];
               vec2[1]= vec[1] - dot * axis[1];
               vec2[2]= vec[2] - dot * axis[2];

               unit_vector(&vec2[0]); 

/*************************************************************/
/*** Total angle to rotate by is then obtained from the ******/
/*** vec1, vec2 dot product.                            ******/
/*************************************************************/

               dot = vec_dot(&vec1[0], &vec2[0]);

               (p_molecule+grp_mem_ind)->theta = -acos(dot);

               printf("Atom %d %s will be rotated through %10.6f degrees during interpolation\n",
                                     grp_mem_ind, (p_molecule+grp_mem_ind)->label, 
                                     RAD_TO_DEG*(p_molecule+grp_mem_ind)->theta);

            }
       }
  }
else
  {
    printf("\n");
/********************************************************************/ 
/*** For the case of grouped atoms work out changes of bond length **/
/*** May have two groups, treat seperately.                        **/
/*** Put group centre to p_atom and p_end_atom.                    **/
/********************************************************************/ 
    if (p_groups->group_type1 == CENTRE_TYPE)
      {
         p_atom= p_molecule+p_groups->centre[0];
         p_end_atom= p_end_molecule+p_groups->centre[0];

         for (igrp=0; igrp<=p_groups->num_grp1; igrp++)
           {
            grp_mem_ind = p_groups->group1[igrp];  

/*** p_atom2 for the other group members ***/
            p_atom2=p_molecule+grp_mem_ind;
            p_end_atom2=p_end_molecule+grp_mem_ind;
            printf("PICKED up index %d as atom %s\n",grp_mem_ind, p_atom2->label);
            printf("PICKED up centre index %d \n",p_groups->centre[0]);

            if (grp_mem_ind != p_groups->centre[0])
              {
                 inter_atom_vector(p_atom2, p_atom, &vec[0]); 

                 inter_atom_vector(p_end_atom2, p_end_atom, &vec1[0]); 

                 *p_orig_bond1= size_vector(&vec[0]);

                 *p_delta_bond1= size_vector(&vec1[0])-*p_orig_bond1;

                 printf("centre to atom %d vector is %10.6f %10.6f %10.6f\n",
                               grp_mem_ind, vec[0],vec[1],vec[2]);

             printf("This atom is %s the centre is %s\n", p_atom2->label, 
                                                           p_atom->label);
             printf("Original bond length for atom  %d %s - %s is %10.6f\n", grp_mem_ind, 
                                                                             p_atom->label, 
                                                                             p_atom2->label, 
                                                                             *p_orig_bond1);

             printf("Final    bond length for atom  %d %s - %s is %10.6f\n", grp_mem_ind, 
                                                                             p_end_atom->label, 
                                                                             p_end_atom2->label, 
                                                                             size_vector(&vec1[0]));

             printf("Change in bond length for atom %d is %10.6f\n", grp_mem_ind, *p_delta_bond1);

             }
           else
             {
               *p_orig_bond1=0.0; *p_delta_bond1=0.0; 
               printf("This is the group centre, no need to test!\n");
             }
            p_orig_bond1++; p_delta_bond1++;
          }
      }
/**** If have second group sort that out now ****/
    if (num_groups > 0 && p_groups->group_type2 == CENTRE_TYPE )
      {
         printf("\nHave two groups to work with for the second:\n");
         p_atom= p_molecule+p_groups->centre[1];
         p_end_atom= p_end_molecule+p_groups->centre[1];
         printf("\nSecond group centre index %d is atom %s:\n", p_groups->centre[1], p_atom->label);

         for (igrp=0; igrp<=p_groups->num_grp2; igrp++)
           {
             grp_mem_ind = p_groups->group2[igrp];  
             p_atom2=p_molecule+grp_mem_ind;
             p_end_atom2=p_end_molecule+grp_mem_ind;

             if (grp_mem_ind != p_groups->centre[1])
               {
                  inter_atom_vector(p_atom2, p_atom, &vec[0]); 
                  inter_atom_vector(p_end_atom2, p_end_atom, &vec1[0]); 

                  *p_orig_bond2 = size_vector(&vec[0]);

                  *p_delta_bond2= size_vector(&vec1[0]) - *p_orig_bond2;

                  printf("centre to atom %d vector is %10.6f %10.6f %10.6f\n",
                               grp_mem_ind,
                               vec[0],vec[1],vec[2]);

                  printf("This atom is %s the centre is %s\n", p_atom2->label, 
                                                                             p_atom->label);
                  printf("Original bond length for atom  %d is %10.6f\n", grp_mem_ind, *p_orig_bond2);
                  printf("Change in bond length for atom %d is %10.6f\n", grp_mem_ind, *p_delta_bond2);
              }
           else
             {
               *p_orig_bond2=0.0; *p_delta_bond2=0.0; 
               printf("This is the group centre, no need to test!\n");
             }
            p_orig_bond2++; p_delta_bond2++;
          }
      }
    printf("\n\nData for grouped interpolation set up\n\n");
  }

return;
}
