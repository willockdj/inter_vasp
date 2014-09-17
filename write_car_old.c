/*********************************************************************/
/*** write_car will now just write out the co-ords it is sent   ******/
/*** without alterations. Any messing with the structure should ******/
/*** be done elsewhere! Dave Willock November 2005              ******/
/*********************************************************************/

#include <stdio.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"

/* ------Prototype-list---------------------------------------- */

void write_atom_data(FILE *fp, atom *p_atom, double scale_factor);

void put_string (FILE *fp, int *p_ichar, int length);

void min_image( double *x, double *y, double *z, double *p_recip_latt_vec, double *p_latt_vec);

void cart_to_fract( double cart_x,  double cart_y,  double cart_z, 
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z, 
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

/* ------------------------------------------------------------ */

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, double *p_recip_latt_vec )

{
   int iloop, mol_current, this_atom;
   int iavec, ibvec, icvec;
   double ta[3],tb[3],t[3], x,y,z;
   double fract_a, fract_b, fract_c;
   char *p_this_char;

   atom current_atom;
   atom *p_atom;

/*** Use the first element of header_line equal -1 as an indicator that */
/*** the header line has not been read                                  */
   if (start_frame)
     {
        if (*p_header_line != -1)
          {
             put_string( fp, p_header_line,100);
          }
        else
          {
             fprintf(fp, "!BIOSYM archive 3\n");
          }
     }

/* check if periodic boundaries were set */

   if (pbc) 
    {
      if (start_frame) fprintf(fp, "PBC=ON\n");

/*
      if (*p_title_line != -1)
        {
          put_string(fp, p_title_line,100);
          fprintf(fp, "\n");
        }
      else
        {
          fprintf(fp, "%s", p_c_title_line);
        }
*/

     fprintf(fp, "Inter_vasp generated file\n");

/*** Use the first element of date_line equal -1 as an indicator that */
/*** the header line has not been read                                  */

   if (*p_date_line != -1)
     {
      put_string(fp, p_date_line,100);
     }
   else
     {
      fprintf(fp, "!DATE Mon Oct 12 12:17:26 1998\n");
     }

      fprintf(fp,"PBC");

      for (iloop=0; iloop < 6; iloop++)
       {
         if ( iloop <= 2 )
           {
             fprintf(fp,"%10.4f",*(p_abc+iloop) * scale_factor * *(p_super+iloop));
           }
         else
           {
             fprintf(fp,"%10.4f",*(p_abc+iloop));
           }
       }
      fprintf(fp," (P1)\n");
    }
   else
    {
      if (start_frame) fprintf(fp,"PBC=OFF\n");
      put_string(fp, p_title_line,100);
      put_string(fp, p_date_line,100);
    }

   mol_current=0;


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

   for (iavec=0; iavec<= (*p_super)-1; iavec++)
     {
       ta[0]=iavec * *p_latt_vec;
       ta[1]=iavec * *(p_latt_vec+1);
       ta[2]=iavec * *(p_latt_vec+2);

       for (ibvec=0; ibvec<=*(p_super+1)-1; ibvec++)
         {
           tb[0]=ibvec * *(p_latt_vec+3);
           tb[1]=ibvec * *(p_latt_vec+4);
           tb[2]=ibvec * *(p_latt_vec+5);

           for (icvec=0; icvec<=*(p_super+2)-1; icvec++)
             {
               t[0]= ta[0] + tb[0] + icvec * *(p_latt_vec+6);
               t[1]= ta[1] + tb[1] + icvec * *(p_latt_vec+7);
               t[2]= ta[2] + tb[2] + icvec * *(p_latt_vec+8);

               p_atom=p_molecule;
               for (this_atom=0; this_atom < num_atoms; this_atom++)
                 {

                    if (*(p_mol_number+this_atom) != mol_current)
                      {
                         mol_current++;
                         fprintf(fp, "end\n");
                      }

                    current_atom = *p_atom;
                    current_atom.x += t[0];
                    current_atom.y += t[1];
                    current_atom.z += t[2];

                    write_atom_data(fp, &current_atom, scale_factor );
                    p_atom++;
                }
           }
       }
   }

fprintf(fp, "end\nend\n");

return;
}

