/*********************************************************************/
/*** write_car will now just write out the co-ords it is sent   ******/
/*** without alterations. Any messing with the structure should ******/
/*** be done elsewhere! Dave Willock November 2005              ******/
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

/* ------Prototype-list---------------------------------------- */

void write_atom_data(FILE *fp, atom *p_atom, double scale_factor, coord_flags *p_fix_flags);

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
		int *p_date_line, atom *p_molecule, int *p_mol_number, int use_mols,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, 
                double *p_recip_latt_vec, coord_flags *p_fix_flags, 
                double *p_magmom, int num_magmom )

{
   int iloop, mol_current, this_atom;
   int iavec, ibvec, icvec;
   int have_magmom;
   double ta[3],tb[3],t[3], x,y,z;
   double fract_a, fract_b, fract_c;
   char *p_this_char;

   atom current_atom;
   atom *p_atom;

   coord_flags *p_this_fix;

/*** Use the first element of header_line equal -1 as an indicator that */
/*** the header line has not been read                                  */

   if (start_frame)
     {
        printf("This is the first frame in write_car\n");
        if (*p_header_line != -1)
          {
             put_string( fp, p_header_line,100);
          }
        else
          {
             fprintf(fp, "!BIOSYM archive 3\n");
          }
     }

/* Check if magmom has been read */
    have_magmom=FALSE;
    if (num_magmom > -1) 
      {
         printf("write_car has magmom data supplied\n");
         have_magmom=TRUE;
         if (num_magmom+1 != num_atoms)
           {
             printf("ERROR: The number of MAGMOM entries (%d) does not match the number of atoms (%d)?\n",
                         num_magmom+1, num_atoms);
             exit(0);
           }
      }

/* check if periodic boundaries were set */

//   printf("DEBUG>> dealing with pbc\n");
   if (pbc) 
    {
      if (start_frame) fprintf(fp, "PBC=ON\n");
    }
   else
    {
      if (start_frame) fprintf(fp,"PBC=OFF\n");
    }


  if (*p_title_line != -1)
    {
      put_string(fp, p_title_line,100);
//          fprintf(fp, "\n");
    }
  else
    {
      fprintf(fp, "%s\n", p_c_title_line);
    }


//     fprintf(fp, "inter_vasp generated file\n");

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

   if (pbc)
     {
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

        mol_current=0;
//        printf("DEBUG>> writing..%d atoms\n", num_atoms);
//        printf("DEBUG>> super %d %d %d.c\n", *p_super, *(p_super+1), *(p_super+2));

        p_atom=p_molecule;
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
                  p_this_fix= p_fix_flags;
                  for (this_atom=0; this_atom < num_atoms; this_atom++)
                    {

                       if (*(p_mol_number+this_atom) != mol_current && use_mols)
                         {
                            mol_current++;
                            fprintf(fp, "end\n");
                         }

                       current_atom = *p_atom;
                       current_atom.x += t[0];
                       current_atom.y += t[1];
                       current_atom.z += t[2];

                       if (have_magmom)
                         {
                           current_atom.part_chge = *(p_magmom+this_atom);
                         }

//                     printf("Writing atom %d scale_factor: %10.6f\n", this_atom, scale_factor);
//                     printf("writing data for atom %s at %10.6f %10.6f %10.6f\n",
//                                        current_atom.label, current_atom.x,current_atom.y,current_atom.z);

                       write_atom_data(fp, &current_atom, scale_factor, p_this_fix );
                       p_this_fix++;
                       p_atom++;
                   }
              }
          }
      }
  }
else
   {
       mol_current=0;
       p_atom=p_molecule;
       p_this_fix= p_fix_flags;
       for (this_atom=0; this_atom < num_atoms; this_atom++)
         {

            if (*(p_mol_number+this_atom) != mol_current && use_mols)
              {
                 mol_current++;
                 fprintf(fp, "end\n");
              }

            write_atom_data(fp, p_atom, scale_factor, p_this_fix );
            p_this_fix++;
            p_atom++;
        }

   }

fprintf(fp,    "end\nend\n");

return;
}

