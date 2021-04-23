/***************************************************************************/
/*** Read in a VASP XDATCAR file for MD trajectories ***********************/
/*** Dave Willock May 06 ***************************************************/
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

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

int read_line(FILE *fp, int *p_ichar);

int read_atom_data_vasp( FILE *fp, atom *p_atom, int *p_mol_number,
                          coord_flags *p_fix_flags );

double get_double(FILE *input_fp, int skip_lines, int *p_error);

int get_integer( FILE *input_fp, int skip_lines, int *p_error );

void put_string (FILE *fp, int *p_ichar, int length);

void find_line(FILE *fp, char *p_key, char *p_key2, int sep, int *p_found, int max_lines);

void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                          double *p_abc );

/*---------------------------------------------------------------------------*/

/* read in a xdatcar vasp formatted file */

int read_xdatcar( FILE *fp, int *p_num_frames, atom *p_molecule, int num_atoms,
                  coord_flags *p_fix_flags, double *p_latt_vec, double *p_recip_latt_vec,
                  double *p_abc, int *p_is_fract, int *p_is_cart, labels *p_atom_names, 
                  int *p_atom_num, int *p_num_labels, int just_count )
{
  int ichar[LINESIZ];
  int i,iatom,idummy, ntype;
  int skip, lower_case=FALSE;
  int error, found, num_of_chars;
  int *p_this_atom_num;

  char *p_key;
  char *tok;

  labels *p_this_atom_name;

/* assume xdatcar file is just a read after atoms have been sorted out */

   printf("Trying to read XDATCAR file format assume direct co-ordinates\n");
   *p_is_fract=TRUE; *p_is_cart=FALSE;
   p_key = "configuration=";

   if (*p_num_frames < 0)
     {
/* Read lattice vectors from top of file ****/

/* Skip header and index lines           ****/

        num_of_chars = read_line( fp, &ichar[0]);
        num_of_chars = read_line( fp, &ichar[0]);

        skip= TRUE;
        *p_latt_vec        = get_double(fp, skip, &error);
        skip= FALSE;
        *(p_latt_vec+1)    = get_double(fp, skip, &error);
        *(p_latt_vec+2)    = get_double(fp, skip, &error);

        skip= TRUE;
        *(p_latt_vec+3)    = get_double(fp, skip, &error);
        skip= FALSE;
        *(p_latt_vec+4)    = get_double(fp, skip, &error);
        *(p_latt_vec+5)    = get_double(fp, skip, &error);

        skip= TRUE;
        *(p_latt_vec+6)    = get_double(fp, skip, &error);
        skip= FALSE;
        *(p_latt_vec+7)    = get_double(fp, skip, &error);
        *(p_latt_vec+8)    = get_double(fp, skip, &error);

        latt_vecs_from_cart( p_latt_vec, p_recip_latt_vec, p_abc );

        skip=TRUE;
/* Read atom labels */
        i=0; p_this_atom_name=p_atom_names;
        while ( (tok=tok_get(fp, skip, FALSE)) )
          {
            i++; 
            if ( i > MAXTYPES ) 
              {
                 printf("ERROR: Too many names in XDATCAR\n");
                 exit(0);
              }
 
            strcpy( &(p_this_atom_name->label[0]),tok);
            p_this_atom_name++;
            skip=FALSE;
          }

/** Read atom numbers **/

        skip=TRUE; *p_num_labels=-1;
        p_this_atom_num= p_atom_num;
        while ( tok=tok_get(fp, skip, FALSE) )
          {
            *p_this_atom_num= atoi(tok);
            p_this_atom_num++;
            ++*p_num_labels;
            skip=FALSE;
          }

        p_this_atom_name=p_atom_names;
        p_this_atom_num= p_atom_num;
        for (i=0; i<=*p_num_labels; i++) 
           {
             printf("%s.......%d\n", p_this_atom_name->label , *p_this_atom_num);
             p_this_atom_name++; p_this_atom_num++;
           }

/* go through file looking for configuration start markers */

        found = TRUE;
        while (found)
          {
             find_line( fp, "none", p_key, 0, &found, -1);
             if (found) ++*p_num_frames; 
          }
        printf("Found %d frames....\n", *p_num_frames);


        return 0;

    }
  else
    {
 /**** Must already have number of frames and am now just reading one in ****/
      skip = TRUE;
      printf("Back in read_xdatcar for another frame\n");
      printf("Expecting %d atoms\n", num_atoms);
 /**** First skip blank line ***/
/*      while (!(tok= tok_get(fp, skip, lower_case))); */
/*      printf("TOK: %s\n",tok);  */

      find_line( fp, "none", p_key, 0, &found, -1);

      if (found)
        {
           printf("Found header...\n");
           ntype=0; 
           p_this_atom_name=p_atom_names; p_this_atom_num= p_atom_num;

           for (iatom=0; iatom < num_atoms; iatom++)
             {
               printf("reading atom data...\n");
               read_atom_data_vasp( fp, p_molecule, &idummy, p_fix_flags );

               printf("Read: %10.6f  %10.6f  %10.6f\n", p_molecule->x,
                                                        p_molecule->y,
                                                        p_molecule->z);


/** Make up label and identity fields for car file **/
               strcpy(p_molecule->elem, p_this_atom_name);
               sprintf(p_molecule->label,"%s%d",p_this_atom_name,ntype+1);
               strcpy(p_molecule->pot, "?");
               strcpy(p_molecule->group, "VASP");

               ntype++;
               if (ntype == *p_this_atom_num)
                 {
                   ntype = 0; p_this_atom_num++;
                   p_this_atom_name++;
                 }

               p_molecule++;
            }
        }
    }

   printf("Frame read back to main.....\n") ;

  return 0;
}
