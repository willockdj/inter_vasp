/********************************************************************/
/*** read_cif reads in a .cif file formatted structure **************/
/*** recognises group labels:                          **************/
/***                                                   **************/
/*** Begun 20th March 2020, Dave Willock               **************/
/********************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "global_values.h"
#include "maxima.h"
#include "debug.h"
#include "structures.h"

/* protype list for this routine */

int read_line(FILE *fp, int *p_ichar);

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int read_atom_data(int *p_ichar, int *p_num_atoms, int at_end, int num_of_chars,
                   int *p_num_of_mols, atom *p_atom, int *p_mol_number,
                   int just_count );

int get_int(int *p_ichar,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

int next_space( int *p_ichar, int start, int num_of_chars );

int next_none_space( int *p_ichar, int start, int num_of_chars );

void copy_int( int *p_ichar1, int *p_ichar2, 
                                  int min_position, int max_position );

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position );

void put_string(FILE *fp, int *p_ichar, int length);
/*---------------------------------------------------------------------------*/

/* read in a .car insight formatted file */

int read_cif( FILE *fp, int *p_title_line, 
              atom *p_molecule, int *p_pbc, int *p_num_atoms, 
              double *p_abc, int *p_been_before, int *p_set_labels,
              is_list *p_is, int just_count)
{
  int ichar[LINESIZ], idum[10];
  int iatom, place, iend;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int at_end, j, atom_lines, loop_line;
  int start,start2,start3,start4;
  int *p_this_index, num_lines;
  int isign, itsanum, ndigi, max_chars;
  int space_group, have_loop;
  int nxt_spc, nxt_non_spc;
  int got_label, got_symbol, symbol_entry, label_entry;

  char *p_key;

  atom *p_atom;

  atom_number *p_this_type;

  *p_pbc = TRUE;

/* just go through file looking for flags */

//_space_group_IT_number           62
//_symmetry_space_group_name_Hall  '-P 2c 2ab'
//_symmetry_space_group_name_H-M   'P b n m'
//_cell_angle_alpha                90
//_cell_angle_beta                 90
//_cell_angle_gamma                90
//_cell_length_a                   5.5202
//_cell_length_b                   5.5202
//_cell_length_c                   7.8067
//
have_loop = FALSE;
atom_lines=-1;
got_label=FALSE;
got_symbol=FALSE;
num_lines=-1;
*p_num_atoms=0;
p_atom= p_molecule;

if (just_count)
  {
    printf("....................just taking a look to count the atoms........\n");
  }
else
  {
    printf("....................here to get the atoms on board this time.....\n");
  }
  

loop_line = num_lines;

 while ( (num_of_chars = read_line(fp, &ichar[0])) != -10)
    {
       num_lines++;
       put_string( stdout, &ichar[0], num_of_chars); 
       printf("  line has %d characters\n", num_of_chars);
//
// Look for space group number
       p_key="_space_group_IT_number";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 22;
           space_group = get_int( &ichar[0], &j, &itsanum,
                                  &ndigi, num_of_chars, &isign);
           printf("Line with space group index: %d\n", space_group);

           if (space_group != 1)
             {
                printf("ERROR: space group is not P1 so cannot continue....\n");
                exit(0);
             }
         }
// Look for a vector              
       p_key="_cell_length_a";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 14;
           *p_abc = get_doub( &ichar[0], num_of_chars, &j, &itsanum );
           printf("Line with a vector: %10.6f\n", *p_abc);
         }
// Look for b vector              
       p_key="_cell_length_b";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 14;
           *(p_abc+1) = get_doub( &ichar[0], num_of_chars, &j, &itsanum );
           printf("Line with b vector: %10.6f\n", *(p_abc+1));
         }
// Look for c vector              
       p_key="_cell_length_c";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 14;
           *(p_abc+2) = get_doub( &ichar[0], num_of_chars, &j, &itsanum );
           printf("Line with c vector: %10.6f\n", *(p_abc+2));
         }
// Look for alpha                 
       p_key="_cell_angle_alpha";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 17;
           *(p_abc+3) = get_doub( &ichar[0], num_of_chars, &j, &itsanum );
           printf("Line with alpha angle: %10.6f\n", *(p_abc+3));
         }
// Look for beta                  
       p_key="_cell_angle_beta";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 15;
           *(p_abc+4) = get_doub( &ichar[0], num_of_chars, &j, &itsanum );
           printf("Line with beta angle: %10.6f\n", *(p_abc+4));
         }
// Look for gamma                 
       p_key="_cell_angle_gamma";
       if (locate_string( p_key, &ichar[0], num_of_chars ))
         {
           j= 16;
           *(p_abc+5) = get_doub( &ichar[0], num_of_chars, &j, &itsanum );
           printf("Line with gamma angle: %10.6f\n", *(p_abc+5));
         }
/* Read in the atomic data */
//loop_
//_atom_site_label
//_atom_site_fract_x
//_atom_site_fract_y
//_atom_site_fract_z
//_atom_site_U_iso_or_equiv
//Sr 0.00000 0.00000 0.25000 0.00608
//Ti 0.50000 0.00000 0.00000 0.00507
//O1 0.00000 0.50000 0.25000 0.01646
//O2 0.75000 0.25000 0.00000 0.02026
// 
// or
//    _atom_site_label
//    _atom_site_type_symbol
//    _atom_site_fract_x
//    _atom_site_fract_y
//    _atom_site_fract_z
//    _atom_site_occupancy
//

       if (have_loop)
         {
            printf("Have loop....\n");
            p_key="_atom_site";
            if (locate_string( p_key, &ichar[0], num_of_chars )) 
              {
                atom_lines++;     
                printf(".....atom_line now %d\n",atom_lines);

                p_key="_atom_site_label";
                if (locate_string( p_key, &ichar[0], num_of_chars ))
                  {
                    got_label=TRUE;
                    label_entry=atom_lines;
                    printf("label is entry %d\n", atom_lines);
                  }

                p_key="_atom_site_type_symbol";
                if (locate_string( p_key, &ichar[0], num_of_chars ))
                  {
                    got_symbol=TRUE;
                    symbol_entry=atom_lines;
                    printf("symbol is entry %d\n", atom_lines);
                  }

                p_key="_atom_site_Cartn";
                if (locate_string( p_key, &ichar[0], num_of_chars )) p_is->cart=TRUE;

                p_key="_atom_site_fract";
                if (locate_string( p_key, &ichar[0], num_of_chars )) p_is->fract=TRUE;

              }
            else if (atom_lines < 0) 
              {
                printf("Have_loop reset\n");
                have_loop = FALSE;
              }
            else
              {
                printf("Atom data line ....\n"); 
                ++*p_num_atoms;

                if (!just_count)
                  {
/* get label and symbol if present elem label */ 
                     printf("Actually filling arrays now....\n");

                     j= 0;
                     if (got_label && got_symbol)
                       {
                          printf("Looking for label and symbol (elem)\n");
                          nxt_non_spc=next_none_space(&ichar[0], j, num_of_chars);
                          nxt_spc= next_space( &ichar[0], nxt_non_spc+1, num_of_chars );
                          
                          printf("Next none space: %d, next_space: %d\n", nxt_non_spc, nxt_spc);
// next_non_spc and next_spc bracket the word
// Use to get label or symbol depending on order
                          if ( nxt_spc > nxt_non_spc && nxt_spc - nxt_non_spc < 5 )
                            {
                              copy_int( &ichar[0], &idum[0], nxt_non_spc, nxt_spc); 

                              if ( label_entry < symbol_entry)
                                {
                                   int_to_string(&idum[0],&(p_atom->label[0]), nxt_spc-nxt_non_spc);
                                }
                              else
                                {
                                   int_to_string(&idum[0],&(p_atom->elem[0]), nxt_spc-nxt_non_spc);
                                }
                            }
                          else
                            {
                               printf("ERROR reading atom label in cif file\n");
                               exit(0);
                            }
// Now get symbol  
                          nxt_non_spc=next_none_space(&ichar[0], nxt_spc, num_of_chars);
                          nxt_spc= next_space( &ichar[0], nxt_non_spc+1, num_of_chars );
                          
                          printf("Next none space: %d, next_space: %d\n", nxt_non_spc, nxt_spc);
// next_non_spc and next_spc bracket the word
// Use to get label or symbol depending on order
                          if ( nxt_spc > nxt_non_spc && nxt_spc - nxt_non_spc < 5 )
                            {
                              copy_int( &ichar[0], &idum[0], nxt_non_spc, nxt_spc); 

                              if ( label_entry < symbol_entry)
                                {
                                  int_to_string(&idum[0],&(p_atom->elem[0]), nxt_spc-nxt_non_spc);
                                }
                              else
                                {
                                   int_to_string(&idum[0],&(p_atom->label[0]), nxt_spc-nxt_non_spc);
                                }
                            }
                          else
                            {
                               printf("ERROR reading atom element in cif file\n");
                               exit(0);
                            }

                          printf("Read label >>%s<< and symbol >>%s<<\n", p_atom->label, p_atom->elem);

// No need for generate_neighbours to do labels
                          *p_set_labels = FALSE;
                       }
                     else if (got_label)
                       {
                          printf("Just looking for label\n");
                          nxt_non_spc=next_none_space(&ichar[0], j, num_of_chars);
                          nxt_spc= next_space( &ichar[0], nxt_non_spc+1, num_of_chars );
                          
                          printf("Next none space: %d, next_space: %d\n", nxt_non_spc, nxt_spc);
// next_non_spc and next_spc bracket the word
// Use to get label
                          if ( nxt_spc > nxt_non_spc && nxt_spc - nxt_non_spc < 5 )
                            {
                              copy_int( &ichar[0], &idum[0], nxt_non_spc, nxt_spc); 

                              int_to_string(&idum[0],&(p_atom->label[0]), nxt_spc-nxt_non_spc);
// Need to strip any numbers from label to get element
                              iend=0;
                              while (iend <= nxt_spc-nxt_non_spc && (idum[iend] < '0' || idum[iend] > '9' )) iend++;
                              int_to_string(&idum[0],&(p_atom->elem[0]), iend);
                            }
                          else
                            {
                               printf("ERROR reading atom label in cif file\n");
                               exit(0);
                            }


                          printf("Read label >>%s<< set symbol >>%s<<\n", p_atom->label, p_atom->elem);
// No need for generate_neighbours to do labels
                          *p_set_labels = FALSE;
                       }
                     else if (got_symbol)
                       {
                          printf("Just looking for symbol (elem)\n");
// Now get symbol  
                          nxt_non_spc=next_none_space(&ichar[0], j, num_of_chars);
                          nxt_spc= next_space( &ichar[0], nxt_non_spc+1, num_of_chars );
                          
                          printf("Next none space: %d, next_space: %d\n", nxt_non_spc, nxt_spc);
// next_non_spc and next_spc bracket the word
// Use to get label
                          if ( nxt_spc > nxt_non_spc && nxt_spc - nxt_non_spc < 5 )
                            {
                              copy_int( &ichar[0], &idum[0], nxt_non_spc, nxt_spc); 

                              int_to_string(&idum[0],&(p_atom->elem[0]), nxt_spc-nxt_non_spc);
                            }
                          else
                            {
                               printf("ERROR reading atom element in cif file\n");
                               exit(0);
                            }

                          printf("Read symbol >>%s<<\n", p_atom->elem);

// Ask generate_neighbours to do labels
                          *p_set_labels = TRUE;
                       }
                     else
                       {
                          printf("ERROR: Neither labels or symbols set in the cif file.");
                          exit(0);
                       }

           
//                  if (nxt_spc > j && nxt_spc-j < 4)
//                  copy_int( &ichar[0], &idum[0], j, nxt_spc); 
//                  int_to_string(&idum[0],&(p_atom->elem[0]), nxt_spc-j);

                    j=nxt_spc;
                    printf("j = nxt_spc = %d\n", j);
                    p_atom->x = get_doub(&ichar[0], num_of_chars, &j, &itsanum);
                    printf("read x= %10.6f, now j = %d\n", p_atom->x, j);
                    p_atom->y = get_doub(&ichar[0], num_of_chars, &j, &itsanum);
                    printf("read y= %10.6f, now j = %d\n", p_atom->y, j);
                    p_atom->z = get_doub(&ichar[0], num_of_chars, &j, &itsanum);
                    printf("read z= %10.6f, now j = %d\n", p_atom->z, j);
/* Make sensible choices for labels and groups */
                    p_atom->part_chge = 0.0;
                    strcpy(p_atom->group, "VASP");
                    strcpy(p_atom->pot, "?");

                    p_atom++;
                  }
              }
         }

       p_key="loop_";
       if (locate_string( p_key, &ichar[0], num_of_chars )) have_loop = TRUE;
    }
 
   printf("abc= %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n", 
               *p_abc, *(p_abc+1), *(p_abc+2), *(p_abc+3), *(p_abc+4), *(p_abc+5));
   printf("cif file contains %d atoms\n", *p_num_atoms);

   if (!just_count)
     {
        p_atom= p_molecule;
        for (iatom=0; iatom < *p_num_atoms; iatom++)
          {
             printf("%s %10.6f %10.6f %10.6f\n", p_atom->elem,
                                                 p_atom->x,
                                                 p_atom->y,
                                                 p_atom->z);
             p_atom++;
          }
     }

  return 0;
}
