/********************************************************************/
/*** read_car reads in a .car file formatted structure **************/
/*** recognises group labels:                          **************/
/*** GRUP : Atoms that are part of the interpolation   **************/
/***        group.                                     **************/
/*** GRP1 : GRP1 and GRP2 labels introduced to allow   **************/
/***        interpolations with a bond between two     **************/
/*** GRP2 : molecular groups being broken.             **************/
/*** FIXX : Atoms that should be fixed: if requested   **************/
/***                                                   **************/
/*** last updated April 30th 2009, Dave Willock & Adam T ************/
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

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

void put_string(FILE *fp, int *p_ichar, int length);
/*---------------------------------------------------------------------------*/

/* read in a .car insight formatted file */

int read_car( FILE *fp, int *p_header_line, int *p_title_line, 
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms, 
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number, 
              double *p_abc, int *p_been_before, group_lists *p_groups, 
              int have_grp, int *p_num_groups, int find_fixed, 
              coord_flags *p_fix_flags, int just_count)
{
  int ichar[LINESIZ];
  int place,itsanum;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int at_end;
  int start,start2,start3,start4;
  int *p_this_index;

  char *p_key;

  atom *p_atom;

/* deal with top of file on first pass only */


 if (!*p_been_before)
   {

/* read in first line and check for BIOSYM header */

     *p_been_before= 1;

     num_of_chars = read_line( fp, p_header_line);
     p_key= "BIOSYM";
     is_biosym= locate_string( p_key, p_header_line, num_of_chars);

     if (!is_biosym)
       {
          printf("The BIOSYM header is missing from this file check it really is a .car file \n");
          return 1;
       }
     else
       {
          printf("This is a BIOSYM file\n");
       }

/* check for periodic boundary conditions */

     num_of_chars = read_line( fp, &ichar[0]);
     p_key = "PBC=ON";
     *p_pbc= locate_string( p_key, &ichar[0], num_of_chars);

   }

/* Read in lines common to all frame entries */

/* get title and date lines */

     num_of_chars= read_line ( fp, p_title_line);

     if (num_of_chars < 0) *p_title_line = -1;

     num_of_chars= read_line ( fp, p_date_line); 


/* if pbc set get abc alpha beta gamma */

   if (*p_pbc == 1) 
     {
        num_of_chars= read_line(fp, &ichar[0]);

        place= -1;
        for (iloop=0; iloop <= 5; iloop++) 
          {
             *(p_abc+iloop)= get_doub(&ichar[0], num_of_chars, &place, &itsanum); 
             if (!itsanum)
	       {
                 printf("failure whilst trying to read in lattice vectors from line: \n");
                 put_string( stdout, &ichar[0],100);
                 return -15;
               }
          }
        printf("abc= %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n", 
                    *p_abc, *(p_abc+1), *(p_abc+2), *(p_abc+3), *(p_abc+4), *(p_abc+5));
     }

/* Read in the atomic data */

   at_end=0;
   *p_num_atoms=0;
   *p_num_of_mols=0;
   p_atom= p_molecule;
   *p_num_mol_members=0;

   debug=FALSE;
   while (at_end != 2) 
     {
       num_of_chars= read_line( fp, &ichar[0]);

       at_end= read_atom_data(&ichar[0],p_num_atoms, at_end, num_of_chars,
                              p_num_of_mols, p_atom, 
                              p_mol_number, just_count);
      
       if (at_end == 0)
         {
           ++*p_num_mol_members;
           if (!just_count) p_atom++;
//           p_mol_number++;
         }
     
       if (at_end == 2)
         {
           p_num_mol_members++;
           *p_num_mol_members=0;
         }
    }

   if (just_count) printf("All atoms counted\n");
             else  printf("All atoms read\n");

/***********************************************/
/*** Look for GRUP group for intervasp *********/
/*** Added May 04 DJW                  *********/
/***********************************************/

   if (have_grp && !just_count)
     {
       p_atom= p_molecule;

       p_groups->num_grp1 = -1;
       p_groups->num_grp2 = -1;
       *p_num_groups = -1;
       
       for ( iloop = 0; iloop < *p_num_atoms; iloop++)
         {
           if (strcmp(p_atom->group, "GRP1") == 0 || strcmp(p_atom->group,"TOR1") == 0 
                                                           || strcmp(p_atom->group, "GRUP") == 0)
              {
                if (p_groups->num_grp1 == -1) ++(*p_num_groups);
                (p_groups->num_grp1)++;
                p_groups->group1[p_groups->num_grp1] = iloop;
              }
           if (strcmp(p_atom->group, "GRP2") == 0  || strcmp(p_atom->group,"TOR2") == 0 )
              {
                if (p_groups->num_grp2 == -1) ++(*p_num_groups);
                (p_groups->num_grp2)++;
                p_groups->group2[p_groups->num_grp2] = iloop;
              }
           p_atom++;
         }

       if (*p_num_groups == -1)
         {
           printf("ERROR : No GRUP,GRP1 or GRP2 group members found in car");
           printf(" file when they were expected!\n");
           exit(0);
         }
     }
/***********************************************/
/*** Look for fixed  atoms if requested    *****/
/*** These are flagged by the first letter *****/
/*** of the group label being "F", then    *****/
/*** the other three letters are T/F flags *****/
/*** for a,b,c direction fixing.           *****/
/*** So that FTFT would fix only the       *****/
/*** b-vector fractional co-ordinate in    *****/
/*** VASP.                                 *****/
/*** Added Nov. 06 Dave Willock            *****/
/***********************************************/

   if (find_fixed && !just_count)
     {
       p_atom= p_molecule;
       for ( iloop = 0; iloop < *p_num_atoms; iloop++)
         {
           if (strncmp(&(p_atom->group[0]), "F", 1) == 0)
              {
                p_fix_flags->fx =FALSE;
                p_fix_flags->fy =FALSE;
                p_fix_flags->fz =FALSE;
                if (strncmp(&(p_atom->group[1]), "F", 1) == 0) p_fix_flags->fx = TRUE;
                if (strncmp(&(p_atom->group[2]), "F", 1) == 0) p_fix_flags->fy = TRUE;
                if (strncmp(&(p_atom->group[3]), "F", 1) == 0) p_fix_flags->fz = TRUE;
              }
           else
              {
                p_fix_flags->fx =FALSE;
                p_fix_flags->fy =FALSE;
                p_fix_flags->fz =FALSE;
              }
           p_atom++; 
           p_fix_flags++;
         }
     }

//  --*p_num_atoms;
//  *p_num_of_mols--;

   debug=FALSE;
   printf("Returning from read_car.. with %d atoms %d mols..\n", *p_num_atoms, *p_num_of_mols);
  return 0;
}
