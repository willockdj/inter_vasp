/********************************************************************/
/*** read_pdb reads in a .pdb file formatted structure **************/
/*** recognises group labels:                          **************/
/*** GRUP : Atoms that are part of the interpolation   **************/
/***        group.                                     **************/
/*** GRP1 : GRP1 and GRP2 labels introduced to allow   **************/
/***        interpolations with a bond between two     **************/
/*** GRP2 : molecular groups being broken.             **************/
/*** FIXX : Atoms that should be fixed: if requested   **************/
/***                                                   **************/
/*** last updated April 3rd 2020, Dave Willock         **************/
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

void read_atom_data_pdb(int *p_ichar, int *p_num_atoms, int num_of_chars, atom *p_atom);

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

void put_string(FILE *fp, int *p_ichar, int length);
/*---------------------------------------------------------------------------*/

/* read in a .car insight formatted file */

int read_pdb( FILE *fp, int *p_header_line, int *p_title_line, 
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms, 
              double *p_abc, int *p_been_before, int just_count)
{
  int ichar[LINESIZ];
  int place,itsanum;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int at_end;
  int start,start2,start3,start4;
  int *p_this_index;
  int have_atom, num_lines;

  char *p_key;

  atom *p_atom;

  *p_num_atoms=0;
  p_atom= p_molecule;

  while ( (num_of_chars = read_line(fp, &ichar[0])) != -10)
    {
       num_lines++;
       put_string( stdout, &ichar[0], num_of_chars); 
       printf("  line has %d characters\n", num_of_chars);

// Look for crystal lattice information
       p_key = "CRYST1";
       if (locate_string( p_key, &ichar[0], num_of_chars))
         {
            *p_pbc= TRUE;
            place= 7;
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

// Check for atom lines
       p_key = "HETATM";
       have_atom = locate_string( p_key, &ichar[0], num_of_chars);
       p_key = "ATOM";
       have_atom = have_atom || locate_string( p_key, &ichar[0], num_of_chars);

       if (have_atom)
         {
            if (just_count) (*p_num_atoms)++;
            else
              { 
                 (*p_num_atoms)++;
                 read_atom_data_pdb(&ichar[0],p_num_atoms, num_of_chars, p_atom);
                 p_atom++;
              }
         }
    }


  return 0;
}
