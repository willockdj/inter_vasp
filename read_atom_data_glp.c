/* routine to decipher a .glp atom data line */

#include <stdio.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"
#include "global_values.h"

/* ------------------------------------------------------------ */

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int next_none_space( int *p_ichar, int start, int num_of_chars );

int next_space( int *p_ichar, int start, int num_of_chars );

int get_int( int *p_ichar, int *point_j,
                               int *itsanum, int *ndigi,int i, int *sign);

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

void  copy_int( int *p_ichar1, int *p_ichar2, int min, int max);

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position );

/* ------------------------------------------------------------ */

int read_atom_data_glp( int *p_ichar, int num_of_chars, atom *p_atom, int *p_mol_number )

{
int is_end,idummy,place,itsanum,iloop,ndigi,sign;
int idum[10], igroup, good_line;
char *p_key;

/* read in atom data from a single normal atom data line */

   good_line= FALSE;

/* copy over atom name */

   copy_int( p_ichar, &idum[0], 0, 3); 
   int_to_string(&idum[0], &(p_atom->label[0]), 4);

/* get a,b,c fractional co-ordinates */

   place= 1;
   place= next_space (p_ichar, place, num_of_chars );

/* get core or shell label not working yet! */
   place= next_none_space( p_ichar, place, num_of_chars );
   place= next_space (p_ichar, place, num_of_chars );

   p_atom->x = get_doub(p_ichar, num_of_chars, &place, &itsanum);
   p_atom->y = get_doub(p_ichar, num_of_chars, &place, &itsanum);
   p_atom->z = get_doub(p_ichar, num_of_chars, &place, &itsanum);

   if (itsanum) good_line = TRUE;

/* Read in charges *****************/
/* if present      *****************/

   p_atom->part_chge  = get_doub(p_ichar, num_of_chars, &place, &itsanum);
   printf("Read %10.6f  %10.6f  %10.6f  %10.6f \n",
            p_atom->x, p_atom->y, p_atom->z, p_atom->part_chge);

/* Fill in all extra data required */

   *p_mol_number= 0;

   if (good_line) return TRUE;
   return FALSE;
}
