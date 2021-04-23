/* routine to decipher a .car atom data line */

#include <stdio.h>
#include <string.h>
#include <limits.h>
#include "maxima.h"
#include "structures.h"

/* ------------------------------------------------------------ */

int locate_string( char *p_key, int *p_ichar, int num_of_chars );

int next_none_space( int *p_ichar, int start, int num_of_chars );

int next_space( int *p_ichar, int start, int num_of_chars );

int get_int( int *p_ichar, int *point_j,
                               int *itsanum, int *ndigi,int i, int *sign);

double get_doub(int *p_ichar, int num_of_chars, int *p_place, int *p_itsanum);

void  copy_int( int *p_ichar1, int *p_ichar2, int min, int max);

void put_string( FILE *fp, int *p_ichar, int length);

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position );

/* ------------------------------------------------------------ */

void read_atom_data_pdb(int *p_ichar, int *p_num_atoms, int num_of_chars, atom *p_atom)
{
int is_end,idummy,place,itsanum,iloop,ndigi,sign;
int idum[10], igroup;
char *p_key;

/* a normal atom data line */

/* copy over atom name */

        copy_int( p_ichar, &idum[0], 0, 3); 
        int_to_string(&idum[0], &(p_atom->label[0]), 4);

/* get x,y,z co-ordinates */
/* in a pdb file the x co-ordinate should be after column 28 */

        place= 28;
        place= next_space (p_ichar, place, num_of_chars );

        p_atom->x = get_doub(p_ichar, num_of_chars, &place, &itsanum);
        p_atom->y = get_doub(p_ichar, num_of_chars, &place, &itsanum);
        p_atom->z = get_doub(p_ichar, num_of_chars, &place, &itsanum);

/* get group label which start column 18 */ 

        place=17;

        copy_int( p_ichar, &idum[0], place, place+3); 
        int_to_string(&idum[0],&(p_atom->group[0]), 4);


/******* Get the group number, which can contain a string!! ***********/

        place= 6;

        copy_int( p_ichar, &idum[0], place, place+3);
        int_to_string(&idum[0],&(p_atom->group_no[0]), 4);

/* get element type */

        place= 72;

        place= next_none_space( p_ichar, place, num_of_chars );

        if (place == -1 )
          {
            printf("Read error whilst trying to get the element type from pdb file line:");
            put_string(stdout,p_ichar,200);
            return -1;
          }

        copy_int( p_ichar, &idum[0], place, place+1); 
        int_to_string( &idum[0], &(p_atom->elem[0]), 2);

        if (p_atom->elem[1] == ' ') p_atom->elem[1] = '\0';
 
        place= next_space (p_ichar, place, num_of_chars );

/*** Make up the rest **/
        p_atom->part_chge= get_doub(p_ichar, num_of_chars, &place, &itsanum);
        strcpy(p_atom->pot, "?");

        return;
  }
 
     
