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

   int read_atom_data(int *p_ichar, int *p_num_atoms, int at_end,
                      int num_of_chars, int *p_next_mol, atom *p_atom,
                      int *p_mol_number, int just_count)

{
int is_end,idummy,place,itsanum,iloop,ndigi,sign;
int idum[10], igroup;
char *p_key;

/* check for end of molecule or file */
 
 p_key= "end";
 is_end= locate_string( p_key, p_ichar, num_of_chars);

 if (is_end && !at_end )
   {
/* first end encountered for a while, just end of molecule */
 
//     (*p_next_mol)++;
     return 1;
   }
 else if (is_end && at_end )
   {
/* the end of the atom list */

     return 2;
   }
 else
   {
/* a normal atom data line */

    if (just_count)
      {
        (*p_num_atoms)++;
        *p_mol_number= *p_next_mol; 
      }
    else
      {
/* copy over atom name */

        copy_int( p_ichar, &idum[0], 0, 3); 
        int_to_string(&idum[0], &(p_atom->label[0]), 4);

/* get x,y,z co-ordinates */

        place= 1;
        place= next_space (p_ichar, place, num_of_chars );

        p_atom->x = get_doub(p_ichar, num_of_chars, &place, &itsanum);
        p_atom->y = get_doub(p_ichar, num_of_chars, &place, &itsanum);
        p_atom->z = get_doub(p_ichar, num_of_chars, &place, &itsanum);

/* get group label */ 

        place= next_none_space( p_ichar, place, num_of_chars );

        if (place == -1 )
          {
            printf("Read error whilst trying to get the group from car file line:");
            put_string(stdout, p_ichar,200);
            return -1;
          }
      
        copy_int( p_ichar, &idum[0], place, place+3); 
        int_to_string(&idum[0],&(p_atom->group[0]), 4);

        place= next_space (p_ichar, place, num_of_chars );

/******* Get the group number, which can contain a string!! ***********/

        place= next_none_space( p_ichar, place, num_of_chars );

        if (place == -1 )
          {
             printf("Read error whilst trying to get the group number from car file line:");
             put_string(stdout, p_ichar,200);
             return -1;
          }

        copy_int( p_ichar, &idum[0], place, place+3);
        int_to_string(&idum[0],&(p_atom->group_no[0]), 4);

        place= next_space (p_ichar, place, num_of_chars );

/* get potential type */

        place= next_none_space( p_ichar, place, num_of_chars );

        if (place == -1 )
          {
            printf("Read error whilst trying to get the potential type from car file line:");
            put_string(stdout, p_ichar,200);
            return -1;
          }

        copy_int( p_ichar, &idum[0], place, place+2);
        int_to_string( &idum[0], &(p_atom->pot[0]), 2);

        place= next_space (p_ichar, place, num_of_chars );

/* get element type */

        place= next_none_space( p_ichar, place, num_of_chars );

        if (place == -1 )
          {
            printf("Read error whilst trying to get the element type from car file line:");
            put_string(stdout,p_ichar,200);
            return -1;
          }

        copy_int( p_ichar, &idum[0], place, place+1); 
        int_to_string( &idum[0], &(p_atom->elem[0]), 2);

        if (p_atom->elem[1] == ' ') p_atom->elem[1] = '\0';
 
        place= next_space (p_ichar, place, num_of_chars );

/* get partial charge */

        p_atom->part_chge= get_doub(p_ichar, num_of_chars, &place, &itsanum);

        if (!itsanum)
          {
            printf("Failed to find partial charge on car file line: \n");
            put_string(stdout, p_ichar,200);
          }

/* sort out molecule number and increment atom counter */

        *p_mol_number= *p_next_mol; 
        (*p_num_atoms)++;
       }
     }
 
  return 0;
     
}
