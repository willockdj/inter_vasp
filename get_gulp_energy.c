#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include "structures.h"
#include "maxima.h"
#include "global_values.h"

/* protype list for this routine */

int read_line(FILE *fp, int *p_char);

int locate_string( char *p_key, int *p_char, int num_of_chars );

int read_atom_data_glp( int *p_ichar, int num_of_chars, atom *p_atom, int *p_mol_number);

double get_doub(int *p_char, int num_of_chars, int *p_place, int *p_itsanum);

int get_int(int *p_char,int *point_j, int *point_itsa,
                             int *point_ndigi, int max_chars, int *sign);

void string_to_int(char *p_char2, int *p_char1, int max_position );

int next_none_space( int *p_ichar, int start, int num_of_chars );

int next_space( int *p_ichar, int start, int num_of_chars );

/*---------------------------------------------------------------------------*/

/* read in the energy from a gulp output file */

double get_gulp_energy( FILE *fp, int is_opt)
{
  int ichar[LINESIZ];
  int place,itsanum;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int idum[10], ichar2;
  int found_all, at_end, sign, ndigi;
  int start,start2,start3,start4;
  int is_energy, is_units, is_cart;
  int num_keys, good_line, iatom;
  int got_title, is_title;
  int ishell, is_spec, ispec;
  int num_occurances, occurances;

  double energy;

  char *p_key_energy;
  char *p_key_units;
  char *p_key_cart;
  char *p_key_title;
  char *p_key_species;
  char *p_char;
  char *p_key_shell;
  char *p_Al_elem;
  char *p_Mg_elem;
  char *p_Si_elem;
  char *p_O_elem;
  char *p_Na_elem;
  char  cdummy[5];

  atom *p_atom;
  atom *p_this_shell;

  charge_list *p_this_spec;

/* read in first line and check for BIOSYM header */

     num_keys=1;
     p_key_energy= "Total lattice energy";
     p_key_units = "eV";

     found_all=0;
     num_of_chars=0;
     num_occurances=1;
     if (is_opt) num_occurances=2;

     occurances=0;
     while ( (found_all < num_keys || occurances < num_occurances) && num_of_chars != -10 )
       { 

/*******************************************************************/
/***** Avoid reading a new line if the last thing done requires ****/
/***** a read failure                                           ****/
/*******************************************************************/

          num_of_chars = read_line( fp, &ichar[0]);
      
          is_energy = locate_string( p_key_energy,  &ichar[0], num_of_chars);
          is_units  = locate_string( p_key_units ,  &ichar[0], num_of_chars);

          if (is_energy && is_units)
            {
               place=0;
               energy= get_doub(&ichar[0], num_of_chars, &place, &itsanum);
               occurances++;
               found_all++;
            }
       }

  if (num_of_chars == -10) 
    {
      printf ("Failed to get energy from gulp.out\n");
      exit(0);
    }
  return energy;
}
