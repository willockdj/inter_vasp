/*******************************************************************/
/*** Gulp file reader **********************************************/
/*** added ignorance of hashes feature June 03 *********************/
/*******************************************************************/

#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include <string.h>
#include "maxima.h"
#include "structures.h"
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

void int_to_string(int *p_ichar1, char *p_ichar2, int max_position );

/*---------------------------------------------------------------------------*/

/* read in a .gulp file */

int read_gulp( FILE *fp, int *p_header_line, int *p_title_line, 
              atom *p_molecule, atom *p_shells, int *p_date_line, int *p_pbc, 
              int *p_num_atoms, int *p_num_shells,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number, 
              double *p_abc, int *p_been_before, int *p_top_bit,
              int *p_num_top_chars, int *p_bottom_bit, int *p_num_bottom_chars,
              int *p_space_group, charge_list *p_spec_charges, int *p_num_species)

{
  int ichar[LINESIZ];
  int place,itsanum;
  int idave,num_of_chars,is_biosym,idummy,iloop;
  int idum[10], ichar2;
  int found_all, at_end, sign, ndigi;
  int start,start2,start3,start4;
  int is_cell, is_fract, is_cart, is_hash;
  int num_keys, good_line, iatom;
  int got_title, is_title;
  int ishell, is_spec, ispec;

  double *p_this_abc;

  char *p_key_cell;
  char *p_key_fract;
  char *p_key_cart;
  char *p_key_title;
  char *p_key_species;
  char *p_key_hash;
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

     printf("DEBUG>> In read_gulp\n");
     *p_been_before= 1;

     num_keys=3;
     p_key_cell= "cell";
     p_key_fract= "fractional";
     p_key_cart= "cartesian";
     p_key_title= "title";
     p_key_species= "species";
     p_key_shell= "shel";
     p_key_hash= "#";

     found_all=0;
     got_title= FALSE;
     num_of_chars=0;

     while ( found_all < num_keys && num_of_chars != -10 )
       { 

/*******************************************************************/
/***** Avoid reading a new line if the last thing done requires ****/
/***** a read failure                                           ****/
/*******************************************************************/

          if (!is_fract && !is_cart) num_of_chars = read_line( fp, &ichar[0]);
          printf("Read line in GULP file\n");
      
          is_hash = locate_string( p_key_hash, &ichar[0], num_of_chars);
          is_title= locate_string( p_key_title, &ichar[0], num_of_chars);
          is_cell = locate_string( p_key_cell,  &ichar[0], num_of_chars)   && !is_hash;
          is_fract= locate_string( p_key_fract, &ichar[0], num_of_chars)   && !is_hash;
          is_cart = locate_string( p_key_cart, &ichar[0], num_of_chars)    && !is_hash;
          is_spec = locate_string( p_key_species, &ichar[0], num_of_chars) && !is_hash;

     if (is_cell)
       {
         printf("Getting cell data\n"); 
         found_all++;
         num_of_chars = read_line( fp, &ichar[0]);
   
         p_this_abc=p_abc;
         place= 0;
         *p_this_abc = get_doub(&ichar[0], num_of_chars, &place, &itsanum);
         p_this_abc++;
         *p_this_abc = get_doub(&ichar[0], num_of_chars, &place, &itsanum);
         p_this_abc++;
         *p_this_abc = get_doub(&ichar[0], num_of_chars, &place, &itsanum);
         p_this_abc++;
         *p_this_abc = get_doub(&ichar[0], num_of_chars, &place, &itsanum);
         p_this_abc++;
         *p_this_abc = get_doub(&ichar[0], num_of_chars, &place, &itsanum);
         p_this_abc++;
         *p_this_abc = get_doub(&ichar[0], num_of_chars, &place, &itsanum);

         printf("Read in this_abc alpha beta gamma as:");
         printf("%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f \n",
                                       *p_abc, *(p_abc+1), *(p_abc+2), 
                                       *(p_abc+3), *(p_abc+4), *(p_abc+5));
    
       }
     else if (is_fract || is_cart )
       {
         found_all++;
         place=0; 
         if (is_fract)
            {  
               *p_space_group= get_int(&ichar[0], &place, &itsanum, 
                                                    &ndigi, num_of_chars, &sign);
               printf("Space group read as %d\n", *p_space_group);
            }

         p_atom= p_molecule;
         p_this_shell= p_shells;
         *p_num_atoms = 0;
         *p_num_shells = 0;
         good_line= TRUE;

         while (good_line)
           {
             num_of_chars = read_line( fp, &ichar[0]);
             printf("Getting atom data\n"); 

/*** put shells in separate array ****/

             if (locate_string( p_key_hash,  &ichar[0], num_of_chars))
               {
                 continue;
               }
             else if (locate_string( p_key_shell,  &ichar[0], num_of_chars)) 
               {
                 good_line= read_atom_data_glp(&ichar[0], num_of_chars, p_this_shell,
                                                                    p_mol_number );
                 p_this_shell++;
                 ++*p_num_shells; 
               }
             else
               {
                 good_line= read_atom_data_glp(&ichar[0], num_of_chars, p_atom, 
                                                                    p_mol_number ); 
                 p_atom++;
                 ++*p_num_atoms;
               }
           }

         --*p_num_atoms;
         --*p_num_shells;

         p_atom= p_molecule;
         p_this_shell= p_shells;

         p_Al_elem= "al";
         p_Si_elem= "si";
         p_O_elem= "o";
         p_Na_elem= "na";
         p_Mg_elem= "mg";
         for (iatom=0; iatom <= *p_num_atoms; iatom++)
           {
              strcpy(&cdummy[0], &(p_atom->label[0]));
              for (ichar2=0; cdummy[ichar2] != '\0'; ichar2++) cdummy[ichar2]=tolower(cdummy[ichar2]);

              if (strncmp(p_Mg_elem, &cdummy[0], 2) == 0) 
                                                    sscanf("Mg","%s", &(p_atom->elem[0]));
              if (strncmp(p_Al_elem, &cdummy[0], 2) == 0) 
                                                    sscanf("Al","%s", &(p_atom->elem[0]));
              if (strncmp(p_Si_elem, &cdummy[0], 2) == 0) 
                                                    sscanf("Si","%s", &(p_atom->elem[0]));
              if (strncmp(p_Na_elem, &cdummy[0], 1) == 0) 
                                                     sscanf("Na","%s", &(p_atom->elem[0]));
              if (strncmp(p_O_elem, &cdummy[0], 1) == 0) 
                                                     sscanf("O","%s", &(p_atom->elem[0]));
              p_atom++;
           }
         for (iatom=0; iatom <= *p_num_shells; iatom++)
           {
              strcpy(&cdummy[0], &(p_this_shell->label[0]));
              for (ichar2=0; cdummy[ichar2] != '\0'; ichar2++) cdummy[ichar2]=tolower(cdummy[ichar2]);

              if (strncmp(p_Al_elem, &cdummy[0], 2) == 0) 
                                                    sscanf("Al","%s", &(p_this_shell->elem[0]));
              if (strncmp(p_Mg_elem, &cdummy[0], 2) == 0) 
                                                    sscanf("Mg","%s", &(p_this_shell->elem[0]));
              if (strncmp(p_Si_elem, &cdummy[0], 2) == 0) 
                                                    sscanf("Si","%s", &(p_this_shell->elem[0]));
              if (strncmp(p_Na_elem, &cdummy[0], 1) == 0) 
                                                     sscanf("Na","%s", &(p_this_shell->elem[0]));
              if (strncmp(p_O_elem, &cdummy[0], 1) == 0) 
                                                     sscanf("O","%s", &(p_this_shell->elem[0]));
              p_this_shell++;
           }
       }

     else if (is_spec)
       {
         printf("Getting species data\n"); 
         found_all++;
         place=0;
         *p_num_species= get_int(&ichar[0], &place, &itsanum,
                                               &ndigi, num_of_chars, &sign);
         printf("Will read in charges for %d species\n", *p_num_species);
         --*p_num_species;

         p_this_spec= p_spec_charges;
         for (ispec=0; ispec<= *p_num_species; ispec++)
           {
              num_of_chars = read_line( fp, &ichar[0]);

              place=0;
              place= next_none_space( &ichar[0], place, num_of_chars );

              int_to_string(&ichar[place], &(p_this_spec->label[0]), 4);

              p_this_spec->is_core = !locate_string( p_key_shell,  
                                                            &ichar[0], num_of_chars);
         
              place= next_space( &ichar[0], place, num_of_chars );

              p_this_spec->part_chge= get_doub( &ichar[0], num_of_chars, 
                                                &place, &itsanum);
              
              p_this_spec++;
           }
       }
     else if (is_title)
       {
         printf("Getting title\n");
         got_title= TRUE;
         num_of_chars = read_line( fp, p_title_line );
       }
     else if (num_of_chars == -10)
       {
         printf("End of file found before gulp lines 'cell', 'species' or 'fractional'\n");
         printf("Are you sure this is a gulp file?\nStopping\n");
         exit(0);
       }
     else if (is_hash) 
       {
         continue;
       }
     else
       {
/**** Anything else needs to be remembered for regurgatation *******/
         for (iloop=0; iloop <= num_of_chars; iloop++)
           {
              *p_top_bit= ichar[iloop];
              ++*p_num_top_chars;
              p_top_bit++;
           }
         *p_top_bit= '\n';
         p_top_bit++;
         ++*p_num_top_chars;
         *p_top_bit= '\0';
       }
     }

/****************************************************/
/***** read remainder of file into bottom bit *******/
/****************************************************/

     *p_num_bottom_chars=0;

/****************************************************/
/*** if the last thing read is frac will have *******/
/*** the line that the error occured on left  *******/
/****************************************************/

     if (!is_fract && !is_cart) num_of_chars = read_line( fp, &ichar[0] );

     num_of_chars = read_line( fp, &ichar[0] );

     while (num_of_chars != -10)
       {
         for(iloop=0; iloop <= num_of_chars; iloop++)
           {
             *p_bottom_bit = ichar[iloop];
             p_bottom_bit++;
             ++*p_num_bottom_chars;
           }

         *p_bottom_bit='\n';
         p_bottom_bit++;
         ++*p_num_bottom_chars;

         num_of_chars = read_line( fp, &ichar[0] );

       }

     p_bottom_bit++;
     *p_bottom_bit='\0';
     (*p_num_bottom_chars)++;

/****************************************************/
/*** Fill in atom charges from spec list ************/
/****************************************************/

     p_atom= p_molecule;
     for (iatom=0; iatom < *p_num_atoms; iatom++)
       {
         p_this_spec= p_spec_charges;
         for (ispec=0; ispec<= *p_num_species; ispec++)
           {
              if (  strcmp(&(p_this_spec->label[0]), &(p_atom->label[0])) == 0
                 && p_this_spec->is_core) 
                  {
                    p_atom->part_chge = p_this_spec->part_chge;
                  }
              p_this_spec++;
           }
         p_atom++;
       }

     p_this_shell= p_shells; 
     for (ishell=0; ishell < *p_num_shells; ishell++)
       {
         p_this_spec= p_spec_charges;
         for (ispec=0; ispec<= *p_num_species; ispec++)
           {
              if (  strcmp(&(p_this_spec->label[0]), &(p_this_shell->label[0])) == 0
                 && !p_this_spec->is_core)
                  {
                    p_this_shell->part_chge = p_this_spec->part_chge;
                  }
              p_this_spec++;
           }
         p_this_shell++;
       }

/****************************************************/
/*** Fill in any blanks *****************************/
/****************************************************/

     if (!got_title)
       {
         p_char= "No title given\n";
         string_to_int(p_char, p_title_line, 80); 
       }

     p_atom= p_molecule;
     for (iatom=0; iatom < *p_num_atoms; iatom++)
       {
          sscanf("GRUP","%s", &(p_atom->group[0]));
          sscanf("?","%s", &(p_atom->pot[0]));
          sscanf("1","%s", &(p_atom->group_no[0]));
          p_atom++;
       }

     p_this_shell= p_shells;
     for (ishell=0; ishell < *p_num_shells; ishell++)
       {
          sscanf("GRUP","%s", &(p_this_shell->group[0]));
          sscanf("?","%s", &(p_this_shell->pot[0]));
          sscanf("1","%s", &(p_this_shell->group_no[0]));
          p_this_shell++;
       }
   
    *p_num_mol_members= *p_num_atoms;

  if (is_fract) return GULP_FRACT;
  if (is_cart)  return GULP_CART;
}
