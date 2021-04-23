/***************************************************************************/
/*** Read in a VASP OUTCAR *************************************************/
/*** Dave Willock & Ed Jeffery May 04 **************************************/
/***                                  **************************************/
/*** Updated Feb./Mar. 09 Dave Willock and Kara Howard *********************/
/***                                                                     ***/
/*** Kara's update for reading energy of converged structures.           ***/
/***                                                                     ***/
/*** Altered recognition of end of labels list in OUTCAR by reading all  ***/
/*** POTCAR lines until a VRHFIN input is encountered. This was required ***/
/*** because simulation of deuterated vibrational spectra required an    ***/
/*** additional POTCAR file with type H but mass D. Original method of   ***/
/*** looking for repeating labels fails for this case.                   ***/
/***************************************************************************/
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

void cut_to_uscore( char *p_word );

void find_line(FILE *fp, char *p_key, char *p_key2, int sep, int *p_found, int max_lines);

void find_line_with_stopper(FILE *fp, char *p_key, char *p_key2, char *p_stopper, int sep, int *p_found, int max_lines);

void latt_vecs_from_cart( double *p_latt_vec, double *p_recip_latt_vec,
                            double *p_abc);

int read_force_block( FILE *fp, e_vec *p_forces, int num_atoms);
/*---------------------------------------------------------------------------*/

/* read in a outcar file for VASP movies */

int read_outcar( FILE *fp, atom *p_molecule, double *p_latt_vec,
                 double *p_recip_latt_vec, double *p_abc, double *p_eigenvals, 
                 e_vec *p_eigenvecs, e_vec *p_forces, e_vec *p_chain_forces,
                 int *p_num_atoms, atom_number *p_types, int *p_num_types,
                 int *p_num_modes, int need_freq, int need_force,
                 int need_energy, int *p_have_band, double *p_energy,
                 int need_fermi, double *p_fermi, int *p_num_energies, int just_count )

{
  int iloop, jloop, skip, iatom;
  int noskip=FALSE;
  int num_to_skip;
  int itype, imode;
  int have_all_labels;
  int sep,masssep, found;

  double *p_this_latt_vec;
  double *p_this_eigenval;
  double freq_sign;
  double type_mass[MAXTYPES];
 
  char *vasp_Ver;
  char *p_key, *p_key2, *p_key3, *p_stopper;
  char *p_masskey, *p_masskey2, *p_massstopper;
  char *tok;
  char *p_label;

  atom *p_atom;

  e_vec *p_this_eigenvec;

/* DEBUG */
  printf("Entered read_outcar routine ");

  if (just_count) printf("just counting this time\n");
  else            printf("will read data this time\n");
  printf("Will see end of file from >>%s<<\n", END_OF_INPUT);
/* get all atom labels from POTCAR lines */
/* up to stopper                         */

  p_key= "VRHFIN";
  p_key2= "none";
  p_stopper= "POSCAR";

  p_masskey = "POMASS";
  p_masskey2 = "none";

  sep = 0;
  masssep = 0;

  skip = TRUE;
  *p_num_types=-1;

  have_all_labels = FALSE;
  tok= tok_get( fp, skip, FALSE);
  vasp_Ver=tok;
  printf("running on VASP version >>%s<<\n", vasp_Ver);

/* read in set of labels until we see the stopper keyword    */
/* which will return found=FALSE fron find_line_with_stopper */
/*  */
/* one additional label so decrement num_types at end.       */

  while (!have_all_labels)
    {
      found=FALSE;
      find_line_with_stopper( fp, p_key, p_key2, p_stopper, sep, &found, -1 );
          
      if (found)
        {   
          tok= tok_get( fp, FALSE, FALSE); 
          p_label = tok; 
          printf("Debug --- p_label >>%s<< \n", p_label);
/* remove leading equal sign if present */

	  if (*p_label == '=')
            {
              printf("Found equals sign in label moving pointer\n");
              p_label++;
              printf("label=>>%s<<\n",p_label);
            }
          (*p_num_types)++;
          strcpy((p_types+*p_num_types)->atom_type, p_label);
/*	  if ((strcmp(p_label, "PAW")) == 0)
             {
               tok= tok_get( fp, FALSE, FALSE); 
               p_label = tok; 
               printf("Debug2 --- p_label >>%s<< \n", p_label);
             }
          tok= tok_get( fp, FALSE, FALSE);
*/
        }
      else
        {
/*        (*p_num_types)--; */
          printf("number of atoms types >>%d<<\n",*p_num_types);
          have_all_labels=TRUE;          
        }
/****************************************************************/
/*** Find mass of atoms *****************************************/
/****************************************************************/

  if (!have_all_labels)
    {
      found = FALSE;
      find_line( fp, p_masskey, p_masskey2, masssep, &found, -1 );  

      if (found)
        {
          tok=tok_get( fp, noskip, FALSE);
          tok=tok_get( fp, noskip, FALSE);
          type_mass[*p_num_types] = atof(tok);
          printf("Mass = %f\n", type_mass[*p_num_types]);
        }
    }
  }

/************/
/*find VRHFIN*/
/*************/
/*    
  p_key = "VRHFIN";
  p_key2 = "none";
  sep = 0;
  p_stopper = "TITEL";
  skip = TRUE;
  *p_num_types=-1;
  have_all_labels = FALSE;
  tok= tok_get( fp, skip, FALSE);
  while (!have_all_labels)
  {
    found=FALSE;
//    find_line( fp, p_key, p_key2, sep, &found, -1 );  
//      if (found)
//        {
//          tok= tok_get( fp, FALSE, FALSE); 
//          tok= tok_get( fp, FALSE, FALSE); 
      find_line_with_stopper( fp, p_key, p_key2, p_stopper, sep, &found, -1 );

         if (found)
         {
          tok= tok_get( fp, FALSE, FALSE); 
          p_label = tok;
          printf("debug --- >>%s<< \n",p_label); 

*/
/* remove leading equal sign if present */
/*
	  if (*p_label == '=')
            {
              printf("Found equals sign in label moving pointer\n");
              p_label++;
              printf("label=>>%s<<\n",p_label);
            }
*/        
/* strip off any underscores */
/*
          cut_to_uscore(p_label); 

          strcpy((p_types+*p_num_types)->atom_type, p_label);
          printf("DEBUG --- label >>%s<<\n", p_label);
          (*p_num_types)++;
         }
        else
         { 
          (*p_num_types)--;
          have_all_labels=TRUE;
         }
    }
*/
/****************************************************************/
/*** Find ions per type info ************************************/
/****************************************************************/

  p_key  = "ions";
  p_key2 = "per";
  sep = 0;

  found=FALSE;
  find_line( fp, p_key, p_key2, sep, &found, -1 );

/*** Skip to numbers on line ****/

  for ( iloop=0; iloop <=1; iloop++) tok= tok_get( fp, FALSE, FALSE);

  printf("In OUTCAR read the following ion types with the following numbers:\n");

/*** When reading number of types note that OUTCAR has real count of types ****/
/*** while inter_vasp wants the highest index for types list, i.e. knock   ****/
/*** one off here.                                                         ****/

  *p_num_atoms=0;
  for ( iloop=0; iloop <= *p_num_types; iloop++)
    {
      (p_types+iloop)->num = atoi(tok_get(fp, FALSE, FALSE))-1;
      printf("%2d %2s %3d\n", iloop, (p_types+iloop)->atom_type, (p_types+iloop)->num);
      *p_num_atoms += (p_types+iloop)->num+1;
    }
//  Note malloced version requires num_atoms to be counter rather than upper index
//  Dave Willock Nov 2017
//  (*p_num_atoms)--;
  printf("num_atoms = %d\n", *p_num_atoms);


/* get lattice vectors */

  p_key= "direct";
  p_key2= "lattice";
  sep=0;

  found=FALSE;
  find_line( fp, p_key, p_key2, sep, &found, -1 );
  
/* get lattice vector components */

   p_this_latt_vec = p_latt_vec;
   for ( iloop = 0; iloop <= 2; iloop++ )
     {  
       skip = TRUE;
/* read numbers */
       for ( jloop = 0; jloop <= 2; jloop++ )
          {
/* Read next line */
            tok= tok_get( fp, skip, FALSE);
            skip=FALSE;

            *p_this_latt_vec = atof(tok);
            printf("component %d %d %10.6f\n", iloop, jloop, *p_this_latt_vec);
            p_this_latt_vec++;
          }
     }

skip = TRUE;

latt_vecs_from_cart( p_latt_vec, p_recip_latt_vec, p_abc );

printf("Got all Lattice vectors\n");

/************************************************************/
/*** Find cartesian co-ordinates of ions  *******************/
/************************************************************/
  p_key= "position";
  p_key2= "cartesian";
  sep = 3;

  found=FALSE;
  find_line( fp, p_key, p_key2, sep, &found, -1 );

/*****************************************/
/*** Read in and label atoms *************/
/*****************************************/

  if (!just_count)
    {
      p_atom=p_molecule;
      for ( itype=0; itype <= *p_num_types; itype++)
        {
          for ( iatom=0; iatom <= (p_types+itype)->num; iatom++)
            {
              p_atom->part_chge= 0.0;
              strcpy(p_atom->pot, "?"); 
              strcpy(p_atom->group, "VASP"); 
              strcpy(p_atom->group_no, "X"); 
              strcpy(p_atom->elem, (p_types+itype)->atom_type);
              sprintf(p_atom->label,"%s%d",  (p_types+itype)->atom_type, iatom+1 );
        
              p_atom->x = atof(tok_get(fp, skip, FALSE));
              p_atom->y = atof(tok_get(fp, FALSE, FALSE));
              p_atom->z = atof(tok_get(fp, FALSE, FALSE));

              p_atom->mass = type_mass[itype];

              printf("%s %10.6f %10.6f %10.6f %s\n", p_atom->label,
                                                     p_atom->x,
                                                     p_atom->y,
                                                     p_atom->z,
                                                     p_atom->elem);
              p_atom++;
           }
      }
  }
/*****************************************************/
/*** Find Normal mode data ***************************/
/*****************************************************/

 if (need_freq)
   {
     p_key= "Degrees";
     p_key2= "freedom";
     sep = 1;

     found=FALSE;
     find_line( fp, p_key, p_key2, sep, &found, -1 );

/* move along line to numbers */
     for ( jloop = 0; jloop <= 2; jloop++ ) tok= tok_get( fp, FALSE, FALSE);

     *p_num_modes = atoi(tok);
     printf("degrees of freedom = %d\n",*p_num_modes);

     if (*p_num_modes > MAXMODES)
        {
          printf("ERROR: Number of modes in OUTCAR file ");
          printf("exceeds maximum of %d for inter_vasp\n", MAXMODES);
          exit(0);
        }

//  printf("vasp vers >>%s<\n", vasp_Ver); 
//  printf("value = >>%d<<\n", strcmp(vasp_Ver, "vasp")); 
//     if ((strncmp(VV, "vasp.4.", 6)) == 0)
       if ((strncmp(vasp_Ver, "vasp.4.6", 6)) == 0)
        {
          p_key= "Eigenvectors";
          p_key2= "SQRT(mass)";
          sep = 3;
        }
     else
        {
          p_key= "Eigenvectors";
          p_key2= "the";
          sep = 3;
        }

//debug
     printf("pkey >>%s<< \n", p_key);
     printf("pkey2 >>%s<< \n", p_key2);

     found=FALSE;
     find_line( fp, p_key, p_key2, sep, &found, -1 );

/* jump lines to just before first mode */
     for ( jloop = 0; jloop <= 1; jloop++ ) tok= tok_get( fp, skip, FALSE);

     p_this_eigenval = p_eigenvals;
     p_this_eigenvec = p_eigenvecs;
     for (imode = 0; imode < *p_num_modes; imode++)
       {
         tok= tok_get( fp, skip, FALSE);
         tok= tok_get( fp, skip, FALSE);
         printf("header : %s\n", tok);
         tok= tok_get( fp, noskip, FALSE);
  
/************************************************/
/*** look out for imaginary mode flag ***********/
/************************************************/
         if (strcmp(tok,"f/i=") == 0)
           {
              num_to_skip=4;
              freq_sign = -1.0;
           }
         else
           {
              num_to_skip=5;
              freq_sign = 1.0;
           }
 
         for ( jloop = 0; jloop <= num_to_skip; jloop++ ) tok= tok_get( fp, noskip, FALSE);
         *p_this_eigenval= freq_sign * atof(tok);
         printf("Mode frequency %10.6f\n", *p_this_eigenval);

         tok= tok_get( fp, skip, FALSE);
         for (iatom = 0; iatom < *p_num_atoms; iatom++)
           {
              tok_get( fp, skip, FALSE);
              for ( jloop = 0; jloop <= 1; jloop++ ) tok= tok_get( fp, noskip, FALSE);

              p_this_eigenvec->dx[iatom] = atof( tok_get( fp, noskip, FALSE ));
              p_this_eigenvec->dy[iatom] = atof( tok_get( fp, noskip, FALSE ));
              p_this_eigenvec->dz[iatom] = atof( tok_get( fp, noskip, FALSE ));

              printf("Atom %d %10.6f %10.6f %10.6f\n", iatom, 
                                                    p_this_eigenvec->dx[iatom],
                                                    p_this_eigenvec->dy[iatom],
                                                    p_this_eigenvec->dz[iatom] );
                                        
           } 
        p_this_eigenval++;
        p_this_eigenvec++;
       }
    }

  if (need_force)
    {
/************************************************************************************/
/*** Read file entries repeatedly so we end with final forces ***********************/
/************************************************************************************/
	 
      p_key= "POSITION";
      p_key2= "TOTAL-FORCE";
      p_key3= "TANGENT";
      sep = 0;

      found=TRUE;
      while (found)
	{
           printf("Looking for forces\n");
           find_line( fp, p_key, p_key2, sep, &found, -1 );

	   if (found)
            {
/***************************************************************/
/*** In the elastic band OUTCARs the vector for the band      **/
/*** tangent is reported in the same space as the POSITION    **/
/*** data in the TOTAL-FORCE output.                          **/
/*** The data reported under CHAIN-FORCE is actually the      **/
/*** force due to the chain springs along the chain direction **/
/***************************************************************/
                *p_have_band = FALSE;

                read_force_block( fp, p_forces, *p_num_atoms);

                printf("First force component : %10.6f\n", p_forces->dx[0]);

/****** Look for CHAIN forces to detect if we have an elastic band calc *****/
                find_line( fp, p_key3, "none", sep, p_have_band, 20 );

                if ( *p_have_band )
                  {
          printf("This was an elastic band calculation, reading CHAIN-FORCES\n");
                    read_force_block( fp, p_chain_forces, *p_num_atoms);

                  }
            }
       }
    }
/************************************************************************************/
/*** Read file entries repeatedly so we end with final energy ***********************/
/************************************************************************************/
    
/*      fseek(fp,0,SEEK_SET);

      p_key= "energy";
      p_key2= "energy(sigma->0)";
      sep = 3;

      found=TRUE;
      while (found)
        {
           printf("Looking for energy\n");
           find_line( fp, p_key, p_key2, sep, &found, -1 );

           if (found)
            {
              tok_get( fp, noskip, FALSE);
              *p_energy = atof( tok_get( fp, noskip, FALSE ));
            }

        }
printf("Energy = %f\n", *p_energy);*/

  if (need_energy)
    {
      p_key= "FREE";
      p_key2= "ENERGIE";
      sep = 0;
      skip = TRUE;

      found=TRUE;
      while (found)
        {
           printf("Looking for energy....\n");
           find_line( fp, p_key, p_key2, sep, &found, -1 );

           if (found)
            {
              printf("EDIFF reached\n");
              /*** Skip lines to get to the energy ***/
              for ( jloop = 0; jloop <= 3; jloop++ ) tok= tok_get( fp, skip, FALSE);
              for ( jloop = 0; jloop <= 4; jloop++ ) tok= tok_get( fp, FALSE, FALSE);

              if (just_count)
                {
                   printf("...waiting...1..\n");
                   ++*p_num_energies;
                   printf("Just counting, now have %d energies.\n", *p_num_energies);
                }
              else
                {
                   printf("...waiting...2..\n");
                   *p_energy = atof( tok_get( fp, noskip, FALSE ));
                   printf("Energy = %f\n", *p_energy);
                   p_energy++;
                }
            }
        }
              

    }

/*** Debug Jan 2012 introduced while loop to find final instance of E-fermi rather than first ***/
/*** Dave Willock ***/
  printf("Deciding on Fermi search.....\n");
  if (need_fermi)
    {
      p_key= "E-fermi";
      p_key2= "none";
      sep = 0;

      found=FALSE;
      while (!found)
        {
           printf("Looking for Fermi energy\n");
           find_line( fp, p_key, p_key2, sep, &found, -1 );

           if (found)
            {
              printf("Fermi energy found\n");
              *p_fermi = atof( tok_get( fp, noskip, FALSE ));
            }

           printf("Fermi energy = %f\n", *p_fermi);
       }

    }
  else
    {
       printf("Fermi energy read from OUTCAR not requested\n");
    }


/*** For mallocing version need to go back to main for memory assignment before filling arrays ***/

  printf("Returning from read_outcar.c\n");
  if (just_count) return *p_num_atoms;
   else return 0;
}
