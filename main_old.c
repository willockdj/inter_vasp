/******************************************************************************/
/*                                                                            */
/* Inter_vasp is a general purpose program for setting up VASP jobs           */
/* and for converting the results for analysis in Materials Studio.           */
/*                                                                            */
/*                                      Begun June 03, Dave Willock           */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "maxima.h"
#include "global_values.h"
#include "structures.h"
#include "constants.h"

#define MAIN 0
#include "data.h"
#include "header.h"
#undef MAIN

void open_file(FILE **p_file, char *p_filename, char *p_status);

int read_input( FILE *input_fp, char *p_title, char *p_master_input,
                int *p_is_gulp, int *p_is_car, int *p_is_vasp, 
                char *p_variable_label, char *p_end_input, 
                int *p_end_is_gulp, int *p_end_is_car, int *p_end_is_vasp, 
                double *p_temperature, double *p_min_weight, int *p_assess, 
                int *p_analyse, int *p_num_to_set, int *p_num_per_formula, 
                char *p_potcar_input, char *p_outcar_input, int *p_have_out,
                int *p_need_car, char *p_car_output, int *p_need_poscar, 
                int *p_need_arc, char *p_arc_output,
                char *p_poscar_output, int *p_need_freq, int *p_need_force,
                int *p_num_inter, char *p_grp_cnt_lab, int *p_have_grp,
                char *p_mol_cnt_lab, int *p_have_mol,
                int *p_num_steps, double *p_amplitude, int *p_linear,
                int *p_mode, int *p_need_pdb, char *p_pdb_output,
                int *p_need_late, double *p_traj_switch, int *p_need_morph );

int read_car( FILE *fp, int *p_header_line, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_abc, int *p_been_before, int *p_group_indices, 
              int *p_num_in_group, int have_group);

int read_gulp( FILE *fp, int *p_header_line, int *p_title_line,
              atom *p_molecule, atom *p_shells, int *p_date_line, int *p_pbc,
              int *p_num_atoms, int *p_num_shells,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_abc, int *p_been_before, int *p_top_bit,
              int *p_num_top_chars, int *p_bottom_bit, int *p_num_bottom_chars,
              int *p_space_group, charge_list *p_spec_charges, int *p_num_species);

int read_poscar( FILE *fp, int *p_title_line,
                 atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
                 int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
                 double *p_latt_vec, double *p_recip_latt_vec,  
                 double *p_abc, int *p_been_before, int *p_is_fract,
		 int *p_is_cart, int *p_ion_number, int *p_num_types, 
                 double *p_scale_factor);

int read_outcar( FILE *fp, atom *p_molecule, double *p_latt_vec,
                 double *p_recip_latt_vec, double *p_abc, double *p_eigenvals, 
                 e_vec *p_eigenvecs, e_vec *p_forces, int *p_num_atoms, 
                 int *p_num_types, int *p_num_modes, int need_freq, 
                 int need_force);

void generate_neighbours( atom *p_molecule, int num_atoms,
                          atom_number *p_types, int *p_num_types,
                          int use_pbc,  double *p_recip_latt_vec, 
                          double *p_latt_vec);

void sort_by_elem( atom *p_molecule, int num_atoms, atom_number *p_types, 
                   int num_types );

void cart_latt_vecs( double *p_abc, double *p_latt_vec, 
                                         double *p_recip_latt_vec);

int compare_strings( char *p_ichar1, char *p_ichar2 );

void string_to_int(char *p_ichar2, int *p_ichar1, int max_position );

void write_car( FILE *fp, int *p_header_line, int *p_title_line, 
                char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number,
                int pbc, double *p_abc, int num_atoms, double scale_factor,
                int start_frame);

void write_gulp(FILE *fp, int *p_title_line, atom *p_molecule, int num_atoms,
                atom *p_shells, int num_shells, int fract_or_cart, 
                double *p_abc,
                int *p_top_bit, int num_top_chars, int *p_bottom_bit,
                int num_bottom_chars, int space_group,
                double *p_recip_latt_vec, int *p_state,
                charge_list *p_spec_charges, int *p_num_species);

void write_poscar( FILE *fp, atom *p_molecule, double *p_fract_coords,
                   atom_number *p_types, int num_types,
                   double *p_latt_vec, double *p_scale_factor, int num_atoms,
                   int *p_title_line, char *p_c_title_line, int pbc, 
                   int is_fract);

void rotate_with_flags(atom *p_molecule, double *p_axis, double *p_origin,
                       double theta, int *p_flag_list, int num_atoms);

void vec_cross(double *p_A, double *p_B, double *p_cross);

void unit_vector(double *p_vector);

double size_vector(double *p_vector);
 
/******************************************************************/
/*** Functions specific to the job in hand ************************/
/******************************************************************/

void cart_to_fract( double cart_x,  double cart_y,  double cart_z,
                    double *fract_a, double *fract_b, double *fract_c,
                    double *p_recip_latt_vec );

void fract_to_cart( double *cart_x, double *cart_y, double *cart_z,
                    double  frac_a, double  frac_b, double  frac_c,
                    double *p_latt_vec );

int find_mol(atom *p_molecule, int num_atoms );

double atomic_mass_list( char *element );

void centre_of_mass(double *p_c_of_m, double *p_total_mass, 
                    atom *p_molecule, int num_atoms );

void move_molecule(atom *p_molecule, int num_atoms, double *move_vec);

int find_atom_label(atom *p_molecule, int num_atoms, char *p_ref_lab,
                    int *p_mol_index);

void flag_chosen_atoms( atom *p_molecule, atom *p_end_molecule,
                        int num_atoms, int *p_chosen_indices,
                        int *p_num_chosen_atoms, int mol_ind, 
                        double *p_start_cofm, 
                        double *p_end_cofm, double *p_inter_cofm);

void write_pdb(FILE *pdb_fp, atom *p_molecule, 
               double *p_abc, int num_atoms);

/******************************************************************/
/*** Start of MAIN ************************************************/
/******************************************************************/

int main(argc, argv)
  int   argc;
  char  **argv;
{

  atom molecule[MAXATOMS];
  atom end_molecule[MAXATOMS];
  atom step_molecule[MAXATOMS];
  atom chosen_one[MAXATOMS];
  atom chosen_one_end[MAXATOMS];
  atom shells[MAXATOMS];
  atom end_shells[MAXATOMS];

  atom_number types[MAXATOMS], end_types[MAXATOMS];

  bonds old_bonds[MAXBONDS], new_bonds[MAXBONDS];
  int num_old_bonds, num_new_bonds;

  e_vec eigenvecs[MAXMODES], forces;

  double eigenvalues[MAXMODES];
  double rms_x, rms_y, rms_z, rms_t;

  int num_atoms, num_shells, good_read, num_in_group, num_modes;
  int num_chosen_atoms;
  int end_num_atoms, end_num_shells, start_frame, imode;
  int header_line[LINESIZ], title_line[LINESIZ], date_line[LINESIZ];
  int end_header_line[LINESIZ], end_title_line[LINESIZ], end_date_line[LINESIZ];
  int pbc, num_of_mols, num_mol_members[MAXMOL], mol_number[MAXATOMS]; 
  int end_num_of_mols, end_num_mol_members[MAXMOL], end_mol_number[MAXATOMS]; 
  int been_before, grp_mem_ind, grp_cnt_ind, igrp, imol;
  int have_mol, mol_cnt_ind, mol_ind;
  int iii,iloop, iatom, iatom2, icomp, icoord, iletter;
  int ineigh, neigh_index, num_types, end_num_types;
  int ibond, ineigh2, still_there;
  int start_mol, num_in_this_mol;
  int ion_number[MAXTYPES];
  int num_steps, linear;
  int just_hit, now_close, num_steps_left;
  int num_steps_switch;

  double abc[6], pi, theta;
  double recip_latt_vec[9], latt_vec[9];
  double vec[3], vec1[3], size;
  double transfer_vec[3];
  double unit[3], move_dist;
  double delta_bond[MAXATOMS], new_bond;
  double orig_bond[MAXATOMS];
  double amplitude, new_bond_length;
  double latest_bond;

  double start_cofm[3];
  double end_cofm[3];
  double inter_cofm[3];
  double last_cofm[3];
  double cofm_step[3];
  double total_mass;
  double traj_switch;

  char *p_key;
  char c_title_line[LINESIZ];
  char c_header_line[LINESIZ];
  char grp_cnt_lab[7];
  char mol_cnt_lab[7];
  char variable_label[60];
  char master_input[60];
  char end_input[60];
  char potcar_input[60];
  char outcar_input[60];
  char poscar_output[60];
  char step_output[60];
  char car_output[60];
  char arc_output[60];
  char pdb_output[60];
  char frame_output[60];

/******************************************************************/
/*** Variables specific to the job in hand ************************/
/******************************************************************/

  FILE *fp_input_struct;
  FILE *fp_end_input_struct;
  FILE *fp_list_file;
  FILE *fp_control_file;
  FILE *fp_master_struct;
  FILE *fp_gulp_input;
  FILE *fp_gulp_output;
  FILE *fp_car_input;
  FILE *fp_outcar;
  FILE *fp_car_output;
  FILE *fp_arc_output;
  FILE *fp_pdb_output;
  FILE *fp_vasp_input;
  FILE *fp_vasp_output;

double fract_coords[3*MAXATOMS], end_fract_coords[3*MAXATOMS];
double inter_vec[3*MAXATOMS];
double temperature,  min_weight;
double scale_factor, step, cstep, step_inc;

int is_fract, is_cart, is_gulp, is_car, is_vasp, assess, analyse;
int end_is_gulp, end_is_car, end_is_vasp;
int num_to_set, num_per_formula;
int top_bit[500], bottom_bit[500];
int num_top_chars, num_bottom_chars, num_species;
int space_group, state;
int need_car, need_arc, need_poscar, need_interpolate, need_freq;
int need_pdb, need_late, need_morph;
int have_transfer, atom_transfered;
int need_force, which_mode;
int num_inter, have_grp, have_out, check;
int group_indices[MAXATOMS];
int chosen_indices[MAXATOMS];

charge_list spec_charges[MAXATOMS];

/*******************************************************************/
/*********** read input file (name input by user on the ************/
/*********** command line).                             ************/
/*******************************************************************/

printf("Starting\n");
if (argc >= 1)
  {
    printf("Command file is : >>%s<<\n", argv[1]);
    open_file( &fp_control_file, argv[1], "r");
  }
else
  {
    printf("Error      : No input file name given on command line\n");
    printf("Use Syntax : substitute <control_file>\n");
  }

strcpy(c_title_line, "No Title given");
strcpy(c_header_line, "No Header given");
strcpy(master_input, "No Master file Supplied");
strcpy(end_input, "No End file Supplied");
strcpy(potcar_input, "No POTCAR file Supplied");
strcpy(variable_label, "No Variable label Supplied");
title_line[0]  = -1;
header_line[0] = -1;
  date_line[0] = -1;
scale_factor = -1.0;
num_inter = 1;
num_steps = 10;
linear = FALSE;
which_mode= -1;
amplitude = 0.5;

pi = acos(-1.0);
printf("Pi is %10.6f\n", pi);

end_is_gulp= FALSE;
end_is_vasp= FALSE;
end_is_car = FALSE;
is_gulp= FALSE;
is_vasp= FALSE;
is_car = FALSE;
need_car = FALSE;
need_pdb = FALSE;
need_arc = FALSE;
need_poscar = FALSE;
need_freq = FALSE;
need_late = FALSE;
need_force = FALSE;
need_morph = FALSE;
num_to_set = 0;
num_per_formula = 0;
traj_switch = 1.5;

temperature = 0.0;
min_weight  = 0.0;

read_input( fp_control_file, &c_title_line[0], &master_input[0],
            &is_gulp, &is_car, &is_vasp, &variable_label[0],
            &end_input[0], &end_is_gulp, &end_is_car, &end_is_vasp, 
            &temperature, &min_weight, &assess, 
            &analyse, &num_to_set, &num_per_formula, 
            &potcar_input[0], &outcar_input[0], &have_out, 
            &need_car, &car_output[0], 
	    &need_poscar, &need_arc, &arc_output[0],
            &poscar_output[0], &need_freq, &need_force,
            &num_inter, &grp_cnt_lab[0], &have_grp, 
            &mol_cnt_lab[0], &have_mol,
            &num_steps, &amplitude, &linear, &which_mode,
            &need_pdb, &pdb_output[0], &need_late, &traj_switch,
            &need_morph );

need_interpolate = end_is_gulp || end_is_car || end_is_vasp;

fclose(fp_control_file);

printf("Job Title        : %s\n\n", c_title_line);
printf("Master file      : %s\n", master_input);
if (is_gulp)
   {
     printf("                           this is a gulp structure file\n");
   }
if (is_car)
   {
     printf("                           this is a car structure file\n");
   }
if (is_vasp)
   {
     printf("                           this is a VASP POSCAR structure file\n");
     printf(" The POTCAR file is : %s\n", potcar_input);
   }

if (end_is_gulp)
   {
     printf("End point file will be read from %s, this is a gulp structure file\n", end_input);
   }
if (end_is_car)
   {
     printf("End point file will be read from %s, this is a car structure file\n", end_input);
   }
if (end_is_vasp)
   {
     printf("End point file will be read from %s, this is a VASP POSCAR structure file\n", end_input);
   }

if ( need_interpolate ) 
  {
    printf("Will provide %d interpolated structures\n", num_inter);
    if (have_grp) 
      {
        printf("Will interpolate GRUP group centred of >>%s<< with special geometry rules\n"
                      ,grp_cnt_lab);
        if (!is_car) 
           {
              printf("ERROR : group centre defined but structure not given as car file\n");
              exit(0);
           }
      }
    if (have_mol)
      {
        if ( need_late )
          {
        printf("Will interpolate molecule containing");
        printf(" >>%s<< to look for late transition state.\n",mol_cnt_lab);
        printf("Using %10.6f as the multiplier of the final bondlength at which the", 
                                                                          traj_switch);
        printf(" transfer occurs\n");
          }
        else
          {
        printf("Will interpolate molecule containing >>%s<< as rigid body.\n"
                      ,mol_cnt_lab);
            if (need_morph)
              {
                printf("Internal geometry of molecule will be altered smoothly to end point\n");
              }
          }

        if (!is_car) 
           {
              printf("ERROR : molecule centre defined but structure not given as car file\n");
              exit(0);
           }
      }

  }
if (have_out) printf("OUTCAR file is   : %s\n", outcar_input);
if (need_force)
  {
     if (have_out)
       {
         printf("Will analyse OUTCAR file for atomic forces.\n");
       }
     else 
       {
         printf("ERROR : Forces requested but no OUTCAR supplied\n");
       }
  }
if (need_freq)
  {
     if (have_out)
       {
         printf("Will analyse OUTCAR file for vibrational frequencies\n");
         if (which_mode > 0)
              printf("Animation files for mode %d only requested\n", 
                            which_mode);
       }
     else 
       {
         printf("ERROR : Frequencies requested but no OUTCAR supplied\n");
       }
  }

if (need_car)
	{
	  printf("Will produce a car file called %s\n", car_output);
	}
if (need_arc)
	{
	  printf("Will produce an arc file with general stem >>%s<<\n", arc_output);
          printf("Using %d frames and oscillation amplitude %10.6f\n",
                                   num_steps, amplitude);
          if (linear) printf("frames will be linearly related, for eigen-following\n");
	}
if (need_pdb)
	{
	  printf("Will produce a pdb file named >>%s<<\n", pdb_output);
	}
if (need_poscar)
	{
	  printf("Will produce a POSCAR file called %s\n", poscar_output);
	}

/******************************************************************/
/*********** read data file (name input by user on the ************/
/*********** command line).                            ************/
/******************************************************************/

been_before= FALSE;

if (is_car)
  {
   open_file( &fp_car_input, master_input, "r");

    good_read=  read_car( fp_car_input, &header_line[0], &title_line[0],
                          &molecule[0], &date_line[0], &pbc, &num_atoms,
                          &num_of_mols, &num_mol_members[0], &mol_number[0],
                          &abc[0], &been_before, &group_indices[0], 
                          &num_in_group, have_grp);

    if (have_grp)
      {
         printf("GRUP defined as containing %d atoms:\n", num_in_group+1);
         check = FALSE;
         for (iloop=0; iloop<=num_in_group; iloop++)
           {
             grp_mem_ind = group_indices[iloop]; 
             printf("%d %d >>%s<<",iloop, grp_mem_ind, molecule[grp_mem_ind].label );
             if (strcmp(molecule[grp_mem_ind].label, grp_cnt_lab) == 0)
                {
                  check=TRUE;
                  grp_cnt_ind = grp_mem_ind;
                  printf("   the central atom\n");
                }
             else
                {
                  printf("\n");
                }
           }
         if (!check)
           {
             printf("ERROR: The central atom given in the input does not occur in GRUP\n");
             exit(0);
           }
      }
    scale_factor=1.0;

    is_cart = TRUE;
    is_fract= FALSE;
    fclose(fp_car_input);
  }
else if (is_gulp)
  {
   open_file( &fp_gulp_input, master_input, "r");

   read_gulp( fp_gulp_input, &header_line[0], &title_line[0],
              &molecule[0], &shells[0], &date_line[0], &pbc,
              &num_atoms, &num_shells, &num_of_mols, 
              &num_mol_members[0], &mol_number[0],
              &abc[0], &been_before, &top_bit[0],
              &num_top_chars, &bottom_bit[0], &num_bottom_chars,
              &space_group, &spec_charges[0], &num_species);

   scale_factor=1.0;

   printf("For GULP inputs need to set is_cart or is_fract\n");

    fclose(fp_gulp_input);
  }
else if (is_vasp)
  {
    open_file( &fp_vasp_input, master_input, "r");
    
    good_read = read_poscar( fp_vasp_input, &title_line[0],
                             &molecule[0], &date_line[0], &pbc, &num_atoms,
                             &num_of_mols, &num_mol_members[0], &mol_number[0],
                             &latt_vec[0], &recip_latt_vec[0], &abc[0], 
                             &been_before, &is_fract, &is_cart,
                             &ion_number[0], &num_types, &scale_factor);

    num_mol_members[0]=num_atoms;
    pbc=TRUE;
    fclose(fp_vasp_input);

/***** Get atom label information from the POTCAR file for the job *****/

    open_file( &fp_vasp_input, potcar_input, "r");

    good_read = read_potcar( fp_vasp_input, &molecule[0], &ion_number[0], &num_types); 
    fclose(fp_vasp_input);
  }
/********************************************************************************/
/*** Process the OUTCAR file to make convergence movie DJW & EJ May 04 **********/
/********************************************************************************/
    else if (have_out)
      {
        open_file( &fp_outcar, master_input, "r");

        read_outcar( fp_outcar, &molecule[0], &latt_vec[0], 
                     &recip_latt_vec[0], &abc[0], &eigenvalues[0],
                     &eigenvecs[0],
                     &forces, &num_atoms, &num_types, &num_modes, need_freq,
                     need_force);

        num_of_mols = 1;
        num_mol_members[0]=num_atoms+1;
        scale_factor=1.0;
        is_cart = TRUE;
        is_fract= FALSE;
        pbc=TRUE;
        fclose(fp_outcar);
      }
/********************************************************************************/
/**** End of special case stuff !! **********************************************/
/********************************************************************************/

else
  {
    printf("ERROR: Do not understand format of master file\n");
    exit(0);
  }

/*****************************************************************************************/
/*** Convert fractional co-ordinates to cartessian ***************************************/
/*****************************************************************************************/

if (is_fract)
  {
     icoord=0;

     printf("Will calculate cartessian using\n");
     printf("%10.6f %10.6f %10.6f\n", latt_vec[0], latt_vec[1], latt_vec[2]);
     printf("%10.6f %10.6f %10.6f\n", latt_vec[3], latt_vec[4], latt_vec[5]);
     printf("%10.6f %10.6f %10.6f\n\n", latt_vec[6], latt_vec[7], latt_vec[8]);

     printf("Have reciprocal space vectors\n");
     printf("%10.6f %10.6f %10.6f\n", recip_latt_vec[0], recip_latt_vec[1], recip_latt_vec[2]);
     printf("%10.6f %10.6f %10.6f\n", recip_latt_vec[3], recip_latt_vec[4], recip_latt_vec[5]);
     printf("%10.6f %10.6f %10.6f\n\n", recip_latt_vec[6], recip_latt_vec[7], recip_latt_vec[8]);

     for (iloop=0; iloop < num_atoms; iloop++)
         {

            fract_coords[icoord] = molecule[iloop].x;
            icoord++;
            fract_coords[icoord] = molecule[iloop].y;
            icoord++;
            fract_coords[icoord] = molecule[iloop].z;
            icoord++;

            printf("Sending frac to fract_to_cart: %10.6f %10.6f %10.6f\n", molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);

            fract_to_cart( &(molecule[iloop].x), &(molecule[iloop].y), &(molecule[iloop].z),
                           fract_coords[icoord-3],fract_coords[icoord-2],fract_coords[icoord-1],
                           &latt_vec[0] );

            printf("Now get cartessian: %10.6f %10.6f %10.6f\n\n", molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);
         }
   }

if (good_read == 0)
  {

/**** For car files generate cartessian lattice vectors from abc alpha beta gamma ****/

     if (pbc && is_car) cart_latt_vecs( &abc[0], &latt_vec[0], &recip_latt_vec[0]);

     printf("Found:\n %d atoms \n %d molecules \n",num_atoms,num_of_mols);

     start_mol = 0;
     for (iloop = 0; iloop < num_of_mols; iloop++)
       {
          printf("Second off to generate neighbours with start_mol= %d, num_mol_members= %d\n",
                                                 start_mol,num_mol_members[iloop]);

          generate_neighbours( &molecule[start_mol], num_mol_members[iloop]-1,
                               &types[0], &num_types,
                               pbc, &recip_latt_vec[0],  &latt_vec[0]);

          printf("Molecule %d has %d members starts at %d\n",iloop,num_mol_members[iloop], start_mol);

          for (iatom= start_mol; iatom < start_mol+num_mol_members[iloop]; iatom++)
             {
                printf("%s (elem= %s) with %d neighbours : ", 
                                molecule[iatom].label, 
                                molecule[iatom].elem, 
                                molecule[iatom].num_neigh); 

                for (ineigh=0; ineigh< molecule[iatom].num_neigh; ineigh++)
                  {
                     neigh_index= start_mol+molecule[iatom].neighb[ineigh];
                     printf("%s ",molecule[neigh_index].label);
                  }
                printf("\n");
             }
          start_mol = start_mol + num_mol_members[iloop];
       }

/*****************************************************/
/*** Put molecule numbers on each atom according to **/
/*** connectivity found.                            **/
/*****************************************************/

num_of_mols = find_mol( &molecule[0], num_mol_members[0]-1 );

for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {

     molecule[iatom].mass= atomic_mass_list(molecule[iatom].elem);

     printf("%s (elem= %s, mass=%10.6f) is in molecule %d with %d neighbours : ", 
                       molecule[iatom].label, 
                       molecule[iatom].elem, 
                       molecule[iatom].mass,
                       molecule[iatom].mol,
                       molecule[iatom].num_neigh); 

     for (ineigh=0; ineigh< molecule[iatom].num_neigh; ineigh++)
        {
          neigh_index= molecule[iatom].neighb[ineigh];
          printf("%s ",molecule[neigh_index].label);
        }
     printf("\n");
   }

/*****************************************************/
/***** Read End point file for interpolation cases ***/
/*****************************************************/
if ( need_interpolate ) 
  {
if (end_is_car)
  {
    open_file( &fp_car_input, end_input, "r");
    been_before=FALSE;

/*******************************************************/
/*** No need to re-establish group atoms at end point **/
/*** hence FALSE at end of call.                      **/
/*** This allows group defined in master only         **/
/*******************************************************/

    good_read=  read_car( fp_car_input, &end_header_line[0], &end_title_line[0],
                          &end_molecule[0], &end_date_line[0], &pbc, &end_num_atoms,
                          &end_num_of_mols, &end_num_mol_members[0], &end_mol_number[0],
                          &abc[0], &been_before, &group_indices[0], &num_in_group,
                          FALSE );

    is_cart = TRUE;
    is_fract= FALSE;
    fclose(fp_car_input);
    
  }
else
  {
    printf("Can only cope with end point car files at the moment....\n");
  }

/********************************************************************************/
/** Check consistency of end point and start point ******************************/
/********************************************************************************/

   if (end_num_atoms != num_atoms)
     {
        printf("ERROR: Trying to interpolate structures with different numbers of atoms\n");
        printf("       Counted %d at start point but %d at end.\n", num_atoms, end_num_atoms);
        exit(0);
     }

/********************************************************************************/
/** Assemble neighbour info for end point ***************************************/
/********************************************************************************/

     start_mol = 0;
     for (iloop = 0; iloop < end_num_of_mols; iloop++)
       {
          printf("First off to generate neighbours for end point\n");

          generate_neighbours( &end_molecule[start_mol], 
                               end_num_mol_members[iloop]-1,
                               &end_types[0], &end_num_types,
                               pbc, &recip_latt_vec[0],  &latt_vec[0]);

          printf("Molecule %d in end point structure has %d members starts at %d\n",
                                         iloop,end_num_mol_members[iloop], start_mol);

          for (iatom= start_mol; iatom < start_mol+end_num_mol_members[iloop]; iatom++)
             {

                end_molecule[iatom].mass= 
                                 atomic_mass_list(end_molecule[iatom].elem);

                printf("%s (elem= %s, mass %10.6f) with %d neighbours : ",
                                                   end_molecule[iatom].label, 
                                                   end_molecule[iatom].elem, 
                                                   end_molecule[iatom].mass,
                                                   end_molecule[iatom].num_neigh);

                for (ineigh=0; ineigh< end_molecule[iatom].num_neigh; ineigh++)
                  {
                     neigh_index= start_mol+end_molecule[iatom].neighb[ineigh];
                     printf("%s ",end_molecule[neigh_index].label);
                  }
                printf("\n");
             }
          start_mol = start_mol + end_num_mol_members[iloop];
       }

/***************************************************************************/
/*** Look for molecules in end structure ***********************************/
/***************************************************************************/

end_num_of_mols = find_mol( &end_molecule[0], end_num_mol_members[0]-1 );

for (iatom= 0; iatom < end_num_mol_members[0]; iatom++)
   {
     printf("%s (elem= %s) is in molecule %d with %d neighbours : ", 
                       end_molecule[iatom].label, 
                       end_molecule[iatom].elem, 
                       end_molecule[iatom].mol,
                       end_molecule[iatom].num_neigh); 

     for (ineigh=0; ineigh< end_molecule[iatom].num_neigh; ineigh++)
        {
          neigh_index= end_molecule[iatom].neighb[ineigh];
          printf("%s ",end_molecule[neigh_index].label);
        }
     printf("\n");
   }

   if (end_num_of_mols != num_of_mols)
     {
        printf("Warning: Trying to interpolate structures with different numbers of molecules\n");
        printf("         This may not be a disaster since molecules may have formed or       \n");
        printf("         fragmented so you should know if it is expected.                    \n");
     }

/*************************************************************************/
/*** Check types for start and end points ********************************/
/*************************************************************************/
   
    if ( end_num_types != num_types)
      {
         printf("ERROR: Structures for interpolation contain ");
         printf("different numbers of atom types\n");
         printf("       Have %d types at start but %d at end\n", 
                                                num_types, end_num_types);
         exit(0);
      }
    for (iloop = 0; iloop < end_num_types ; iloop++)
       {
          if ( types[iloop].num != end_types[iloop].num)
             {
                printf("Start structure has %d of atom type %s",
                          types[iloop].num, types[iloop].atom_type);

                printf(" which does not match with %d in end point\n",
                                                  end_types[iloop].num);
                exit(0);
             }
       }

/******************************************************************************/
/*** Convert fractional co-ordinates to cartessian for end point **************/
/******************************************************************************/

if (is_fract)
  {
     icoord=0;
     for (iloop=0; iloop < num_atoms; iloop++)
         {

            end_fract_coords[icoord] = end_molecule[iloop].x;
            icoord++;
            end_fract_coords[icoord] = end_molecule[iloop].y;
            icoord++;
            end_fract_coords[icoord] = end_molecule[iloop].z;
            icoord++;

            fract_to_cart( &(end_molecule[iloop].x), 
                           &(end_molecule[iloop].y), 
                           &(end_molecule[iloop].z),
                           end_fract_coords[icoord-3],
                           end_fract_coords[icoord-2],
                           end_fract_coords[icoord-1],
                           &latt_vec[0] );
         }
   }

/********************************************************************************/
/*** Carry out interpolation : output car files and POSCAR files ****************/
/********************************************************************************/

/*****************************************/
/** Sort start and end point   ***********/
/*****************************************/

         sort_by_elem( &molecule[0], num_mol_members[0], &types[0], num_types);
         sort_by_elem( &end_molecule[0], num_mol_members[0], &types[0], num_types);

/*****************************************************/
/*** locate chosen molecule if present **************/
/*****************************************************/

  if (have_mol)
    {
       check= find_atom_label(&molecule[0], num_mol_members[0]-1, 
                              &mol_cnt_lab[0], &mol_ind);

       if (!check)
         {
            printf("ERROR : Atom label %s provided to define molecule is absent\n",
                                              mol_cnt_lab);
            exit(0);
         }
       else
         {

/*****************************************************/
/*** flag chosen atoms *******************************/
/*** and re-define chosen atoms relative to the cofm */
/*****************************************************/

           flag_chosen_atoms( &molecule[0], &end_molecule[0],
                              num_mol_members[0]-1, &chosen_indices[0],
                              &num_chosen_atoms, mol_ind, 
                              &start_cofm[0], &end_cofm[0], &inter_cofm[0]);

/*****************************************************/
/*** For late transition states get bonding info *****/
/*****************************************************/

           if (need_late) 
             {
               printf("Analysing bonding before and after for late transition\n");
               printf("Bonds at start:\n");
               just_hit=TRUE;
               now_close= FALSE;
               for ( iloop=0; iloop<=num_chosen_atoms; iloop++ )
                 {
                    iatom=chosen_indices[iloop];
                    printf("%d %s : ", iatom, molecule[iatom].label);

                    for (ineigh=0; ineigh< molecule[iatom].num_neigh; ineigh++)
                      {
                         iatom2 = molecule[iatom].neighb[ineigh];
                         printf("%d %s ",iatom2, molecule[iatom2].label);
                      }
                    printf("\n");
                 } 

               num_old_bonds=-1;
               num_new_bonds=-1;
               printf("\nBonds at end:\n");
               for ( iloop=0; iloop<=num_chosen_atoms; iloop++ )
                 {
                    iatom=chosen_indices[iloop];
                    printf("%d %s : ", iatom, end_molecule[iatom].label);

                    for (ineigh=0; ineigh< end_molecule[iatom].num_neigh; ineigh++)
                      {
                         iatom2 = end_molecule[iatom].neighb[ineigh];
                         printf("%d %s ",iatom2, end_molecule[iatom2].label);
                      }

                    if ( molecule[iatom].num_neigh > end_molecule[iatom].num_neigh )
                     {
/*****************************************************/
/*** This atom has lost bonds  ***********************/
/*** Work out which ones       ***********************/
/*****************************************************/
                        for (ineigh=0; ineigh< molecule[iatom].num_neigh; ineigh++)
                          {
                            iatom2 = molecule[iatom].neighb[ineigh];

                            still_there=FALSE;
                            for (ineigh2=0; ineigh2< end_molecule[iatom].num_neigh; ineigh2++)
                              {
                                still_there = iatom2 == end_molecule[iatom].neighb[ineigh2]
                                              || still_there;
                              }

                            if (!still_there)
                              {
                                num_old_bonds++;
                                if ( num_old_bonds < MAXBONDS) 
                                  {
                                    old_bonds[num_old_bonds].atom1 = iatom;
                                    old_bonds[num_old_bonds].atom2 = iatom2;
                                  }
                                else
                                  {
                            printf("ERROR: number of old bonds exceeds maximum allowed\n");
                            exit(0);
                                  }
                              }
                          }
                       
                     }
                    else if ( molecule[iatom].num_neigh < end_molecule[iatom].num_neigh )
                     {
/*****************************************************/
/*** This atom has gained a bond *********************/
/*****************************************************/
                        for (ineigh=0; ineigh< end_molecule[iatom].num_neigh; ineigh++)
                          {
                            iatom2 = end_molecule[iatom].neighb[ineigh];

                            still_there=FALSE;
                            for (ineigh2=0; ineigh2< molecule[iatom].num_neigh; ineigh2++)
                              {
                                still_there = iatom2 == molecule[iatom].neighb[ineigh2]
                                              || still_there;
                              }

                            if (!still_there)
                              {
                                num_new_bonds++;
                                if ( num_new_bonds < MAXBONDS) 
                                  {
                                    new_bonds[num_new_bonds].atom1 = iatom;
                                    new_bonds[num_new_bonds].atom2 = iatom2;
                                  }
                                else
                                  {
                            printf("ERROR: number of new bonds exceeds maximum allowed\n");
                            exit(0);
                                  }
                              }
                          }
                     }
                    else
                     {
/*****************************************************/
/*** This atom has the same number of bonds but ******/
/*** may have changed partners !                ******/
/*****************************************************/
                        for (ineigh=0; ineigh< end_molecule[iatom].num_neigh; ineigh++)
                          {
                            iatom2 = end_molecule[iatom].neighb[ineigh];

                            for (ineigh2=0; ineigh2< molecule[iatom].num_neigh; ineigh2++)
                              {
                                if (iatom2 != molecule[iatom].neighb[ineigh2])
                                   {
                                     num_old_bonds++;
                                     num_new_bonds++;

                                     old_bonds[num_old_bonds].atom1 = iatom;
                                     old_bonds[num_old_bonds].atom2 =
                                                      molecule[iatom].neighb[ineigh2];

/**** Check not already got old bond ****************/
                                     for (ibond=0; ibond < num_old_bonds; ibond++)
                                       {
                                         if ( ( old_bonds[ibond].atom1 == iatom
                                             && old_bonds[ibond].atom2 == 
                                                            molecule[iatom].neighb[ineigh2] )
                                             || ( old_bonds[ibond].atom2 == iatom
                                             && old_bonds[ibond].atom1 == 
                                                            molecule[iatom].neighb[ineigh2] ))
                                           {
                                             num_old_bonds--;
                                           }
                                       } 

                                     new_bonds[num_new_bonds].atom1 = iatom;
                                     new_bonds[num_new_bonds].atom2 = iatom2;
                                   }
                              }
                          }
                     }

                    printf("\n");
                 } 

               printf("\nLost %d old bonds\n", num_old_bonds+1);
               for (iloop=0; iloop <= num_old_bonds; iloop++)
                 {
                   iatom = old_bonds[iloop].atom1;
                   iatom2= old_bonds[iloop].atom2;
                   printf("%d %s to %d %s\n", iatom, molecule[iatom].label,
                                              iatom2, molecule[iatom2].label);
                 }

               printf("\nGained %d new bonds\n", num_new_bonds+1);
               for (iloop=0; iloop <= num_new_bonds; iloop++)
                 {
                   iatom = new_bonds[iloop].atom1;
                   iatom2= new_bonds[iloop].atom2;
                   printf("%d %s to %d %s\n", iatom, molecule[iatom].label,
                                              iatom2, molecule[iatom2].label);

                   have_transfer=FALSE;
                   if    (iatom == old_bonds[iloop].atom1
                       || iatom == old_bonds[iloop].atom2)
                     {
                       have_transfer=TRUE;
                       atom_transfered = iatom;
                     }
                   if    (iatom2 == old_bonds[iloop].atom1
                       || iatom2 == old_bonds[iloop].atom2)
                     {
                       have_transfer=TRUE;
                       atom_transfered = iatom2;
                     }
                      
                 }

               if (have_transfer)
                 {
                   printf("This is a transfer reaction with atom %d %s transfered\n",
                                 atom_transfered, molecule[atom_transfered].label);

/**********************************************/
/** Work out the vector for the transfered ****/
/** atom                                   ****/
/**********************************************/
 
                   check=FALSE;
                   for ( iloop=0; iloop<=num_chosen_atoms; iloop++ )
                     {
                       if ( atom_transfered == chosen_indices[iloop] ) check=TRUE;
                     }
                   
                   if (check)
                     {
/**** transfered atom is from molecule     ******/
                        printf("transfer is from the molecule\n");
                        transfer_vec[0] = end_cofm[0] + end_molecule[atom_transfered].x
                                        - start_cofm[0] - molecule[atom_transfered].x;

                        transfer_vec[1] = end_cofm[1] + end_molecule[atom_transfered].y
                                        - start_cofm[1] - molecule[atom_transfered].y;

                        transfer_vec[2] = end_cofm[2] + end_molecule[atom_transfered].z
                                        - start_cofm[2] - molecule[atom_transfered].z;

/**** Get final bond length ****/
                        if ( new_bonds[0].atom1 == atom_transfered )
                          {
                            iatom = new_bonds[0].atom2;
                          }
                        else
                          {
                            iatom = new_bonds[0].atom1;
                          }

                        vec[0] = end_cofm[0]+end_molecule[atom_transfered].x
                                  - end_molecule[iatom].x;

                        vec[1] = end_cofm[1]+end_molecule[atom_transfered].y
                                  - end_molecule[iatom].y;

                        vec[2] = end_cofm[2]+end_molecule[atom_transfered].z
                                  - end_molecule[iatom].z;
                     }
                   else
                     {
/**** transfered atom is not from molecule ******/
                        printf("transfer is to the molecule\n");
                        transfer_vec[0] = end_molecule[atom_transfered].x
                                        - molecule[atom_transfered].x;

                        transfer_vec[1] = end_molecule[atom_transfered].y
                                        - molecule[atom_transfered].y;

                        transfer_vec[2] = end_molecule[atom_transfered].z
                                        - molecule[atom_transfered].z;
/**** Get final bond length ****/
                        if ( new_bonds[0].atom1 == atom_transfered )
                          {
                            iatom = new_bonds[0].atom2;
                          }
                        else
                          {
                            iatom = new_bonds[0].atom1;
                          }

                        vec[0] = end_molecule[atom_transfered].x
                                  - end_molecule[iatom].x;

                        vec[1] = end_molecule[atom_transfered].y
                                  - end_molecule[iatom].y;

                        vec[2] = end_molecule[atom_transfered].z
                                  - end_molecule[iatom].z;
                     }

                   printf("Have transfer vector %10.6f  %10.6f  %10.6f\n",
                                                transfer_vec[0],
                                                transfer_vec[1],
                                                transfer_vec[2]);

                   new_bond_length= sqrt( vec[0] * vec[0] 
                                         +vec[1] * vec[1] 
                                         +vec[2] * vec[2] );

                   printf("Final length of bond formed = %10.6f Angstroms\n",
                                new_bond_length);

                 }
               else
                 {
                   printf("late can only deal with transfer reactions at the moment\n");
                   exit(0);
                 }
             }
         }
    }
/*****************************************/
/** Get inter structure vector ***********/
/*****************************************/

         for (iloop=0; iloop < num_atoms; iloop++)
           {
             inter_vec[iloop*3]   = end_molecule[iloop].x - molecule[iloop].x; 
             inter_vec[iloop*3+1] = end_molecule[iloop].y - molecule[iloop].y; 
             inter_vec[iloop*3+2] = end_molecule[iloop].z - molecule[iloop].z; 
           }

/***********************************************/
/** Generate intermediates and output files ****/
/***********************************************/

         step= 1.0 / ( num_inter + 1 );
         printf("Interpolation step : %10.6f\n", step);

/********************************************************************/ 
/*** For the case of grouped atoms work out changes of bond length **/
/********************************************************************/ 
         if (have_grp)
           {
              for (igrp=0; igrp<=num_in_group; igrp++)
                {
                  grp_mem_ind = group_indices[igrp]; 
                  if (grp_mem_ind != grp_cnt_ind)
                    {
                       vec[0] =  molecule[grp_cnt_ind].x 
                               - molecule[grp_mem_ind].x; 

                       vec[1] =  molecule[grp_cnt_ind].y 
                               - molecule[grp_mem_ind].y; 

                       vec[2] =  molecule[grp_cnt_ind].z 
                               - molecule[grp_mem_ind].z; 

                       vec1[0] = end_molecule[grp_cnt_ind].x 
                                  - end_molecule[grp_mem_ind].x; 

                       vec1[1] = end_molecule[grp_cnt_ind].y 
                                  - end_molecule[grp_mem_ind].y; 

                       vec1[2] = end_molecule[grp_cnt_ind].z 
                                  - end_molecule[grp_mem_ind].z; 

                       orig_bond[igrp]= size_vector(&vec[0]);

                       delta_bond[igrp]= size_vector(&vec1[0])
                                        -size_vector(&vec[0]);

                   printf("centre to atom %d vector is %10.6f %10.6f %10.6f\n",
                                         grp_mem_ind,
                                         vec[0],vec[1],vec[2]);
 
                       printf("This atom is %s the centre is %s\n", molecule[grp_mem_ind].label, 
                                                                    molecule[grp_cnt_ind].label);
                       printf("Original bond length for atom  %d is %10.6f\n", grp_mem_ind, orig_bond[igrp]);
                       printf("Change in bond length for atom %d is %10.6f\n", grp_mem_ind, delta_bond[igrp]);
                    }
                }
           }

/******************************************************************/

         for (iloop=1; iloop <= num_inter+1; iloop++)
            {
               if ( have_grp )
                 {
                   for (iatom=0; iatom <= num_atoms; iatom++)
                     {
/************************************************************/
/**** Check if the atom is a member of the group defined ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
                       check = FALSE;
                       for (igrp=0; igrp<=num_in_group; igrp++)
                          {
                             grp_mem_ind = group_indices[igrp]; 
                             if ( grp_mem_ind == iatom ) check = TRUE;
                          }

                       if ( !check || iatom == grp_cnt_ind)
                         {
                            step_molecule[iatom] = molecule[iatom];
                            step_molecule[iatom].x = step_molecule[iatom].x 
                                                + step * iloop * inter_vec[iatom*3];
                            step_molecule[iatom].y = step_molecule[iatom].y 
                                                + step * iloop * inter_vec[iatom*3+1];
                            step_molecule[iatom].z = step_molecule[iatom].z 
                                                + step * iloop * inter_vec[iatom*3+2];
                         }
                     }
/****************************************************************/ 
/*** Now deal with the special geometry for the grouped atoms ***/ 
/****************************************************************/ 
                   for (igrp=0; igrp<=num_in_group; igrp++)
                     {
                       grp_mem_ind = group_indices[igrp]; 

                       printf("Using different interpolation for atom %d\n",grp_mem_ind);
                       if (grp_mem_ind != grp_cnt_ind)
                         {
                            step_molecule[grp_mem_ind] = molecule[grp_mem_ind];
/***********************************/
/*** First make the normal step ****/
/***********************************/
                            step_molecule[grp_mem_ind].x = step_molecule[grp_mem_ind].x 
                                                + step * iloop * inter_vec[grp_mem_ind*3];
                            step_molecule[grp_mem_ind].y = step_molecule[grp_mem_ind].y 
                                                + step * iloop * inter_vec[grp_mem_ind*3+1];
                            step_molecule[grp_mem_ind].z = step_molecule[grp_mem_ind].z 
                                                + step * iloop * inter_vec[grp_mem_ind*3+2];
/******************************************/
/*** Now check the vector to the centre ***/
/******************************************/
                            new_bond = orig_bond[igrp] + step * iloop * delta_bond[igrp];

                            vec[0] = step_molecule[grp_mem_ind].x 
                                             - step_molecule[grp_cnt_ind].x;
                            vec[1] = step_molecule[grp_mem_ind].y 
                                             - step_molecule[grp_cnt_ind].y;
                            vec[2] = step_molecule[grp_mem_ind].z 
                                             - step_molecule[grp_cnt_ind].z;

                            unit_vector(&vec[0]); 

                            step_molecule[grp_mem_ind].x = step_molecule[grp_cnt_ind].x 
                                                          + new_bond * vec[0]; 

                            step_molecule[grp_mem_ind].y = step_molecule[grp_cnt_ind].y 
                                                          + new_bond * vec[1]; 

                            step_molecule[grp_mem_ind].z = step_molecule[grp_cnt_ind].z 
                                                          + new_bond * vec[2]; 
                         }
                     }
                 }
               else if ( have_mol && !need_late )
                 {

                   cofm_step[0] = start_cofm[0] + step * iloop * inter_cofm[0];
                   cofm_step[1] = start_cofm[1] + step * iloop * inter_cofm[1];
                   cofm_step[2] = start_cofm[2] + step * iloop * inter_cofm[2];

                   for (iatom=0; iatom <= num_atoms; iatom++)
                     {
/************************************************************/
/**** Check if the atom is a member of the group defined ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
                       check = FALSE;
                       for (imol=0; imol<=num_chosen_atoms; imol++)
                          {
                             if ( chosen_indices[imol] == iatom ) check = TRUE;
                          }

                       if ( !check )
                         {
                            step_molecule[iatom] = molecule[iatom];
                            step_molecule[iatom].x = step_molecule[iatom].x 
                                          + step * iloop * inter_vec[iatom*3];
                            step_molecule[iatom].y = step_molecule[iatom].y 
                                          + step * iloop * inter_vec[iatom*3+1];
                            step_molecule[iatom].z = step_molecule[iatom].z 
                                          + step * iloop * inter_vec[iatom*3+2];
                         }
                     }
/************************************************************/
/*** Now interpolate the rigid molecule atoms ***************/
/************************************************************/

                   if ( need_morph )
                     {
                       for (imol=0; imol<=num_chosen_atoms; imol++)
                         {
                           iatom = chosen_indices[imol];
                           step_molecule[iatom].x = cofm_step[0] + molecule[iatom].x
                                          + step * iloop * inter_vec[iatom*3];
                           step_molecule[iatom].y = cofm_step[1] + molecule[iatom].y
                                          + step * iloop * inter_vec[iatom*3+1];
                           step_molecule[iatom].z = cofm_step[2] + molecule[iatom].z
                                          + step * iloop * inter_vec[iatom*3+2];
                         }
                     }
                   else
                     {
                       for (imol=0; imol<=num_chosen_atoms; imol++)
                         {
                           iatom = chosen_indices[imol];
                           step_molecule[iatom].x = cofm_step[0] + molecule[iatom].x;  
                           step_molecule[iatom].y = cofm_step[1] + molecule[iatom].y;  
                           step_molecule[iatom].z = cofm_step[2] + molecule[iatom].z;  
                         }
                     }
                 }
               else if ( have_mol && need_late )
                 {
/************************************************************/
/*** Use late transition state algorithms *******************/
/************************************************************/

                   if (have_transfer)
                     {
/*** get vector from transfered atom change of position *****/
                       
                       if (iloop != 1)
                         {
                           last_cofm[0] = cofm_step[0];
                           last_cofm[1] = cofm_step[1];
                           last_cofm[2] = cofm_step[2];
                         }
                       else
                         {
                           latest_bond = -1.0;
                         }

                       if ( !now_close &&
                            ( latest_bond < 0.0 || latest_bond > traj_switch*new_bond_length) )
                         { 
                           cofm_step[0] = start_cofm[0] + step * iloop * transfer_vec[0];
                           cofm_step[1] = start_cofm[1] + step * iloop * transfer_vec[1];
                           cofm_step[2] = start_cofm[2] + step * iloop * transfer_vec[2];
                         }
                       else
                         {
                           if (just_hit)  
                             {
                               just_hit=FALSE;
                               now_close= TRUE;

                               start_cofm[0]=last_cofm[0];
                               start_cofm[1]=last_cofm[1];
                               start_cofm[2]=last_cofm[2];
                               
                               transfer_vec[0] = end_cofm[0] - last_cofm[0];
                               transfer_vec[1] = end_cofm[1] - last_cofm[1];
                               transfer_vec[2] = end_cofm[2] - last_cofm[2];

                               num_steps_left = num_inter+1 - iloop;
                               num_steps_switch = iloop;

                               cstep = 1.0/ (num_steps_left+1);
                             }
                          
                           iii = iloop-num_steps_switch+1;

                           printf("iloop %d num_steps_left %d iii %d\n",
                                      iloop, num_steps_left, iii);

                           cofm_step[0] = start_cofm[0] + cstep * iii * transfer_vec[0];
                           cofm_step[1] = start_cofm[1] + cstep * iii * transfer_vec[1];
                           cofm_step[2] = start_cofm[2] + cstep * iii * transfer_vec[2];
                         }

                     }

                   for (iatom=0; iatom <= num_atoms; iatom++)
                     {
/************************************************************/
/**** Check if the atom is a member of the group defined ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
                       check = FALSE;
                       for (imol=0; imol<=num_chosen_atoms; imol++)
                          {
                             if ( chosen_indices[imol] == iatom ) check = TRUE;
                          }

                       if ( !check )
                         {
                            step_molecule[iatom] = molecule[iatom];
                            step_molecule[iatom].x = step_molecule[iatom].x 
                                          + step * iloop * inter_vec[iatom*3];
                            step_molecule[iatom].y = step_molecule[iatom].y 
                                          + step * iloop * inter_vec[iatom*3+1];
                            step_molecule[iatom].z = step_molecule[iatom].z 
                                          + step * iloop * inter_vec[iatom*3+2];
                         }
                     }
/************************************************************/
/*** Now interpolate the rigid molecule according to late ***/
/************************************************************/

                   if (now_close)
                    {
                      for (imol=0; imol<=num_chosen_atoms; imol++)
                        {
                          iatom = chosen_indices[imol];
                          step_molecule[iatom].x = cofm_step[0] + molecule[iatom].x
                                          + cstep * iii * inter_vec[iatom*3];
                          step_molecule[iatom].y = cofm_step[1] + molecule[iatom].y
                                          + cstep * iii * inter_vec[iatom*3+1];
                          step_molecule[iatom].z = cofm_step[2] + molecule[iatom].z
                                          + cstep * iii * inter_vec[iatom*3+2];
                        }
                    }
                   else
                    {
                      for (imol=0; imol<=num_chosen_atoms; imol++)
                        {
                          iatom = chosen_indices[imol];
                          step_molecule[iatom].x = cofm_step[0] + molecule[iatom].x;
                          step_molecule[iatom].y = cofm_step[1] + molecule[iatom].y;
                          step_molecule[iatom].z = cofm_step[2] + molecule[iatom].z;
                        }
                    }

                   iatom2 = end_molecule[atom_transfered].neighb[0];
 
                   vec[0] = step_molecule[atom_transfered].x
                           -step_molecule[iatom2].x;
 
                   vec[1] = step_molecule[atom_transfered].y
                           -step_molecule[iatom2].y;
 
                   vec[2] = step_molecule[atom_transfered].z
                           -step_molecule[iatom2].z;

                   latest_bond = sqrt( vec[0] * vec[0]
                                      +vec[1] * vec[1]
                                      +vec[2] * vec[2] );

                   printf("latest bond length %10.6f\n", latest_bond);
                 }
/************************************************************/
/*** Otherwise this is a straight linear interpolation ******/
/************************************************************/
               else
                 {
                   for (iatom=0; iatom <= num_atoms; iatom++)
                     {
                       step_molecule[iatom] = molecule[iatom];
                       step_molecule[iatom].x = step_molecule[iatom].x 
                                       + step * iloop * inter_vec[iatom*3];
                       step_molecule[iatom].y = step_molecule[iatom].y 
                                       + step * iloop * inter_vec[iatom*3+1];
                       step_molecule[iatom].z = step_molecule[iatom].z 
                                       + step * iloop * inter_vec[iatom*3+2];
                     }
                 }

               sprintf(step_output,"%s%d", poscar_output, iloop);

               open_file( &fp_vasp_output, step_output, "w");

               write_poscar( fp_vasp_output, &step_molecule[0],  
                             &fract_coords[0], &types[0], num_types, 
                             &latt_vec[0], &scale_factor, num_atoms,
                             &title_line[0], &c_title_line[0], pbc, is_fract);

               fclose(fp_vasp_output);

/***** write pdb file if requested *************************/

               if (need_pdb)
                 {
                   if (iloop==1)
                          open_file( &fp_pdb_output, pdb_output, "w");

                   write_pdb(fp_pdb_output, &step_molecule[0], 
                             &abc[0], num_atoms);
                      
                   if (iloop==num_inter+1)
                          fclose(fp_pdb_output);
                 }
/***** Strip .car from file to allow insertion of number ***/

               iletter=0;
               strcpy(step_output, car_output);
               while (strncmp(&step_output[iletter],".",1) != 0) iletter++;
               step_output[iletter]= '\0';

               sprintf(step_output,"%s%d.car", step_output, iloop);

               open_file( &fp_car_output, step_output, "w");

               start_frame=TRUE;
               write_car( fp_car_output, &header_line[0], &title_line[0], 
                          &c_title_line[0], &date_line[0], &step_molecule[0], 
                          &mol_number[0], pbc, &abc[0], 
                          num_mol_members[0], scale_factor, start_frame );

               fclose(fp_car_output);
           }
   }

    if (need_poscar)
       {
         open_file( &fp_vasp_output, poscar_output, "w");

         sort_by_elem( &molecule[0], num_mol_members[0], &types[0], num_types);

         write_poscar( fp_vasp_output, &molecule[0],  &fract_coords[0],
                       &types[0], num_types, &latt_vec[0], &scale_factor, num_atoms,
                       &title_line[0], &c_title_line[0], pbc, is_fract);

         fclose(fp_vasp_output);
       }
    if (need_car)
       {
          open_file( &fp_car_output, car_output, "w");

          start_frame = TRUE;
          write_car( fp_car_output, &header_line[0], &title_line[0], 
                     &c_title_line[0], &date_line[0], &molecule[0], 
                     &mol_number[0], pbc, &abc[0], 
                     num_mol_members[0], scale_factor, start_frame);

          fclose(fp_car_output);
       }
    if (need_arc)
      {
          printf("Will write %d modes\n", num_modes);

          for ( imode = 0; imode < num_modes; imode++)
            {
              if (which_mode < 0 || which_mode == imode+1)
                {
              sprintf(frame_output,"%s%d.arc",arc_output,imode+1);

              sprintf(c_title_line,"frame 1 Mode %d frequency %10.6f cm-1",
                                             imode+1, eigenvalues[imode]);

              open_file( &fp_arc_output, frame_output, "w");

              start_frame = TRUE;
              theta = 2.0 * pi / num_steps;
              if (linear) step_inc = amplitude/num_steps;

              step = -step_inc;
              for (iloop = 0; iloop < num_steps; iloop++)
                {
                  if (linear) 
                    {
                       step += step_inc;
                    }
                  else
                    {
                       step = amplitude * sin(iloop * theta); 
                    }

                  for (iatom=0; iatom <= num_atoms; iatom++)
                    {
                      step_molecule[iatom] = molecule[iatom];
                      step_molecule[iatom].x = step_molecule[iatom].x 
                                            + step * eigenvecs[imode].dx[iatom];
                      step_molecule[iatom].y = step_molecule[iatom].y 
                                            + step * eigenvecs[imode].dy[iatom];
                      step_molecule[iatom].z = step_molecule[iatom].z 
                                            + step * eigenvecs[imode].dz[iatom];
                    }

                  write_car( fp_arc_output, &header_line[0], 
                             &title_line[0], &c_title_line[0],
                             &date_line[0], &step_molecule[0], 
                             &mol_number[0], pbc, &abc[0], 
                             num_mol_members[0], scale_factor, start_frame);

                  if (need_poscar)
                    {
                      sprintf(frame_output,"%sMode_%d_step%d",
                                      poscar_output, imode+1, iloop);

                      open_file( &fp_vasp_output, frame_output, "w");

                      sort_by_elem( &step_molecule[0], num_mol_members[0], 
                                    &types[0], num_types);

                      write_poscar( fp_vasp_output, &step_molecule[0],  
                                    &fract_coords[0], &types[0], 
                                    num_types, &latt_vec[0], 
                                    &scale_factor, num_mol_members[0],
                                    &title_line[0], &c_title_line[0], 
                                    pbc, is_fract);

                      fclose(fp_vasp_output);
                    }
/*********************************************/
/** Remove repeats of frequency information **/
/*********************************************/

                  sprintf(c_title_line,"frame %d ", iloop+2); 

                  start_frame = FALSE;
                }
              fclose(fp_arc_output);
            }
         }
       }

/*********************************************/
/*** Give force data *************************/
/*********************************************/
     if (need_force)
       {
          printf("Forces read from OUTCAR file (eV/Angstrom)\n");
          rms_x=0.0;
          rms_y=0.0;
          rms_z=0.0;
          rms_t=0.0;
          for (iatom=0; iatom <= num_atoms; iatom++)
            {
               printf("%4s : %10.6f %10.6f %10.6f\n", molecule[iatom].label,
                                                      forces.dx[iatom],
                                                      forces.dy[iatom],
                                                      forces.dz[iatom]);
              
               rms_x += forces.dx[iatom] * forces.dx[iatom];
               rms_y += forces.dy[iatom] * forces.dy[iatom];
               rms_z += forces.dz[iatom] * forces.dz[iatom];
            }
          rms_x = sqrt( rms_x / (num_atoms+1));
          rms_y = sqrt( rms_y / (num_atoms+1));
          rms_z = sqrt( rms_z / (num_atoms+1));

          rms_t = sqrt( rms_x * rms_x + rms_y * rms_y + rms_z * rms_z );

          printf("\n RMS : %10.6f %10.6f %10.6f\n", rms_x, rms_y, rms_z);
          printf("\n RMS magnitude : %10.6f\n", rms_t);
       }
   }
 else
   {
     printf("Bad read!!! returned value: %d\n",good_read);
   }

return 0;
}
