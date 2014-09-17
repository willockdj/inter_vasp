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
#include "reader.h"

#define MAIN 0
#include "data.h"
#include "header.h"
#include "debug.h"
#include "own_maths.h"
#undef MAIN

void open_file(FILE **p_file, char *p_filename, char *p_status);

int read_input( FILE *input_fp, char *p_title, char *p_master_input,
                int *p_is_gulp, int *p_is_car, int *p_is_vasp, int *p_is_siesta,
                char *p_variable_label, char *p_end_input, 
                int *p_end_is_gulp, int *p_end_is_car, int *p_end_is_vasp, 
                double *p_temperature, double *p_min_weight, int *p_assess, 
                int *p_analyse, int *p_num_to_set, int *p_num_per_formula, 
                char *p_potcar_input, char *p_outcar_input, int *p_have_out,
                char *p_incar_input, int *p_have_incar,
                int *p_need_car, char *p_car_output, int *p_need_poscar, 
                int *p_need_arc, char *p_arc_output,
                char *p_poscar_output, int *p_need_freq, int *p_need_force,
                int *p_need_energy, char *p_energy_output,
                int *p_num_inter, char *p_grp_cnt_lab, char *p_grp_cnt_lab2,
                int *p_have_grp, int *p_num_grps, group_lists *p_groups,
                char *p_mol_cnt_lab, int *p_have_mol,
                int *p_num_steps, double *p_amplitude, int *p_linear, 
                int *p_mode, int *p_need_pdb, char *p_pdb_output,
                int *p_need_gulp, char *p_gulp_output,
                int *p_need_late, int *p_need_shells, int *p_num_shell_species,
                labels *p_shell_species, double *p_switch, 
                int *p_need_morph, int *p_need_angle, char *p_axis1_lab,
                char *p_axis2_lab, double *p_mag_after_switch, int *p_super,
                int *p_need_shift, char_list *p_image_files, int *p_num_images,
                char_list *p_outcar_files, int *p_num_outcars,
                char *p_doscar_input, int *p_need_dos, char *p_dos_output, 
                double *p_dos_smear, int *p_need_partdos, int *p_part_dos_list,
                int *p_num_atoms_pdos, int *p_spd, int *p_need_MDtraj,
                char *p_mdtraj_input, int *p_need_multi_dos,  char_list *p_dos_files,
                double *p_dos_weights, int *p_num_dos_files, int *p_read_restart,
                char *p_restart_file, int *p_have_miller, int *p_miller,
                int *p_is_siestados, int *p_is_vaspdos, int *p_need_md_run,  
                int *p_compare_modes, char *p_mode_compare_input, int *p_end_min_image);

int read_car( FILE *fp, int *p_header_line, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_abc, int *p_been_before, group_lists *p_groups,
              int have_grp, int *p_num_grps,
              int find_fixed, coord_flags *p_fix_flags);

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
                 double *p_scale_factor, coord_flags *p_fix_flags);

int read_fdf( FILE *fp, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_latt_vec, double *p_recip_latt_vec,  
              double *p_abc, int *p_been_before, int *p_is_fract,
              int *p_is_cart, int *p_ion_number, int *p_num_types, 
              double *p_scale_factor, int *p_num_free, coord_flags *p_fix_flags);

int read_outcar( FILE *fp, atom *p_molecule, double *p_latt_vec,
                 double *p_recip_latt_vec, double *p_abc, double *p_eigenvals,
                 e_vec *p_eigenvecs, e_vec *p_forces, e_vec *p_chain_forces,
                 int *p_num_atoms, atom_number *p_types, int *p_num_types,
                 int *p_num_modes, int need_freq, int need_force,
                 int need_energy, int *p_have_band, double *p_energy,
                 int need_fermi, double *p_fermi );

int read_incar( FILE *fp, double *p_magmom, int *p_num_magmom );

int read_potcar( FILE *fp, atom *p_molecule,
                 int *p_ion_number, int *p_num_types);

int read_xdatcar(FILE *fp, int *p_num_frames, atom *p_molecule, int num_atoms,
                 coord_flags *p_fix_flags );

void generate_neighbours( atom *p_molecule, int num_atoms, 
                          atom_number *p_types, int *p_num_types,
                          int use_pbc,  double *p_recip_latt_vec, double *p_latt_vec,
                          charge_list *p_spec_charges);

void sort_by_elem( atom *p_molecule, int num_atoms, atom_number *p_types, 
                   int num_types );

void cart_latt_vecs( double *p_abc, double *p_latt_vec, 
                                         double *p_recip_latt_vec);

int compare_strings( char *p_ichar1, char *p_ichar2 );

void string_to_int(char *p_ichar2, int *p_ichar1, int max_position );

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, 
                double *p_recip_latt_vec, coord_flags *p_fix_flags,
                double *p_magmom, int num_magmom );

void write_gulp(FILE *fp, int *p_title_line, atom *p_molecule, int num_atoms,
                atom *p_shells, int num_shells, int fract_or_cart, double *p_abc,
                int space_group, double *p_recip_latt_vec, double *p_latt_vec,
                charge_list *p_spec_charges, int num_species, int *p_super, 
                int need_shells, labels *p_shell_species, int num_shell_species);

void write_poscar( FILE *fp, atom *p_molecule, double *p_fract_coords,
                   atom_number *p_types, int num_types,
                   double *p_latt_vec, double *p_scale_factor, int num_atoms,
                   int *p_title_line, char *p_c_title_line, int pbc, 
                   int is_fract, coord_flags *p_fix_flags);

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
                        char *p_mol_cnt_lab,
                        double *p_start_cofm, 
                        double *p_end_cofm, double *p_inter_cofm);

void react_coords(atom *p_molecule, atom *p_end_molecule, int num_atoms,
                  char_list *p_image_files, int num_images,
                  double *p_recip_latt_vec, double *p_latt_vec,
                  int pbc);

void write_pdb(FILE *pdb_fp, atom *p_molecule, double *p_abc, int num_atoms,
               int *p_super, double *p_latt_vec );

void min_image( double *x, double *y, double *z,
                           double *p_recip_latt_vec, double *p_latt_vec);

double atom_separation_squared(atom *p_A, atom *p_B, int pbc,
                           double *p_recip_latt_vec, double *p_latt_vec);

/******************************************************************/
/******************************************************************/
/******************************************************************/

int read_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                 int num_atoms_pdos, dos *p_pdos, int *p_spd, int num_columns );

int count_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                 int num_atoms_pdos, dos *p_pdos, int *p_spd );

void smear_dos( dos *p_dos, int num_dos, double smear );

void d_centre( dos *p_dos, int num_dos, double fermi, double *p_dcentre, double *p_dfilling );

void write_doscsv(FILE *fp, char *p_title_x, char *p_title_y, 
                  char *p_title_z,
                  dos *p_dos,
                  int have_tot, int num, int num_columns, double fermi);

void cut_to_dot( char *p_word );

int read_siesta_vectors( FILE *fp, double *p_eigenvals, 
                         e_vec *p_eigenvecs, int *p_num_atoms, int num_free,
                         int *p_num_modes);

void reorientate_cell(atom *p_molecule, int *p_num_atoms, double *p_latt, 
                      double *p_recip_latt, double *p_abc, int *p_miller,
                      double *p_new_latt, double *p_new_recip_latt, 
                      double *p_new_abc, atom *p_slab_mol, int *p_num_slab_atoms);

/******************************************************************/
/*** Start of MAIN ************************************************/
/******************************************************************/

int main(argc, argv)
  int   argc;
  char  **argv;
{

  atom molecule[MAXATOMS];
  atom end_molecule[MAXATOMS];
  atom slab_mol[MAXATOMS];
  atom step_molecule[MAXATOMS];
  atom chosen_one[MAXATOMS];
  atom chosen_one_end[MAXATOMS];
  atom shells[MAXATOMS];
  atom end_shells[MAXATOMS];

  atom_number types[MAXATOMS], end_types[MAXATOMS];

  bonds old_bonds[MAXBONDS], new_bonds[MAXBONDS];
  int num_old_bonds, num_new_bonds;

  e_vec eigenvecs[MAXMODES], forces, chain_forces;
  e_vec new_eigenvecs[MAXMODES];

  char_list image_files[MAXIMAGES];
  char_list outcar_files[MAXOUTCARS];
  char_list dos_files[MAX_DOS_FILES];

  double eigenvalues[MAXMODES];
  double magmom[MAXATOMS];
  double new_eigenvalues[MAXMODES];
  double dist, dmin, diffx, diffy, diffz, max_freq_shift;
  double rms_x, rms_y, rms_z, rms_t;
  double energy_vasp;
  double fermi;

  int num_atoms, num_shells, good_read, num_in_group, num_groups, num_modes;
  int num_new_atoms, num_new_types, num_new_modes;
  int num_chosen_atoms, num_free, match, isign;
  int end_num_atoms, end_num_shells, start_frame, imode, jmode;
  int done_start_frame = FALSE;
  int header_line[LINESIZ], title_line[LINESIZ], date_line[LINESIZ];
  int end_header_line[LINESIZ], end_title_line[LINESIZ], end_date_line[LINESIZ];
  int pbc, num_of_mols, num_mol_members[MAXMOL], mol_number[MAXATOMS]; 
  int end_num_of_mols, end_num_mol_members[MAXMOL], end_mol_number[MAXATOMS]; 
  int map_atoms[MAXATOMS];
  int been_before, grp_mem_ind, grp_cnt_ind, grp_cnt_ind2, igrp, jgrp, imol;
  int is_centre, grp_cnt;
  int have_band, have_mol, mol_cnt_ind, mol_ind;
  int iii,iloop, jloop, iatom, jatom, iatom2, icomp, icoord, iletter;
  int idos, index;
  int ineigh, neigh_index, num_types, end_num_types;
  int ibond, ineigh2, still_there;
  int start_mol, num_in_this_mol;
  int ion_number[MAXTYPES];
  int num_steps, linear;
  int just_hit, now_close, num_steps_left;
  int num_steps_switch;
  int num_steps_after_switch;
  int iframe, num_frames;
  int flags[MAXATOMS];
  int num_images, num_outcars;
  int have_tot, end_min_image;
  int is_siestados, is_vaspdos;
  int compare_modes;
  
  coord_flags fix_flags[MAXATOMS];

  double abc[6], new_abc[6], theta, dot, tot_dot, dot_list[MAXMODES];
  double recip_latt_vec[9], latt_vec[9];
  double new_recip_latt_vec[9], new_latt_vec[9];
  double vec[3], vec1[3], vec2[3], size;
  double transfer_vec[3];
  double unit[3], move_dist;
  double delta_bond1[MAXATOMS], new_bond;
  double delta_bond2[MAXATOMS];
  double orig_bond1[MAXATOMS];
  double orig_bond2[MAXATOMS];
  double this_delta, this_orig;
  double amplitude, new_bond_length;
  double latest_bond, temp_bond;
  double d, dx, dy, dz;
  double fa, fb, fc;
  double timestep;
  double start_cofm[3];
  double end_cofm[3];
  double sum_pos[3], sum_pos1[3], sum_check[3];
  double shift_centre[3];
  double inter_cofm[3];
  double last_cofm[3];
  double cofm_step[3];
  double total_mass;
  double traj_switch;
  double distance_gone, distance_left;
  double set_temp, ke, speed;
  double dcentre, dfilling;

  char *p_key;
  char c_title_line[LINESIZ];
  char c_header_line[LINESIZ];
  char title_x[LINESIZ];
  char title_y[LINESIZ];
  char title_z[LINESIZ];
  char grp_cnt_lab[7];
  char grp_cnt_lab2[7];
  char mol_cnt_lab[7];
  char axis1_lab[7];
  char axis2_lab[7];
  char elem_rec[7];
  char variable_label[60];
  char master_input[60];
  char end_input[60];
  char potcar_input[60];
  char incar_input[60];
  char restart_file[60];
  char mdtraj_input[60];
  char outcar_input[60];
  char doscar_input[60];
  char filename[60];
  char dos_output[60];
  char poscar_output[60];
  char step_output[60];
  char car_output[60];
  char arc_output[60];
  char pdb_output[60];
  char gulp_output[60];
  char frame_output[60];
  char frame_output_pdb[60];
  char mode_compare_input[60];
  char command[60];
  char energy_output[60];

/******************************************************************/
/*** Variables specific to the job in hand ************************/
/******************************************************************/

  FILE *fp_input_frame;
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
  FILE *fp_siesta_input;
  FILE *fp_vasp_output;
  FILE *fp_doscar_input;
  FILE *fp_dos_output;
  FILE *fp_input;

double fract_coords[3*MAXATOMS], end_fract_coords[3*MAXATOMS];
double inter_vec[3*MAXATOMS];
double temperature,  min_weight;
double scale_factor, step, cstep, step_inc;
double mag_after_switch;
double axis[3], origin[3];
double dos_smear, dos_weights[MAX_DOS_FILES];
double dos_norm;

int is_fract, is_cart, is_gulp, is_car, is_vasp, assess, analyse;
int is_siesta;
int end_is_gulp, end_is_car, end_is_vasp;
int num_to_set, num_per_formula;
int axis1_ind, axis2_ind;
int top_bit[500], bottom_bit[500];
int num_top_chars, num_bottom_chars, num_species;
int space_group, state;
int need_car, need_arc, need_poscar, need_interpolate, need_freq;
int need_mdtraj, need_multi_dos, num_dos_files;
int need_shift, need_react_coord;
int need_pdb, need_gulp, need_late, need_morph, need_shells;
int need_angle, need_dos, num_dos_points, spd[4];
int need_partdos, part_dos_list[MAXATOMS], num_atoms_pdos;
int need_md_run;
int lowpdos, hipdos, this_pdos, nrang;
int num_shell_species;
int have_transfer, atom_transfered;
int need_force, which_mode;
int need_energy;
int need_fermi;
int num_inter, have_grp, num_grps, have_out;
int have_incar;
int num, num_magmom;
int check, check1, check2;
int have_miller, miller[3], look_for_layers;
int chosen_indices[MAXATOMS];
int super[3], surf_super[3];
int any_rec, matched_neigh;
int read_restart;
int num_columns;
int num_slab_atoms;

group_lists groups;

charge_list spec_charges[MAXATOMS];

labels shell_species[10];

dos ndos[MAX_DOS], pdos[MAX_DOS];
dos temp_ndos[MAX_DOS], temp_pdos[MAX_DOS];

/*******************************************************************/
/*********** read input file (name input by user on the ************/
/*********** command line).                             ************/
/*******************************************************************/

printf("Starting inter_vasp version 2.1.0\n");
printf("Last Built : 10th January 2007\n");

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
strcpy(mdtraj_input, "No XDATCAR file Supplied");
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
dos_smear = 0.2;
num_magmom=-1;

pi = acos(-1.0);
printf("Pi is %10.6f\n", pi);

end_is_gulp= FALSE;
end_is_vasp= FALSE;
end_is_car = FALSE;
is_gulp= FALSE;
is_vasp= FALSE;
is_car = FALSE;
is_siesta = FALSE;
need_car = FALSE;
need_shift = FALSE;
need_pdb = FALSE;
need_gulp = FALSE;
need_shells = FALSE;
num_shell_species=-1;
need_arc = FALSE;
need_poscar = FALSE;
need_mdtraj = FALSE;
need_dos = FALSE;
need_multi_dos = FALSE;
spd[0]= TRUE;
spd[1]= TRUE;
spd[2]= TRUE;
need_partdos = FALSE;
need_freq = FALSE;
compare_modes = FALSE;
need_late = FALSE;
any_rec = FALSE;
need_angle = FALSE;
need_force = FALSE;
have_band = FALSE;
have_miller = FALSE;
need_morph = FALSE;
num_to_set = 0;
num_per_formula = 0;
traj_switch = 1.5;
mag_after_switch= 0.5;
num_shells=-1;
num_dos_files=0;
num_images=-1;

super[0] = 1;
super[1] = 1;
super[2] = 1;

surf_super[0] = 1;
surf_super[1] = 1;
surf_super[2] = 1;

temperature = 0.0;
min_weight  = 0.0;

num_grps=-1;

read_input( fp_control_file, &c_title_line[0], &master_input[0],
            &is_gulp, &is_car, &is_vasp, &is_siesta, &variable_label[0],
            &end_input[0], &end_is_gulp, &end_is_car, &end_is_vasp, 
            &temperature, &min_weight, &assess, 
            &analyse, &num_to_set, &num_per_formula, 
            &potcar_input[0], &outcar_input[0], &have_out, 
            &incar_input[0], &have_incar,
            &need_car, &car_output[0], 
	    &need_poscar, &need_arc, &arc_output[0],
            &poscar_output[0], &need_freq, &need_force,
            &need_energy, &energy_output[0],
            &num_inter, &grp_cnt_lab[0], &grp_cnt_lab2[0],
            &have_grp, &num_grps, &groups,
            &mol_cnt_lab[0], &have_mol,
            &num_steps, &amplitude, &linear, &which_mode,
            &need_pdb, &pdb_output[0], &need_gulp, &gulp_output[0],
            &need_late, &need_shells, &num_shell_species, &shell_species[0],
            &traj_switch, &need_morph, &need_angle,
            &axis1_lab[0], &axis2_lab[0], &mag_after_switch, &super[0],
            &need_shift, &image_files[0], &num_images,
            &outcar_files[0], &num_outcars,
            &doscar_input[0], &need_dos, &dos_output[0], 
            &dos_smear, &need_partdos, &part_dos_list[0], &num_atoms_pdos,
            &spd[0], &need_mdtraj, &mdtraj_input[0], &need_multi_dos,
            &dos_files[0], &dos_weights[0], &num_dos_files,  
            &read_restart, &restart_file[0], &have_miller, &miller[0],
            &is_siestados, &is_vaspdos, &need_md_run,
            &compare_modes, &mode_compare_input[0], &end_min_image );

if (need_md_run)
  {
     printf("Will perform our own MD run using VASP forces\n");
  }

if (need_multi_dos)
  {
    printf("Found %d DOSCAR files will expect to amalgamate them\n", num_dos_files);

    for (iloop=0; iloop <= num_dos_files; iloop++)
      {
         printf("%s with weight %10.6f\n", dos_files[iloop].name, dos_weights[iloop]);
      }
  }

need_interpolate = (end_is_gulp || end_is_car || end_is_vasp) && num_images < 0;
need_react_coord = (end_is_gulp || end_is_car || end_is_vasp) && num_images > 0;

if (num_images > 0 && !need_react_coord)
  {
    printf("ERROR: images supplied but no end point reference structure\n");
    exit(0);
  }
else if (num_images==0)
  {
    printf("ERROR: Cannot have zero images in a react coords run\n");
    exit(0);
  }

fclose(fp_control_file);

printf("Job Title        : %s\n\n", c_title_line);
printf("Master file      : %s\n", master_input);
if (read_restart)
   {
      printf("Will take atom co-ordinates from the restart (CONTCAR) file: %s\n", restart_file);
   }
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
if (is_siesta)
   {
     printf("                           this is a SIESTA structure file\n");
     printf(" The fdf file is : %s\n", master_input);
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
if (end_min_image)
   {
     printf("Will use minimum image of end point atoms to corresponding start point for interpolation\n");
   }
else
   {
     printf("Will use end point co-ordinates as supplied\n");
   }


if ( need_interpolate ) 
  {
    printf("Will provide %d interpolated structures\n", num_inter);
    if (need_shift)
      {
         printf("Will shift all structures to elliminate drift during interpolation\n");
      }
    if (need_arc)
      {
        printf("Will generate an arc file of the elastic band set of structures named : %s\n",arc_output);

        open_file( &fp_arc_output, arc_output, "w");
      }

    if (have_grp) 
      {
             
             printf("Have %d groups :\n", num_grps);

                  if (groups.group_type1 == CENTRE_TYPE)
                    {
                       printf("Will interpolate GRP1 group centred on >>%s<<\n", grp_cnt_lab);
                    }
                  else if (groups.group_type1 == ANGLE_TYPE)
                    {
                       printf("Will rotate the GRP1 atoms around the axis defined by the line\n");
                       printf("from atom >>%s<< to atom >>%s<<\n", axis1_lab, axis2_lab);
                    }
                  else
                    {
                       printf("ERROR: No type set for group interpolation\n");
                       exit(0);
                    }

             if (num_grps > 1)
               {
                  if (groups.group_type2 == CENTRE_TYPE)
                    {
                       printf("Will interpolate GRP2 group centred on >>%s<<\n", grp_cnt_lab2);
                    }
                  else if (groups.group_type2 == ANGLE_TYPE)
                    {
                       printf("Will rotate the GRP2 atoms around the axis defined by the line\n");
                       printf("from atom >>%s<< to atom >>%s<<\n", axis1_lab, axis2_lab);
                    }
               }

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

        if (have_grp)
          {
            printf("ERROR: Request for late and group method interpolation, decide on one or the other!\n");
            exit(0);
          }
        printf("Will interpolate molecule containing");
        printf(" >>%s<< to look for late transition state.\n",mol_cnt_lab);
        printf("Using %10.6f as the multiplier of the final bondlength at which the", 
                                                                          traj_switch);
         
        printf("After switch will do increase image density by %10.6f\n", mag_after_switch);
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
else
  {
    printf("This is not an interpolation run\n");
  }

/****** Try and drive VASP MD using our own integration routines ****/

if (need_md_run)
  {
     timestep = 0.001;
     set_temp = 300;
     open_file( &fp_outcar, master_input, "r");

/**** Read in current forces and energy from OUTCAR ****/

     need_fermi=FALSE;
     good_read = read_outcar( fp_outcar, &molecule[0], &latt_vec[0],
                              &recip_latt_vec[0], &abc[0], &eigenvalues[0],
                              &eigenvecs[0], &forces, &chain_forces,
                              &num_atoms, &types[0], &num_types, &num_modes, need_freq,
                              need_force, need_energy, &have_band, &energy_vasp, need_fermi,
                              &fermi);

     pbc = TRUE;   /* This should be set by read_outcar */

     printf("Back from read_outcar with %d atoms\n\n", num_atoms);

     for (iatom= 0; iatom <= num_atoms; iatom++)
       {
         printf("%d %s %10.6f  %10.6f  %10.6f  %10.6f force: %10.6f  %10.6f  %10.6f\n",
                                 iatom+1, molecule[iatom].label,
                                 molecule[iatom].mass,
                                 molecule[iatom].x, 
                                 molecule[iatom].y, 
                                 molecule[iatom].z, 
                                 forces.dx[iatom],
                                 forces.dy[iatom],
                                 forces.dz[iatom]);
       }

/**** Write arc file ****/

     open_file( &fp_arc_output, "md_movie.arc", "w");

     sprintf(c_title_line,
             "Inter_vasp generated arc file.");

     scale_factor = 1.0;

     start_frame=TRUE;

     write_car( fp_arc_output, &header_line[0],
                &title_line[0], &c_title_line[0],
                &date_line[0], &molecule[0],
                &mol_number[0], pbc, &abc[0],
                num_atoms+1, scale_factor, start_frame,
                &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0],
                &magmom[0], num_magmom);

     start_frame=FALSE;

/**** Set all velocities to the mean for temperature ( in 2D ) ****/
/**** using KE = RT and 1/2 m v^2                              ****/

       ke = BOLTZ * set_temp;

     for (iatom= 0; iatom <= num_atoms; iatom++)
       {
       speed = sqrt( 2.0 * ke / molecule[iatom].mass);

       if ( iatom == 0 )
               printf("Initial speed set to %10.6f\n", speed);

/****  ( 1.0 - 2.0*drand48() ) will be distributed +/- 1 **********/

       vec[0] =  1.0 - 2.0*drand48();
       vec[1] =  1.0 - 2.0*drand48();
       vec[2] =  1.0 - 2.0*drand48();

       unit_vector(&vec[0]);

       molecule[iatom].vx = vec[0] * speed;
       molecule[iatom].vy = vec[1] * speed;
       molecule[iatom].vz = vec[2] * speed;

     }

/**** Loop over number of MD runs ****/

    for (iatom=0; iatom <= num_atoms; iatom++)
       {

/*** a(t) = f(t) / m *****/
 
/*** Include conversion between eV and internal units for forces ***/

          molecule[iatom].ax = forces.dx[iatom]/(INTUNITS_TO_EV*molecule[iatom].mass);
          molecule[iatom].ay = forces.dy[iatom]/(INTUNITS_TO_EV*molecule[iatom].mass);
          molecule[iatom].az = forces.dz[iatom]/(INTUNITS_TO_EV*molecule[iatom].mass);

/*** v(t+dt/2) = v(t-dt/2) + a(t) ***/
/*** for first step put in half iteration to offset velocity from positions ****/
          if (iloop==0)
            {
              molecule[iatom].vx += molecule[iatom].ax * 0.5*timestep;
              molecule[iatom].vy += molecule[iatom].ay * 0.5*timestep;
              molecule[iatom].vz += molecule[iatom].az * 0.5*timestep;
            }
          else
            {
              molecule[iatom].vx += molecule[iatom].ax * timestep;
              molecule[iatom].vy += molecule[iatom].ay * timestep;
              molecule[iatom].vz += molecule[iatom].az * timestep;

            }

/*** r(t+dt) = r(t) + v(t+dt/2) t *******/
          molecule[iatom].x += molecule[iatom].vx * timestep;
          molecule[iatom].y += molecule[iatom].vy * timestep;
          molecule[iatom].z += molecule[iatom].vz * timestep;


       }
     printf("\n");

     for (iatom= 0; iatom <= num_atoms; iatom++)
       {
         printf("%d %s %10.6f  %10.6f  %10.6f  %10.6f acc: %10.6f  %10.6f  %10.6f\n",
                                 iatom+1, molecule[iatom].label,
                                 molecule[iatom].mass,
                                 molecule[iatom].x,
                                 molecule[iatom].y,
                                 molecule[iatom].z,
                                 molecule[iatom].ax,
                                 molecule[iatom].ay,
                                 molecule[iatom].az);
       }


/**** Write the new co-ordinates to a POSCAR file ******/

    open_file( &fp_vasp_output, "POSCAR", "w");

    scale_factor = 1.0;
    write_poscar( fp_vasp_output, &molecule[0],  &fract_coords[0],
                  &types[0], num_types, &latt_vec[0], &scale_factor, num_atoms,
                  &title_line[0], &c_title_line[0], pbc, FALSE, &fix_flags[0]);

    fclose(fp_vasp_output);

/**** Run the vasp job to get next set of forces ***/

    fflush(stdout);

    sprintf(command, "./vasp_script");

    system(command);

/**** Write arc file ****/

    write_car( fp_arc_output, &header_line[0],
               &title_line[0], &c_title_line[0],
               &date_line[0], &molecule[0],
               &mol_number[0], pbc, &abc[0],
               num_atoms+1, scale_factor, start_frame,
               &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
               &magmom[0], num_magmom);

/**** End loop over number of MD runs here ****/

    fclose(fp_arc_output); 
    exit(0);
  }

/****** Report on DOS expectations *******/

if (need_dos)
  {
     if (is_vaspdos)
      {
        printf("Will write density of states files based on DOSCAR file: %s\n", doscar_input);
      }
     else if (is_siestados)
      {
        printf("Will write density of states files based on SIESTA .DOS file: %s\n", doscar_input);
        exit(0);
      }
     else
      {
        printf("ERROR: need_dos requested but no DOS data supplied\n");
        exit(0);
      }
     printf("Output csv format files for DOS will be given stem     : %s_ndos\n", dos_output);
     printf("DOS smearing to use                                    : %10.6f\n", dos_smear);

     if (need_partdos)
       {
         printf("Output csv format files for PDOS will be given stem: %s_pdos\n", dos_output);
         printf("Atoms included in PDOS: ");


/**** Sort part_dos_list into order ****/         

         for (iloop=0; iloop<=num_atoms_pdos; iloop++)
           {
             for (jloop=iloop+1; jloop<=num_atoms_pdos; jloop++)
               {
                 if ( part_dos_list[iloop] > part_dos_list[jloop] )
                    {
                       iii = part_dos_list[jloop];
                       part_dos_list[jloop]= part_dos_list[iloop];
                       part_dos_list[iloop]=iii;
                    } 
               } 
           } 

         for (iloop=0; iloop<=num_atoms_pdos; iloop++)
           {
              printf("%d ", part_dos_list[iloop]);
           }
         printf("\n");

         if (!spd[0] && !spd[1] && !spd[2] )  
          {
             printf("ERROR: Additional flag given for PDOS orbitals does not include s p or d\n");
             exit(0);
          }
        printf("Will include ");
        if (spd[0]) printf("s ");
        if (spd[1]) printf("p ");
        if (spd[2]) printf("d ");
        printf("orbitals in PDOS contributions\n");
      }

     if (need_multi_dos)
       {
          printf("Multi-DOS to sum\n");
          for (idos=0; idos <= num_dos_points; idos++)
            {
               ndos[idos].up_dos = 0.0;
               ndos[idos].down_dos = 0.0;
               ndos[idos].up_totdos = 0.0;
               ndos[idos].down_totdos = 0.0;
            }

          if (need_partdos)
            {
              for (idos=0; idos <= num_dos_points; idos++)
                {
                   pdos[idos].up_dos = 0.0;
                   pdos[idos].down_dos = 0.0;
                   pdos[idos].up_totdos = 0.0;
                   pdos[idos].down_totdos = 0.0;
                }
            }

          dos_norm=0.0;
          for (iii=0; iii <= num_dos_files; iii++)
            {
              open_file( &fp_doscar_input, dos_files[iii].name, "r");

               num_dos_points= read_doscar( fp_doscar_input, &temp_ndos[0], need_partdos, 
                                            &part_dos_list[0], num_atoms_pdos, &temp_pdos[0],
                                            &spd[0], num_columns );

               if (iii==0)
                 {
                   for (idos=0; idos < num_dos_points; idos++)
                     {
                        ndos[idos].energy=temp_ndos[idos].energy;
                     }
                 }
               else
                 {
/***** check DOSCAR files were aligned on energy axis ****/
                   for (idos=0; idos < num_dos_points; idos++)
                     {
                        if ( fabs(ndos[idos].energy-temp_ndos[idos].energy) > 0.0001)
                          {
                             printf("ERROR: Energy axes in DOSCAR files not aligned cannot combined them.\n");
                             printf("       Use EMIN and EMAX in INCAR file to define the range you want\n");
                             printf("       use the same values in all single k-point calculations.\n");
                             exit(0);
                          }
                     }
                 }

               for (idos=0; idos < num_dos_points; idos++)
                 {
                    ndos[idos].up_dos += dos_weights[iii]*temp_ndos[idos].up_dos;
                    ndos[idos].down_dos += dos_weights[iii]*temp_ndos[idos].down_dos;
                    ndos[idos].up_totdos += dos_weights[iii]*temp_ndos[idos].up_totdos;
                    ndos[idos].down_totdos += dos_weights[iii]*temp_ndos[idos].down_totdos;

                    dos_norm+= dos_weights[iii]; 
                 }

/* Sum part_dos similarly */
  
               if (need_partdos)
                 {
                    for (idos=0; idos < num_dos_points; idos++)
                      {
                         pdos[idos].up_dos += dos_weights[iii]*temp_pdos[idos].up_dos;
                         pdos[idos].down_dos += dos_weights[iii]*temp_pdos[idos].down_dos;
                         pdos[idos].up_totdos += dos_weights[iii]*temp_pdos[idos].up_totdos;
                         pdos[idos].down_totdos += dos_weights[iii]*temp_pdos[idos].down_totdos;
                      }
                 }

               printf("Found %d dos points in file %d: %s\n", num_dos_points, iii,
                                                              dos_files[iii].name );
               fclose(fp_doscar_input);
           }
         for (idos=0; idos < num_dos_points; idos++)
           {
              ndos[idos].up_dos = ndos[idos].up_dos / dos_norm;
              ndos[idos].down_dos = ndos[idos].down_dos  / dos_norm;
              ndos[idos].up_totdos = ndos[idos].up_totdos / dos_norm;
              ndos[idos].down_totdos = ndos[idos].down_totdos / dos_norm;
           }
  
         if (need_partdos)
           {
              for (idos=0; idos < num_dos_points; idos++)
                {
                   pdos[idos].up_dos = pdos[idos].up_dos / dos_norm;
                   pdos[idos].down_dos = pdos[idos].down_dos  / dos_norm;
                   pdos[idos].up_totdos = pdos[idos].up_totdos / dos_norm;
                   pdos[idos].down_totdos = pdos[idos].down_totdos / dos_norm;
                }
           }
       }
     else
       {
          printf("Reading single DOSCAR file.\n");

          open_file( &fp_doscar_input, doscar_input, "r");
               
          num_columns= count_doscar( fp_doscar_input, &temp_ndos[0], need_partdos, 
                                      &part_dos_list[0], num_atoms_pdos, &temp_pdos[0],
                                      &spd[0] );

          printf("main: number of orbitals for pdos = %d\n", num_columns);
          fclose(fp_doscar_input);
         
          open_file( &fp_doscar_input, doscar_input, "r"); 

          num_dos_points= read_doscar( fp_doscar_input, &ndos[0], need_partdos, 
                                       &part_dos_list[0], num_atoms_pdos, &pdos[0],
                                       &spd[0], num_columns );

          printf("Found %d dos points\n", num_dos_points);
          fclose(fp_doscar_input);
      }

     smear_dos( &ndos[0], num_dos_points, dos_smear );

     sprintf(filename, "%s_ndos.csv",dos_output);
     open_file( &fp_dos_output, filename, "w");

     sprintf(title_x, "Energy (eV)");
     sprintf(title_y, "DOS");
     sprintf(title_z, "TDOS");
     have_tot = TRUE;

     open_file( &fp_outcar, master_input, "r");

     good_read = read_outcar( fp_outcar, &molecule[0], &latt_vec[0],
                   &recip_latt_vec[0], &abc[0], &eigenvalues[0],
                   &eigenvecs[0], &forces, &chain_forces,
                   &num_atoms, &types[0], &num_types, &num_modes, need_freq,
                   need_force, need_energy, &have_band, &energy_vasp,
                   TRUE, &fermi);

     fclose(fp_outcar);

     write_doscsv( fp_dos_output, title_x, title_y, 
                  title_z, &ndos[0], have_tot, num_dos_points-1, num_columns, fermi);

     fclose(fp_dos_output);

     if (need_partdos)
       {
         smear_dos( &pdos[0], num_dos_points, dos_smear );

         sprintf(filename, "%s_pdos.csv",dos_output);
         open_file( &fp_dos_output, filename, "w");

         sprintf(title_x, "Energy (eV)");
         sprintf(title_y, "PDOS");

/**** Sort out title_y to include atoms in pdos including ****/
/**** possibility of a range of atoms being present       ****/

         nrang=-1;
         lowpdos=-1;
         for (iloop=0; iloop<=num_atoms_pdos; iloop++)
           {
               printf("part_dos_list is %d \n", part_dos_list[iloop]);
               if (lowpdos<0) 
                 {
                   lowpdos=part_dos_list[iloop];
                   printf("lowpdos set to %d \n", lowpdos);
                   nrang++; 
                 }
               this_pdos=part_dos_list[iloop];

               if ( this_pdos == lowpdos + nrang )
                 {
                   hipdos= this_pdos;
                   printf("lowpdos %d and hipdos %d define range\n", lowpdos, hipdos);
                   nrang++; 
                 }
               else if (nrang==1)
                 {
                   printf("lowpdos %d Single number\n", lowpdos);
                   sprintf(title_y,"%s_%d", title_y, lowpdos);
                   lowpdos=this_pdos;
                   nrang=1;
                 }
               else 
                 {
                   printf("Writing range %d to %d\n",lowpdos, hipdos);
                   sprintf(title_y,"%s_%d-%d", title_y, lowpdos, hipdos);
                   lowpdos=this_pdos;
                   nrang=1;
                 }

           }

         printf("Left loop with nrang %d lowpdos %d this_pdos %d and hipdos %d\n",
                          nrang, lowpdos, this_pdos, hipdos);
         if (nrang > 1)
           {
               printf("Range to write left over\n");
               sprintf(title_y,"%s_%d-%d", title_y, lowpdos, hipdos);
           }
         else
           {
               printf("One to write left over\n");
               sprintf(title_y,"%s_%d", title_y, this_pdos);
           }

         sprintf(title_y, "%s_", title_y);
         if (spd[0]) sprintf(title_y, "%ss", title_y); 
         if (spd[1]) sprintf(title_y, "%sp", title_y); 
         if (spd[2]) sprintf(title_y, "%sd", title_y); 
         
         sprintf(title_z, "TPDOS");
         have_tot = TRUE;

         write_doscsv( fp_dos_output, title_x, title_y, 
                      title_z, &pdos[0], have_tot, num_dos_points-1, num_columns, fermi);

       }
     exit(0);
  }


/****** Catch request for pdos with no dos file name defined ******/

if (need_partdos && !need_dos)
  {
     printf("ERROR: Partial DOS requested without 'need dos_file' directive.\n");
     exit(0);
  }


if (need_react_coord)
  {
    printf("Will provide reaction co-ordinate values for images with respect the master and end points\n");

  }

if (have_out) printf("OUTCAR file is   : %s\n", outcar_input);
if (have_incar) 
   {
     printf("INCAR file is   : %s\n", incar_input);

/*** For now just look for and read MAGMON values ***/

     open_file( &fp_input, incar_input, "r");

     printf("Reading file...\n");

     read_incar( fp_input, &magmom[0], &num_magmom);

     printf("Found %d magmom entries\n", num_magmom+1);

     iloop=0;
     while (iloop < num_magmom)
       {
         if ( num_magmom - iloop > 10 ) num=10;
         else num = num_magmom - iloop +1;
         for (jloop=0; jloop < num; jloop++)
           {
              printf("%3d : %5.2f  ", iloop+1, magmom[iloop]);
              iloop++;
           }
         printf("\n");
       }

     fclose(fp_input);
     
   }

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

/*************************************************************/
/** Report on frequency run data input ***********************/
/*************************************************************/
if (need_freq)
  {
     if (have_out)
       {
         printf("Will analyse OUTCAR file for vibrational frequencies\n");
         if (which_mode > 0)
              printf("Animation files for mode %d only requested\n", 
                            which_mode);
       }
     else if (is_siesta)
       {
         printf(" Will look for Eigenvalues and vectors in siesta '.vectors' file\n");
       }
     else 
       {
         printf("ERROR : Frequencies requested but no OUTCAR supplied\n");
       }
     if (need_arc)
	{
	  printf("Will produce an arc file with general stem >>%s<<\n", arc_output);
          printf("Using %d frames and oscillation amplitude %10.6f\n",
                                   num_steps, amplitude);
          if (linear) 
           {
             printf("frames will be linearly related, for eigen-following\n");
             printf("separation between steps will be at %10.6f intervals along eigenvector\n",
                           amplitude/num_steps);
           }
	}
    }

if (need_car)
	{
	  printf("Will produce a car file called %s\n", car_output);
          printf("Will generate super cell of %d %d %d extent\n",
                            super[0], super[1], super[2]);
	}
if (need_pdb)
	{
	  printf("Will produce a pdb file named >>%s<<\n", pdb_output);
          printf("pdb files will show super cell of %d %d %d extent\n",
                            super[0], super[1], super[2]);
	}
if (need_gulp)
	{
	  printf("Will produce a gulp file named >>%s<<\n", gulp_output);
          printf("gulp files will show super cell of %d %d %d extent\n",
                            super[0], super[1], super[2]);
          if (need_shells)
            {
               printf("Will include shells for species: ");
               for (iloop=0; iloop<=num_shell_species; iloop++)
                 {
                   printf("%s ", shell_species[iloop].label);
                 }
               printf("\n");
            }
	}
if (need_poscar)
	{
	  printf("Will produce a POSCAR file called %s\n", poscar_output);
	}
if (need_mdtraj)
        {
          printf("Will convert the trajectory file %s to format ", mdtraj_input);
          if (need_arc)
            {
               printf(" arc file.\n");
               printf("Producing file with stem : %s\n", arc_output); 
            }
          else if (need_pdb)
            {
               printf(" pdb file.\n");
               printf("Producing file with stem : %s\n", pdb_output); 
            }
          else
            {
              printf("\nERROR : No format type defined for trajectory output\n");
              exit(0);
            }
        }

/******************************************************************/
/*********** read data file (name input by user on the ************/
/*********** command line).                            ************/
/******************************************************************/

been_before= FALSE;

printf("DEBUG>> Here\n");

if (is_car)
  {
   open_file( &fp_car_input, master_input, "r");

/**** For master input check for fixing flags in group labels ****/

    good_read=  read_car( fp_car_input, &header_line[0], &title_line[0],
                          &molecule[0], &date_line[0], &pbc, &num_atoms,
                          &num_of_mols, &num_mol_members[0], &mol_number[0],
                          &abc[0], &been_before, &groups, have_grp, 
                          &num_groups, TRUE, &fix_flags[0]);

/*    if (read_restart) */
/*      { */
/*        printf("ERROR: Restart directive given but master file is a car format file\n"); */
/*        printf("       Restart is only intended for copying fixed atom flags over \n"); */
/*        printf("       from POTCAR type files to CONTCAR type structures. \n"); */
/*        exit(0); */
/*      } */

    if (have_grp)
      {
         printf("GRUP defined as containing %d atoms:\n", groups.num_grp1+1);
         if (groups.num_grp2 > -1)
           {
              printf("Second group defined as containing %d atoms:\n", groups.num_grp2+1);
           }
         check1 = FALSE;
         check2 = FALSE;

         if (need_angle)
          {
            for (iloop=0; iloop<=groups.num_grp1; iloop++)
              {
                printf("%d %d >>%s<<",iloop, grp_mem_ind, molecule[grp_mem_ind].label );
                if (strcmp(molecule[grp_mem_ind].label, axis1_lab) == 0)
                   {
                     check1=TRUE;
                     axis1_ind = grp_mem_ind;
                     printf("   the first axis definition atom\n");
                   }
                else if (strcmp(molecule[grp_mem_ind].label, axis2_lab) == 0)
                   {
                     check2=TRUE;
                     axis2_ind = grp_mem_ind;
                     printf("   the second axis definition atom\n");
                   }
                else
                   {
                     printf("\n");
                   }
              }
            if (!check1 && !check2)
              {
                printf("ERROR: The axis definition atoms given in the input do not occur in GRUP\n");
                exit(0);
              }
            else if (!check1)
              {
                printf("ERROR: The first axis definition atom given in the input does not occur in GRUP\n");
                exit(0);
              }
            else if (!check2)
              {
                printf("ERROR: The second axis definition atom given in the input does not occur in GRUP\n");
                exit(0);
              }
          }
         else
          {
            for (iloop=0; iloop<=groups.num_grp1; iloop++)
              {
                grp_mem_ind= groups.group1[iloop];
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

    if (read_restart)
      {
        printf("ERROR: Restart directive given but master file is a gulp format file\n");
        printf("       Restart is only intended for copying fixed atom flags over \n");
        printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
        exit(0);
      }

   scale_factor=1.0;

   printf("For GULP inputs need to set is_cart or is_fract\n");

    fclose(fp_gulp_input);
  }
else if (is_vasp)
  {
    open_file( &fp_vasp_input, master_input, "r");

/*** Use been_before flag to pick up fixed atom flags from master ****/
    been_before=TRUE;
    
    good_read = read_poscar( fp_vasp_input, &title_line[0],
                             &molecule[0], &date_line[0], &pbc, &num_atoms,
                             &num_of_mols, &num_mol_members[0], &mol_number[0],
                             &latt_vec[0], &recip_latt_vec[0], &abc[0], 
                             &been_before, &is_fract, &is_cart,
                             &ion_number[0], &num_types, &scale_factor, &fix_flags[0]);

    num_mol_members[0]=num_atoms;
    pbc=TRUE;
    fclose(fp_vasp_input);

/***** Get atom label information from the POTCAR file for the job *****/

    open_file( &fp_vasp_input, potcar_input, "r");

    good_read = read_potcar( fp_vasp_input, &molecule[0], 
                             &ion_number[0], &num_types); 

    fclose(fp_vasp_input);

    if (read_restart)
      {
        printf("Reading co-ordinates from the restart file : %s\n", restart_file);
        printf("While taking flags for fixed atoms from the master file: %s\n", master_input);

        open_file( &fp_vasp_input, restart_file, "r");

/*** Use been_before flag to ignore fixed atom flags from restart file ****/
        been_before=FALSE;
    
        good_read = read_poscar( fp_vasp_input, &title_line[0],
                                 &molecule[0], &date_line[0], &pbc, &num_atoms,
                                 &num_of_mols, &num_mol_members[0], &mol_number[0],
                                 &latt_vec[0], &recip_latt_vec[0], &abc[0], 
                                 &been_before, &is_fract, &is_cart,
                                 &ion_number[0], &num_types, &scale_factor, &fix_flags[0]);

        if (num_atoms != num_mol_members[0])
          {
             printf("ERROR: master file and restart file contain different numbers of atoms\n");
             exit(0);
          }
      }
  }
/****************************************************************************/
/*** Process the OUTCAR file to make convergence movie DJW & EJ May 04 ******/
/****************************************************************************/
    else if (have_out)
      {
        open_file( &fp_outcar, master_input, "r");

        need_fermi=FALSE;
        good_read = read_outcar( fp_outcar, &molecule[0], &latt_vec[0],
                                 &recip_latt_vec[0], &abc[0], &eigenvalues[0],
                                 &eigenvecs[0], &forces, &chain_forces,
                                 &num_atoms, &types[0], &num_types, &num_modes, need_freq,
                                 need_force, need_energy, &have_band, &energy_vasp, need_fermi,
                                 &fermi);

        printf("DEBUG>> Back in main from read_outcar closing file pointer\n");

        fclose(fp_outcar);

       if (compare_modes) 
         {
           open_file( &fp_outcar, mode_compare_input, "r");

           need_fermi=FALSE;
           good_read = read_outcar( fp_outcar, &end_molecule[0], &new_latt_vec[0], 
                                    &new_recip_latt_vec[0], &new_abc[0], 
                                    &new_eigenvalues[0], &new_eigenvecs[0], &forces, 
                                    &chain_forces, &num_new_atoms, &end_types[0], 
                                    &num_new_types, &num_new_modes, need_freq,
                                    need_force, need_energy, &have_band, &energy_vasp, need_fermi,
                                    &fermi);

           fclose(fp_outcar);

           printf("Read reference modes from file : %s\n", mode_compare_input);

/*** simple checks ***/
           if ( num_atoms != num_new_atoms)
             {
               printf("ERROR: Comparing modes but number of atoms in two structures do not match.\n");
               exit(0);
             }
           else
             {
               printf("Comparing atoms by co-ordinates:\n");

               for (iatom=0; iatom<=num_atoms; iatom++)
                 {
                   dmin=-1.0;
                   map_atoms[iatom]=-1;
                   for (jatom=0; jatom<=num_atoms; jatom++)
                     {
                       dx = molecule[iatom].x - end_molecule[jatom].x;
                       dy = molecule[iatom].y - end_molecule[jatom].y;
                       dz = molecule[iatom].z - end_molecule[jatom].z;

                       min_image( &dx, &dy, &dz,
                                  &recip_latt_vec[0], &latt_vec[0]);
 
                       dist = sqrt(dx*dx + dy*dy + dz*dz);

/**** map_atoms matches atom indicies from master to reference i.e. map_atoms[iatom] is the ***/
/**** index for the atom in the reference list that matches iatom in the master list.       ***/
          
                       if (dmin < 0.0 || dist < dmin)
                         {
                            dmin=dist;
                            map_atoms[iatom] = jatom;
                         }
                     }
                   if (map_atoms[iatom] >= 0) 
                     {
                        printf("Closest atom to master atom %d in reference is %d they are %10.6f apart\n",
                                       iatom, map_atoms[iatom], dmin);
                     }
                   else
                     {
                        printf("ERROR: Problem mapping master atoms to reference.\n");
                        exit(0);
                     }
                 }
             }
           if ( num_modes != num_new_modes)
             {
               printf("ERROR: Comparing modes but number of modes in two structures do not match.\n");
               exit(0);
             }

           printf("Can compare eigenvecs as number of atoms and modes match.\n");
  
           printf("Comparing %d modes for %d atoms.\n", num_modes, num_atoms);

/*** Check reference modes are orthonormal ****/

           printf("Checking that reference modes form an orthonormal set.\n");
           for (imode=0; imode<num_modes; imode++)
             {
                tot_dot=0.0;
                for ( jmode=0; jmode < num_modes; jmode++)
                  {
                    dot=0.0;
                    for (iatom= 0; iatom <= num_atoms; iatom++)
                       { 
                         dot += new_eigenvecs[imode].dx[iatom]*new_eigenvecs[jmode].dx[iatom];
                         dot += new_eigenvecs[imode].dy[iatom]*new_eigenvecs[jmode].dy[iatom];
                         dot += new_eigenvecs[imode].dz[iatom]*new_eigenvecs[jmode].dz[iatom];
                       }
                    if (imode == jmode )
                      {
                        printf("Mode %d with self                    : %10.6f\n", 
                                                                         imode,dot);
                      }  
                    else
                      {
                        tot_dot+=fabs(dot);
                      }
                  }
                printf("Abs. total of %d mode with all others: %10.6f\n\n", imode, tot_dot);
             }
                
/*** Open file for csv output ***/
           open_file( &fp_outcar, "mode_match.csv", "w");
           fprintf(fp_outcar,"%s, ,%s\n",master_input, mode_compare_input);
           fprintf(fp_outcar,"Mode Num, WaveNo (cm-1), Mode Num, WaveNo (cm-1), diff, Composition\n"); 
/*** Work out vector distance between eigenvectors of two structures ***/
           for (imode=0; imode<num_modes; imode++)
             {
                dmin=-1.0;
                for ( jmode=0; jmode < num_modes; jmode++)
                  {
/*** Allow for eigenvectors out of phase ***/
                     for (isign=-1; isign <=1; isign+=2)
                       {
                         dist = 0.0;
                         for (iatom= 0; iatom <= num_atoms; iatom++)
                           { 
                             diffx=eigenvecs[imode].dx[iatom]+isign*new_eigenvecs[jmode].dx[map_atoms[iatom]];
                             diffy=eigenvecs[imode].dy[iatom]+isign*new_eigenvecs[jmode].dy[map_atoms[iatom]];
                             diffz=eigenvecs[imode].dz[iatom]+isign*new_eigenvecs[jmode].dz[map_atoms[iatom]];

                             dist += diffx*diffx + diffy*diffy + diffz*diffz;
                           }

/*         printf("Dist mode %d to %d is %10.6f\n",imode+1, jmode+1, sqrt(dist));   */
/*** test if this is the closest ***/
                        
                         if ( dmin < 0.0 || dist < dmin )
                           { 
                             dmin = dist;
                             match = jmode;
                           }
                      }
/*                    printf("\n"); */
                    dot_list[jmode]=0.0;
                    for (iatom= 0; iatom <= num_atoms; iatom++)
                      {
                        dot_list[jmode]+=eigenvecs[imode].dx[iatom]*new_eigenvecs[jmode].dx[map_atoms[iatom]];
                        dot_list[jmode]+=eigenvecs[imode].dy[iatom]*new_eigenvecs[jmode].dy[map_atoms[iatom]];
                        dot_list[jmode]+=eigenvecs[imode].dz[iatom]*new_eigenvecs[jmode].dz[map_atoms[iatom]];
                      }
                  }
               printf("Mode %d ( %7.1f cm-1) in %s matches %d ( %7.1f cm-1) in %s diff %5.3f : comp ",
                                           imode+1, eigenvalues[imode], master_input,
                                           match+1, new_eigenvalues[match], mode_compare_input,
                                           sqrt(dmin)); 

               fprintf(fp_outcar,"%d,%7.1f,%d,%7.1f,%5.3f", imode+1, eigenvalues[imode], 
                                                            match+1, new_eigenvalues[match],
                                                            sqrt(dmin)); 

               for ( jmode=0; jmode < num_modes; jmode++)
                  {
                    if (fabs(dot_list[jmode]) > 0.1) 
                      {
                        printf("%5.3f of %d ", dot_list[jmode], jmode+1);
                        fprintf(fp_outcar,",%5.3f,%d", dot_list[jmode], jmode+1);
                      }
                  }
               printf("\n");
               fprintf(fp_outcar,"\n");
             }
           fclose(fp_outcar);
           exit(0);
         }

    if (read_restart)
      {
        printf("ERROR: Restart directive given but master file is a OUTCAR format file\n");
        printf("       Restart is only intended for copying fixed atom flags over \n");
        printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
        exit(0);
      }

        num_of_mols = 1;
        num_atoms++; 
        num_mol_members[0]=num_atoms;
        scale_factor=1.0;
        is_cart = TRUE;
        is_fract= FALSE;
        pbc=TRUE;
  
      }

else if (is_siesta)
   {
      printf("Will try to get structure from siesta file\n");
      open_file( &fp_siesta_input, master_input, "r");

    if (read_restart)
      {
        printf("ERROR: Restart directive given but master file is a siesta fdf format file\n");
        printf("       Restart is only intended for copying fixed atom flags over \n");
        printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
        exit(0);
      }

      num_free=-1;
    
      good_read = read_fdf( fp_siesta_input, &title_line[0],
                            &molecule[0], &date_line[0], &pbc, &num_atoms,
                            &num_of_mols, &num_mol_members[0], &mol_number[0],
                            &latt_vec[0], &recip_latt_vec[0], &abc[0], 
                            &been_before, &is_fract, &is_cart,
                            &ion_number[0], &num_types, &scale_factor, &num_free,
                            &fix_flags[0]);

      scale_factor=1.0;

 /**** read_fdf returns the highest index of an atom in the list ***/
      num_atoms++;
      num_mol_members[0]=num_atoms;
      pbc=TRUE;
      fclose(fp_siesta_input);

      cut_to_dot(master_input);
      sprintf(master_input, "%s.vectors", master_input);

      if (need_freq)
        {
          printf("Now getting eigenvalues and vectors from %s file\n",master_input);
          open_file( &fp_siesta_input, master_input, "r");

          read_siesta_vectors( fp_siesta_input, &eigenvalues[0], 
                               &eigenvecs[0], &num_atoms, num_free,
                               &num_modes);
          fclose(fp_siesta_input);

          cut_to_dot(master_input);
          sprintf(master_input, "%s.fdf", master_input);
        }
   }

/*************************************************************************/
/**** End of special case stuff !! ***************************************/
/*************************************************************************/

else
  {
    printf("ERROR: Do not understand format of master file\n");
    exit(0);
  }

/*****************************************************************************/
/*** Convert fractional co-ordinates to cartesian  ***************************/
/*****************************************************************************/

printf("DEBUG>> Finished reading structural data\n");

if (is_fract)
  {
     icoord=0;

     printf("Will calculate cartesian using\n");
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

            printf("Sending frac to fract_to_cart: %10.6f %10.6f %10.6f\n", 
                                  molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);

            fract_to_cart( &(molecule[iloop].x), &(molecule[iloop].y), &(molecule[iloop].z),
                           fract_coords[icoord-3],fract_coords[icoord-2],fract_coords[icoord-1],
                           &latt_vec[0] );

            printf("Now get cartesian: %10.6f %10.6f %10.6f\n\n", 
                                    molecule[iloop].x, molecule[iloop].y, molecule[iloop].z);
         }
   }

printf("DEBUG>> report good_read\n");
if (good_read == 0)
  {

if (read_restart)
  {
     printf("Structure as read from restart file: %s\n",restart_file);
  }
else
  {
     printf("Structure as read from master file: %s\n",master_input);
  }

for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {
      printf("%s (%s) %10.6f %10.6f %10.6f\n",molecule[iatom].label,
                                              molecule[iatom].elem,
                                              molecule[iatom].x,
                                              molecule[iatom].y,
                                              molecule[iatom].z);
   }


/**** For car files generate cartesian lattice vectors from abc alpha beta gamma ****/
/**** do same for SIESTA fdf files, added Nov 06, Dave Willock.                   ****/

     if (pbc && (is_car || is_siesta)) 
       {
         cart_latt_vecs( &abc[0], &latt_vec[0], &recip_latt_vec[0]);
 
          printf("\nFractional co-ordinates of atoms in structure:\n");

         for (iatom= 0; iatom < num_mol_members[0]; iatom++)
            {
/**** Print out fractional co-ordinates version of structure ****/
              cart_to_fract( molecule[iatom].x, molecule[iatom].y, molecule[iatom].z,
                             &fa, &fb, &fc,
                             &recip_latt_vec[0] );

               printf("%s %10.6f %10.6f %10.6f\n",molecule[iatom].label,
                                                  fa, fb, fc);
           }
       }

     printf("Found:\n %d atoms \n %d molecules \n",num_atoms,num_of_mols);


     start_mol = 0;
     for (iloop = 0; iloop < num_of_mols; iloop++)
       {
          printf("Second off to generate neighbours with start_mol= %d, num_mol_members= %d\n",
                                                 start_mol,num_mol_members[iloop]);

          generate_neighbours( &molecule[start_mol], num_mol_members[iloop]-1,
                               &types[0], &num_types,
                               pbc, &recip_latt_vec[0],  &latt_vec[0],
                               &spec_charges[0]);

          printf("Molecule %d has %d members starting at %d\n",iloop,num_mol_members[iloop], start_mol);

          for (iatom= start_mol; iatom < start_mol+num_mol_members[iloop]; iatom++)
             {
/**** If no group number make one up ****/
/**** Dave December 2005             ****/

       if (strncmp(&(molecule[iatom].group_no[0])," ",1) == 0) sprintf(&(molecule[iatom].group_no[0]),"X"); 
       if (molecule[iatom].group_no[0] == '\0') sprintf(&(molecule[iatom].group_no[0]),"X"); 

                printf("%s (elem= %s, g_no >>%s<<) with %d neighbours : ", 
                                molecule[iatom].label, 
                                molecule[iatom].elem, 
                                molecule[iatom].group_no, 
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
printf("Back from find_mol\n");

printf("Structure as currently held: %s\n",master_input);
for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {
      printf("%s %10.6f %10.6f %10.6f\n",molecule[iatom].label,
                                         molecule[iatom].x,
                                         molecule[iatom].y,
                                         molecule[iatom].z);
   }


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
/*** Deal with the user requests *********************/
/*** At this point we have the co-ordinates read in  */
/*** from whatever source and neighbour and molecule */
/*** information assigned.                           */
/***                                                 */
/*** molecule[] : Atom co-ordinates and other info   */
/*** num_mol_members[0] : number of atoms in list    */
/*** latt_vec[] : a nine membered array containing   */
/***              the real space lattice vectors     */
/***   [0], [1], [2] : ax ay az                      */
/***   [3], [4], [5] : bx by bz                      */
/***   [6], [7], [8] : cx cy cz                      */
/*****************************************************/

/**** DEBUG Kara to insert new routine for layer finding here ***/

if (have_miller)
 {

   printf("Miller indices supplied as: %d %d %d\n", miller[0], miller[1], miller[2]);

   reorientate_cell(&molecule[0], &num_mol_members[0], &latt_vec[0], &recip_latt_vec[0],
                    &abc[0], &miller[0], &new_latt_vec[0], &new_recip_latt_vec[0],
                    &new_abc[0], &slab_mol[0], &num_slab_atoms);
   
look_for_layers = TRUE;
if (look_for_layers)
  {
/*    NEW ROUTINE TO BE ADDED HERE */
    if (need_car)
      {
         sprintf(step_output,"%s_%d%d%d.car",car_output,miller[0],miller[1],miller[2]);
         printf("Trying to write car file : %s\n", step_output);
         open_file( &fp_car_output, step_output, "w");

         start_frame = TRUE;
         write_car( fp_car_output, &header_line[0], &title_line[0], 
                    &c_title_line[0], &date_line[0], &slab_mol[0], 
                    &mol_number[0], pbc, &new_abc[0], 
                    num_slab_atoms, scale_factor, start_frame,
                    &super[0], &new_latt_vec[0], &new_recip_latt_vec[0],  
                    &fix_flags[0], &magmom[0], num_magmom);

         fclose(fp_car_output);

         sprintf(step_output,"%s_%d%d%d.pdb",car_output,miller[0],miller[1],miller[2]);
         printf("Trying to write pdb file : %s\n", step_output);
         open_file( &fp_pdb_output, step_output, "w");


         write_pdb(fp_pdb_output, &slab_mol[0], 
                      &new_abc[0], num_slab_atoms, &super[0],
                      &new_latt_vec[0]);
         fclose(fp_car_output);

      }
    exit(0);
  }

  }

/**** DEBUG Ends ************************************************/

if (is_fract && !need_mdtraj)
  {
     printf("Converting fractional to cartesian\n");
     is_fract=FALSE;
     is_cart=TRUE;

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


/*****************************************************/
/*** If requested write out start points *************/
/*****************************************************/

if (need_poscar)
  {
    open_file( &fp_vasp_output, poscar_output, "w");

    sort_by_elem( &molecule[0], num_mol_members[0], &types[0], num_types);

    printf("Writing new POSCAR format file: %s\n", poscar_output);

    write_poscar( fp_vasp_output, &molecule[0],  &fract_coords[0],
                  &types[0], num_types, &latt_vec[0], &scale_factor, num_atoms,
                  &title_line[0], &c_title_line[0], pbc, is_fract, &fix_flags[0]);

    fclose(fp_vasp_output);
  }

printf("Structure as currently held: %s\n",master_input);
for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {
      printf("%s %10.6f %10.6f %10.6f\n",molecule[iatom].label,
                                         molecule[iatom].x,
                                         molecule[iatom].y,
                                         molecule[iatom].z);
   }


if (need_car)
   {
      printf("Trying to write car file : %s\n", car_output);
      open_file( &fp_car_output, car_output, "w");

      start_frame = TRUE;
      write_car( fp_car_output, &header_line[0], &title_line[0], 
                 &c_title_line[0], &date_line[0], &molecule[0], 
                 &mol_number[0], pbc, &abc[0], 
                 num_mol_members[0], scale_factor, start_frame,
                 &super[0], &latt_vec[0], &recip_latt_vec[0],  &fix_flags[0], 
                 &magmom[0], num_magmom);

      fclose(fp_car_output);
   }

/**********************************************************************/
/** Process XDATCAR file to create a trajectory from MD run ***********/
/**********************************************************************/

if (need_mdtraj)
  {
    open_file( &fp_input_frame, mdtraj_input, "r");
    
    printf("Processing trajectory from file %s\n", mdtraj_input);

    if (need_pdb)
     {
        open_file( &fp_pdb_output, pdb_output, "w");
     }
    if (need_arc)
     {
        open_file( &fp_arc_output, arc_output, "w");
     }


    num_frames = -1;

/*** On first call read_xdatcar will send back number of frames *****/
/*** from header line of XDATCAR file                           *****/

    read_xdatcar(fp_input_frame, &num_frames, &step_molecule[0], num_atoms,
                 &fix_flags[0]);

    printf("Have %d frames to read\n", num_frames);
    printf("latt_vecs:\n");
    printf("%10.6f %10.6f %10.6f \n", latt_vec[0],  latt_vec[1],  latt_vec[2]);
    printf("%10.6f %10.6f %10.6f \n", latt_vec[3],  latt_vec[4],  latt_vec[5]);
    printf("%10.6f %10.6f %10.6f \n", latt_vec[6],  latt_vec[7],  latt_vec[8]);

/********************************************************************/
/** Read in frame by frame and process to required output format ****/
/********************************************************************/

    for (iatom=0; iatom < num_atoms; iatom++) 
                              step_molecule[iatom]=molecule[iatom];

    for (iframe=0; iframe < num_frames; iframe++)
      {
        read_xdatcar(fp_input_frame, &num_frames, 
                        &step_molecule[0], num_atoms,
                        &fix_flags[0]);

        printf("Read in frame %d of %d\n", iframe+1, num_frames);
        is_fract=TRUE;
        if (is_fract)
          {
             printf("Processing the fractional co-ords\n");
             for (iloop=0; iloop < num_atoms; iloop++)
                 {

                   vec[0] = step_molecule[iloop].x;
                   vec[1] = step_molecule[iloop].y;
                   vec[2] = step_molecule[iloop].z;

                   fract_to_cart( &(step_molecule[iloop].x),
                                  &(step_molecule[iloop].y), 
                                  &(step_molecule[iloop].z),
                                  vec[0], vec[1], vec[2],
                                  &latt_vec[0] );
                 }
          }

        if (need_pdb)
          {

            write_pdb(fp_pdb_output, &step_molecule[0], 
                      &abc[0], num_atoms, &super[0],
                      &latt_vec[0]);
          }

/*** This is how to write an arc file *****/

        if (need_arc)
          {
            printf("Writing frame to arc file....");
            if ( iframe == 0)
              {
                 printf("first frame\n");
                 sprintf(c_title_line,
                         "Inter_vasp generated trajectory file from MD run.");

                 start_frame=TRUE;
                 write_car( fp_arc_output, &header_line[0], 
                            &title_line[0], &c_title_line[0],
                            &date_line[0], &step_molecule[0], 
                            &mol_number[0], pbc, &abc[0], 
                            num_mol_members[0], scale_factor, start_frame,
                            &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
                            &magmom[0], num_magmom);
              }
            else
              {
                 printf("NOT first frame\n");
                 start_frame=FALSE;
                 write_car( fp_arc_output, &header_line[0], 
                            &title_line[0], &c_title_line[0],
                            &date_line[0], &step_molecule[0], 
                            &mol_number[0], pbc, &abc[0], 
                            num_mol_members[0], scale_factor, start_frame,
                            &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0],
                            &magmom[0], num_magmom);
              }
          }
      }

    if (need_pdb)
      {
         fclose(fp_pdb_output);
      }
    if (need_arc)
      {
         fclose(fp_arc_output);
      }
    fclose(fp_input_frame);
    exit(0);
  }
/**********************************************************************/
/**********************************************************************/
/**********************************************************************/

if (need_pdb)
   {
      open_file( &fp_pdb_output, pdb_output, "w");

      write_pdb(fp_pdb_output, &molecule[0], 
                &abc[0], num_atoms, &super[0],
                &latt_vec[0]);

      if (!need_interpolate) fclose(fp_pdb_output);                
   }

if (need_gulp)
   {
      open_file( &fp_gulp_output, gulp_output, "w");

      space_group=1;
      write_gulp(fp_gulp_output, &title_line[0], 
                 &molecule[0], num_atoms,
                 &shells[0], num_shells, GULP_CART, &abc[0],
                 space_group, &recip_latt_vec[0], &latt_vec[0],
                 &spec_charges[0], num_types, &super[0],
                 need_shells, &shell_species[0], num_shell_species);

      if (!need_interpolate) fclose(fp_gulp_output);                
   }

/*****************************************************/
/***** Read End point file for interpolation cases ***/
/*****************************************************/
if ( need_interpolate || need_react_coord ) 
  {
if (end_is_car)
  {
    open_file( &fp_car_input, end_input, "r");
    been_before=FALSE;

/*******************************************************/
/*** No need to re-establish group atoms at end point **/
/*** hence FALSE at end of call.                      **/
/*** Similarly fixed atoms need only be defined in    **/
/*** master.                                          **/
/*** This allows group defined in master only         **/
/*******************************************************/

    good_read=  read_car( fp_car_input, &end_header_line[0], &end_title_line[0],
                          &end_molecule[0], &end_date_line[0], &pbc, &end_num_atoms,
                          &end_num_of_mols, &end_num_mol_members[0], &end_mol_number[0],
                          &abc[0], &been_before, &groups, have_grp, &num_grps,
                          FALSE, &fix_flags[0]);

/****************************************************************************/
/*** Make end a minimum image with master added Dec 05, Dave Willock ********/
/*** Unless asked not to.                 added Nov 09, Dave Willock ********/
/****************************************************************************/

    if (pbc)
     {
         if (end_min_image)
           {
             for (iatom=0; iatom < num_mol_members[0]; iatom++)
              {
                 dx = end_molecule[iatom].x - molecule[iatom].x;
                 dy = end_molecule[iatom].y - molecule[iatom].y;
                 dz = end_molecule[iatom].z - molecule[iatom].z;

                 min_image( &dx, &dy, &dz,
                           &recip_latt_vec[0], &latt_vec[0]);

                 end_molecule[iatom].x = molecule[iatom].x + dx;
                 end_molecule[iatom].y = molecule[iatom].y + dy;
                 end_molecule[iatom].z = molecule[iatom].z + dz;
               }
            }

         for (iatom=0; iatom < num_mol_members[0]; iatom++)
           {

           d = sqrt(atom_separation_squared(&molecule[iatom], 
                                            &end_molecule[iatom], 
                                            FALSE, &recip_latt_vec[0], 
                                            &latt_vec[0]));

           printf("Atoms %d : %s (Start) and %s (End) are %10.6f apart\n",
                          iatom, molecule[iatom].label, 
                          end_molecule[iatom].label, d);
           }
     }

    is_cart = TRUE;
    is_fract= FALSE;
    fclose(fp_car_input);
  }
else
  {
    printf("Can only cope with end point car files ");
    printf("for interpolation or reaction coords at the moment....\n");
  }

/****************************************************************************/
/** Check consistency of end point and start point **************************/
/****************************************************************************/

   if (end_num_atoms != num_atoms)
     {
        printf("ERROR: Trying to interpolate structures with different numbers of atoms\n");
        printf("       Counted %d at start point but %d at end.\n", num_atoms, end_num_atoms);
        exit(0);
     }

/****************************************************************************/
/** Assemble neighbour info for end point ***********************************/
/****************************************************************************/

     start_mol = 0;
     for (iloop = 0; iloop < end_num_of_mols; iloop++)
       {
          printf("First off to generate neighbours for end point\n");

          generate_neighbours( &end_molecule[start_mol], 
                               end_num_mol_members[iloop]-1,
                               &end_types[0], &end_num_types,
                               pbc, &recip_latt_vec[0],  &latt_vec[0],
                               &spec_charges[0]);

          printf("Molecule %d in end point structure has %d members starts at %d\n",
                                         iloop,end_num_mol_members[iloop], start_mol);

          for (iatom= start_mol; iatom < start_mol+end_num_mol_members[iloop]; iatom++)
             {

                end_molecule[iatom].mass= 
                                 atomic_mass_list(end_molecule[iatom].elem);

                printf("%s (elem= %s, mass %10.6f) x: %10.6f y: %10.6f z: %10.6f, with %d neighbours : ",
                                                   end_molecule[iatom].label, 
                                                   end_molecule[iatom].elem, 
                                                   end_molecule[iatom].mass,
                                                   end_molecule[iatom].x,
                                                   end_molecule[iatom].y,
                                                   end_molecule[iatom].z,
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
/*** Convert fractional co-ordinates to cartesian for end point  **************/
/******************************************************************************/

if (is_fract)
  {
     printf("Converting fractional to cartesian\n");
     is_fract=FALSE;
     is_cart=TRUE;

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
/*** Carry out react_coords run, added December 05 Dave Willock *****************/
/********************************************************************************/

 if (need_react_coord)
   {
      react_coords(&molecule[0], &end_molecule[0], num_atoms, &image_files[0],
                   num_images, &recip_latt_vec[0], &latt_vec[0], pbc);
      exit(0);
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
                              num_mol_members[0], &chosen_indices[0],
                              &num_chosen_atoms, mol_ind, &mol_cnt_lab[0],
                              &start_cofm[0], &end_cofm[0], &inter_cofm[0]);

           printf("Back in main have:\n");
           printf("start_cofm %10.6f %10.6f %10.6f\n", start_cofm[0], start_cofm[1], start_cofm[2]);
           printf("end_cofm   %10.6f %10.6f %10.6f\n", end_cofm[0], end_cofm[1], end_cofm[2]);

/*****************************************************/
/*** For late transition states get bonding info *****/
/*****************************************************/

           any_rec =TRUE;

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
                        printf("This atom has lost bonds\n");
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
                        printf("This atom has gained bonds\n");
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
                        printf("This atom has maintained its number of partners\n");
                        for (ineigh=0; ineigh< end_molecule[iatom].num_neigh; ineigh++)
                          {
                            iatom2 = end_molecule[iatom].neighb[ineigh];

                            matched_neigh = FALSE;
                            for (ineigh2=0; ineigh2< molecule[iatom].num_neigh; ineigh2++)
                              {
                                 printf("Testing %s as partner of %s against %s as partner of %s\n",
                                               molecule[iatom].label, molecule[molecule[iatom].neighb[ineigh2]].label,
                                               end_molecule[iatom].label, end_molecule[iatom2].label);

                                 if (iatom2 == molecule[iatom].neighb[ineigh2])
                                    {
                                      printf("Thats a match\n");
                                      matched_neigh = TRUE;
                                    }
                              }

                            if (!matched_neigh)
                               {
                                 num_old_bonds++;
                                 num_new_bonds++;

                                 old_bonds[num_old_bonds].atom1 = iatom;
                                 old_bonds[num_old_bonds].atom2 =
                                                  molecule[iatom].neighb[ineigh];

                                 new_bonds[num_new_bonds].atom1 = iatom;
                                 new_bonds[num_new_bonds].atom2 = iatom2;
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
                   for (jloop=0; jloop <= num_old_bonds; jloop++)
                     {
                        if    (iatom == old_bonds[jloop].atom1
                            || iatom == old_bonds[jloop].atom2)
                          {
                            have_transfer=TRUE;
                            atom_transfered = iatom;
                            strcpy(elem_rec, molecule[iatom2].elem);
                          }
                        if    (iatom2 == old_bonds[jloop].atom1
                            || iatom2 == old_bonds[jloop].atom2)
                          {
                            have_transfer=TRUE;
                            atom_transfered = iatom2;
                            strcpy(elem_rec, molecule[iatom].elem);
                          }
                     }
                 }

               if (have_transfer)
                 {
                   printf("This is a transfer reaction with atom %d %s transfered\n",
                                 atom_transfered, molecule[atom_transfered].label);

                   printf("The receiving atom has element type %s\n", elem_rec);

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

                   min_image( &vec[0], &vec[1], &vec[2],
                              &recip_latt_vec[0], &latt_vec[0]);

                   new_bond_length= sqrt( vec[0] * vec[0] 
                                         +vec[1] * vec[1] 
                                         +vec[2] * vec[2] );

                   printf("Final length of bond formed = %10.6f Angstroms between %s and %s\n",
                                new_bond_length, end_molecule[atom_transfered].label, end_molecule[iatom].label);

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

         printf("Inter-structure vector:\n");
         for (iloop=0; iloop < num_atoms; iloop++)
           {
             inter_vec[iloop*3]   = end_molecule[iloop].x - molecule[iloop].x; 
             inter_vec[iloop*3+1] = end_molecule[iloop].y - molecule[iloop].y; 
             inter_vec[iloop*3+2] = end_molecule[iloop].z - molecule[iloop].z; 

             printf("%10.6f %10.6f %10.6f \n", inter_vec[iloop*3], 
                                     inter_vec[iloop*3+1], inter_vec[iloop*3+2]);
           }
         printf("\n");

/***********************************************/
/** Generate intermediates and output files ****/
/***********************************************/

         step= 1.0 / ( num_inter + 1 );
         printf("Interpolation step : %10.6f\n", step);

         if (have_grp)
           {
              if (need_angle)
                {
/************************************************************/
/*** Use rotation of GRUP atoms around a vector defined *****/
/*** by the two atom labels supplied.                   *****/
/************************************************************/

                   printf("Interpolating by angle about defined axis\n");
                   printf("Axis defined by atoms %d %s and %d %s\n",
                                axis1_ind, molecule[axis1_ind].label, 
                                axis2_ind, molecule[axis2_ind].label );

/************************************************************/
/* Define axis vector ***************************************/
/************************************************************/

                   origin[0] = molecule[axis1_ind].x;
                   origin[1] = molecule[axis1_ind].y;
                   origin[2] = molecule[axis1_ind].z;

                   axis[0] = molecule[axis2_ind].x-molecule[axis1_ind].x;
                   axis[1] = molecule[axis2_ind].y-molecule[axis1_ind].y;
                   axis[2] = molecule[axis2_ind].z-molecule[axis1_ind].z;

                   unit_vector(&axis[0]); 

/*************************************************************/
/* Work out vector from axis to each of the grouped atoms ****/
/* Check first which group has this type of iteration     ****/
/*************************************************************/
                  if (groups.group_type1 == ANGLE_TYPE) num_in_group = groups.num_grp1;
                  if (groups.group_type2 == ANGLE_TYPE) num_in_group = groups.num_grp2;
                     
                   for (igrp=0; igrp<=num_in_group; igrp++)
                     {
                       if (groups.group_type1 == ANGLE_TYPE) grp_mem_ind = groups.group1[igrp]; 
                       if (groups.group_type2 == ANGLE_TYPE) grp_mem_ind = groups.group2[igrp]; 

                       if (grp_mem_ind != axis1_ind && grp_mem_ind != axis2_ind)   
                          {
/*************************************************************/
/** vector from atom to one of axis defining atoms which *****/
/** must be on the rotation axis.                        *****/
/*************************************************************/
                             vec[0] = molecule[grp_mem_ind].x - molecule[axis1_ind].x;
                             vec[1] = molecule[grp_mem_ind].y - molecule[axis1_ind].y;
                             vec[2] = molecule[grp_mem_ind].z - molecule[axis1_ind].z;

                             dot =   vec[0] * axis[0]
                                   + vec[1] * axis[1]
                                   + vec[2] * axis[2];

/*************************************************************/
/** Take component which is along axis and remove to get *****/
/** vector to atom perpendicular to the axis.            *****/
/*************************************************************/
                             vec1[0]= vec[0] - dot * axis[0];
                             vec1[1]= vec[1] - dot * axis[1];
                             vec1[2]= vec[2] - dot * axis[2];

                             unit_vector(&vec1[0]); 

/*************************************************************/
/** Repeat for the end_point molecule                    *****/
/*************************************************************/
                             vec[0] = end_molecule[grp_mem_ind].x - end_molecule[axis1_ind].x;
                             vec[1] = end_molecule[grp_mem_ind].y - end_molecule[axis1_ind].y;
                             vec[2] = end_molecule[grp_mem_ind].z - end_molecule[axis1_ind].z;

                             dot =   vec[0] * axis[0]
                                   + vec[1] * axis[1]
                                   + vec[2] * axis[2];

                             vec2[0]= vec[0] - dot * axis[0];
                             vec2[1]= vec[1] - dot * axis[1];
                             vec2[2]= vec[2] - dot * axis[2];

                             unit_vector(&vec2[0]); 

/*************************************************************/
/*** Total angle to rotate by is then obtained from the ******/
/*** vec1, vec2 dot product.                            ******/
/*************************************************************/

                             dot =   vec1[0] * vec2[0]
                                   + vec1[1] * vec2[1]
                                   + vec1[2] * vec2[2];
                            
                             molecule[grp_mem_ind].theta = -acos(dot);

                             printf("Atom %d %s will be rotated through %10.6f degrees during interpolation\n",
                                                   grp_mem_ind, molecule[grp_mem_ind].label, 
                                                   RAD_TO_DEG*molecule[grp_mem_ind].theta); 

                          }
                     }
                }
              else
                {
                  printf("\n");
/********************************************************************/ 
/*** For the case of grouped atoms work out changes of bond length **/
/*** May have two groups, treat seperately.                        **/
/********************************************************************/ 
                  for (igrp=0; igrp<=groups.num_grp1; igrp++)
                    {
                      if (groups.group_type1 == CENTRE_TYPE)
                        {
                          grp_mem_ind = groups.group1[igrp];  
                          printf("PICKED up index %d as atom %s\n",grp_mem_ind, molecule[grp_mem_ind].label);
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

                               orig_bond1[igrp]= size_vector(&vec[0]);

                               delta_bond1[igrp]= size_vector(&vec1[0])
                                                 -size_vector(&vec[0]);

                               printf("centre to atom %d vector is %10.6f %10.6f %10.6f\n",
                                             grp_mem_ind,
                                             vec[0],vec[1],vec[2]);
 
                           printf("This atom is %s the centre is %s\n", molecule[grp_mem_ind].label, 
                                                                        molecule[grp_cnt_ind].label);
                           printf("Original bond length for atom  %d is %10.6f\n", grp_mem_ind, orig_bond1[igrp]);
                           printf("Change in bond length for atom %d is %10.6f\n", grp_mem_ind, delta_bond1[igrp]);
                           }
                        }
                    }
/**** If have second group sort that out now ****/
                  if (num_grps > 0 && groups.group_type1 == CENTRE_TYPE )
                    {
                       printf("\nHave two groups to work with for the second:\n");
                       for (igrp=0; igrp<=groups.num_grp2; igrp++)
                         {
                           grp_mem_ind = groups.group2[igrp];  
                           if (grp_mem_ind != grp_cnt_ind2)
                             {
                                vec[0] =  molecule[grp_cnt_ind2].x 
                                        - molecule[grp_mem_ind].x; 

                                vec[1] =  molecule[grp_cnt_ind2].y 
                                        - molecule[grp_mem_ind].y; 

                                vec[2] =  molecule[grp_cnt_ind2].z 
                                        - molecule[grp_mem_ind].z; 

                                vec1[0] = end_molecule[grp_cnt_ind2].x 
                                           - end_molecule[grp_mem_ind].x; 

                                vec1[1] = end_molecule[grp_cnt_ind2].y 
                                           - end_molecule[grp_mem_ind].y; 

                                vec1[2] = end_molecule[grp_cnt_ind2].z 
                                           - end_molecule[grp_mem_ind].z; 

                                orig_bond2[igrp]= size_vector(&vec[0]);

                                delta_bond2[igrp]= size_vector(&vec1[0])
                                                  -size_vector(&vec[0]);

                                 printf("centre to atom %d vector is %10.6f %10.6f %10.6f\n",
                                             grp_mem_ind,
                                             vec[0],vec[1],vec[2]);
 
                                printf("This atom is %s the centre is %s\n", molecule[grp_mem_ind].label, 
                                                                        molecule[grp_cnt_ind2].label);
                                printf("Original bond length for atom  %d is %10.6f\n", grp_mem_ind, orig_bond2[igrp]);
                                printf("Change in bond length for atom %d is %10.6f\n", grp_mem_ind, delta_bond2[igrp]);
                            }
                        }
                    }
                  printf("\n\nData for grouped interpolation set up\n\n");
                }
           }

/******************************************************************/

         iframe=0;
         for (iloop=1; iloop <= num_inter+1; iloop++)
            {
               iframe++;
               if ( have_grp && need_angle )
                 {
                   for (iatom=0; iatom < num_atoms; iatom++)
                     {
/************************************************************/
/**** Check if the atom is a member of the group defined ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
                       check = FALSE;
                       for (igrp=0; igrp<=num_in_group; igrp++)
                          {
 /*                            grp_mem_ind = group_indices[igrp];  */
                             if ( grp_mem_ind == iatom ) check = TRUE;
                          }

                       if ( !check || iatom == axis1_ind || iatom == axis2_ind )
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
/*** Rotate each one by an increment of its total theta.      ***/ 
/****************************************************************/ 
                   for (igrp=0; igrp<=num_in_group; igrp++)
                     {
  /*                     grp_mem_ind = group_indices[igrp]; */

                       printf("Using angle interpolation for atom %d\n",grp_mem_ind);
                       if (grp_mem_ind != axis1_ind && grp_mem_ind != axis2_ind )
                         {
                            step_molecule[grp_mem_ind] = molecule[grp_mem_ind];

                       theta= iloop * molecule[grp_mem_ind].theta/ (num_inter+1);
  
                       for (iii=0; iii < num_atoms; iii++) flags[iii]= FALSE;
                       flags[grp_mem_ind]= TRUE;

                       rotate_with_flags(&step_molecule[0], &axis[0], &origin[0],
                                         theta, &flags[0], num_atoms-1);

                         }
                     }
                 }
               else if (have_grp)
                 {
                   for (iatom=0; iatom < num_atoms; iatom++)
                     {
/************************************************************/
/**** Check if the atom is a member of a group defined   ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
                       check = FALSE;
                       is_centre = iatom == grp_cnt_ind;
                       for (igrp=0; igrp<=groups.num_grp1; igrp++)
                          {
                             grp_mem_ind = groups.group1[igrp];  
                             if ( grp_mem_ind == iatom ) check = TRUE;
                          }
                       if (num_grps>0)
                          {
                             is_centre = iatom == grp_cnt_ind || iatom == grp_cnt_ind2;
                             for (igrp=0; igrp<=groups.num_grp2; igrp++)
                                {
                                   grp_mem_ind = groups.group2[igrp];  
                                   if ( grp_mem_ind == iatom ) check = TRUE;
                                }
                          }

                       if ( !check || is_centre)
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
                   num_in_group = groups.num_grp1; 
                   if (num_grps > 0) num_in_group = groups.num_grp1 +groups.num_grp2+1; 
                         
                   for (igrp=0; igrp<=num_in_group; igrp++)
                     {
/** test which group this atom belongs to ***/
                       if (igrp <= groups.num_grp1)
                         {
                           grp_mem_ind = groups.group1[igrp];  
                           is_centre = grp_mem_ind == grp_cnt_ind;
                           grp_cnt= grp_cnt_ind;
                           this_delta = delta_bond1[igrp];
                           this_orig  = orig_bond1[igrp];
                         } 
                       else
                         {
                           jgrp=igrp-groups.num_grp1-1;  
                           grp_mem_ind = groups.group2[jgrp];
                           is_centre = grp_mem_ind == grp_cnt_ind2;
                           grp_cnt= grp_cnt_ind2;
                           this_delta = delta_bond2[jgrp];
                           this_orig  = orig_bond2[jgrp];
                         }

                       printf("Using different interpolation for atom %d\n",grp_mem_ind);
                       if (!is_centre)
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
                            new_bond = this_orig + step * iloop * this_delta;

                            vec[0] = step_molecule[grp_mem_ind].x 
                                             - step_molecule[grp_cnt].x;
                            vec[1] = step_molecule[grp_mem_ind].y 
                                             - step_molecule[grp_cnt].y;
                            vec[2] = step_molecule[grp_mem_ind].z 
                                             - step_molecule[grp_cnt].z;

                            printf("Have member %d and centre %d\n", grp_mem_ind, grp_cnt);
                            printf("Centre to member vector: %10.6f %10.6f %10.6f\n", vec[0], vec[1], vec[2]);
                            printf("orig: %10.6f delta: %10.6f new: %10.6f\n", this_orig, this_delta, new_bond);

                            unit_vector(&vec[0]); 

                            step_molecule[grp_mem_ind].x = step_molecule[grp_cnt].x 
                                                          + new_bond * vec[0]; 

                            step_molecule[grp_mem_ind].y = step_molecule[grp_cnt].y 
                                                          + new_bond * vec[1]; 

                            step_molecule[grp_mem_ind].z = step_molecule[grp_cnt].z 
                                                          + new_bond * vec[2]; 

                            printf("New coords: %10.6f %10.6f %10.6f \n", step_molecule[grp_mem_ind].x,
                                                                          step_molecule[grp_mem_ind].y,
                                                                          step_molecule[grp_mem_ind].z);
                         }
                     }
                 }
               else if ( have_mol && !need_late )
                 {

                   cofm_step[0] = start_cofm[0] + step * iloop * inter_cofm[0];
                   cofm_step[1] = start_cofm[1] + step * iloop * inter_cofm[1];
                   cofm_step[2] = start_cofm[2] + step * iloop * inter_cofm[2];

                   for (iatom=0; iatom < num_atoms; iatom++)
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
                           step_molecule[iatom] = molecule[iatom];
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
                           step_molecule[iatom] = molecule[iatom];
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
                       
                       if (iframe != 1)
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
                /*           cofm_step[0] = start_cofm[0] + step * iloop * inter_cofm[0]; */
                /*           cofm_step[1] = start_cofm[1] + step * iloop * inter_cofm[1]; */
                /*           cofm_step[2] = start_cofm[2] + step * iloop * inter_cofm[2]; */
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

                               distance_gone = (iloop-1)* step;
                               distance_left = (num_inter-iloop+2)* step;

                               printf("At Switch over distance_gone %10.6f distance_left %10.6f\n",
                                            distance_gone, distance_left);

                               if (mag_after_switch > 0.0)
                                {
                                  printf("Changing step by factor %10.6f\n",mag_after_switch);
                                }
                               else
                                {
                                  printf("ERROR: No magnification factor defined for after switch\n");
                                  exit(0);
                                }

/****************************************************/
/** Re-define loop variables and molecule ***********/
/****************************************************/
                               iloop = 1;
                               num_inter= distance_left/(mag_after_switch*step); 

                               printf("Will now do another %d steps\n", num_inter);

                               step= 1.0 / ( num_inter + 1 );
                        
                               for (iatom=0; iatom < num_atoms; iatom++)
                                 {
                                   molecule[iatom] = step_molecule[iatom];
                                 }
                               for (imol=0; imol<=num_chosen_atoms; imol++)
                                 {
                                   iatom = chosen_indices[imol];
                                   molecule[iatom].x = molecule[iatom].x - cofm_step[0];
                                   molecule[iatom].y = molecule[iatom].y - cofm_step[1];
                                   molecule[iatom].z = molecule[iatom].z - cofm_step[2];
                                 }
                             }
                          
                           cofm_step[0] = start_cofm[0] + step * iloop * transfer_vec[0];
                           cofm_step[1] = start_cofm[1] + step * iloop * transfer_vec[1];
                           cofm_step[2] = start_cofm[2] + step * iloop * transfer_vec[2];
                         }

                     }

                   for (iatom=0; iatom < num_atoms; iatom++)
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
                          step_molecule[iatom] = molecule[iatom];
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
                          step_molecule[iatom] = molecule[iatom];
                          step_molecule[iatom].x = cofm_step[0] + molecule[iatom].x;
                          step_molecule[iatom].y = cofm_step[1] + molecule[iatom].y;
                          step_molecule[iatom].z = cofm_step[2] + molecule[iatom].z;
                        }
                    }

                   if ( any_rec )
                     {
/******************************************************/
/** Any surface atom of the right element type ********/
/** will do for the switch, find closest       ********/
/******************************************************/

                        latest_bond = -1.0;
                        for (iatom=0; iatom < num_atoms; iatom++) 
                          {
                            if (strcmp(step_molecule[iatom].elem, elem_rec) == 0)
                              {
                                 
                                 vec[0] = step_molecule[atom_transfered].x
                                         -step_molecule[iatom].x;
 
                                 vec[1] = step_molecule[atom_transfered].y
                                         -step_molecule[iatom].y;
 
                                 vec[2] = step_molecule[atom_transfered].z
                                         -step_molecule[iatom].z;

                                 temp_bond = sqrt( vec[0] * vec[0]
                                                  +vec[1] * vec[1]
                                                  +vec[2] * vec[2] );

                                 if (latest_bond < 0 || temp_bond < latest_bond) 
                                   {
                                     latest_bond = temp_bond;
                                     iatom2 = iatom;
                                   }
                              }
                          }
                        if (latest_bond < 0)
                          {
                             printf("ERROR: Cannot find any suitable receiving element types for late transition state.\n");
                             exit(0);
                          }
                     }
                   else
                     {
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
                    }

                   printf("latest bond length %10.6f receiver %s\n", latest_bond, step_molecule[iatom2].label);
                 }
/************************************************************/
/*** Otherwise this is a straight linear interpolation ******/
/************************************************************/
               else
                 {
                   for (iatom=0; iatom < num_atoms; iatom++)
                     {
                       step_molecule[iatom] = molecule[iatom];
                       printf("Atom at start : %s %10.6f %10.6f %10.6f\n",
                              step_molecule[iatom].label, step_molecule[iatom].x,
                              step_molecule[iatom].y,     step_molecule[iatom].z);
                                    
                       step_molecule[iatom].x = step_molecule[iatom].x 
                                       + step * iloop * inter_vec[iatom*3];
                       step_molecule[iatom].y = step_molecule[iatom].y 
                                       + step * iloop * inter_vec[iatom*3+1];
                       step_molecule[iatom].z = step_molecule[iatom].z 
                                       + step * iloop * inter_vec[iatom*3+2];

                       printf("Atom after move: %s %10.6f %10.6f %10.6f\n\n",
                              step_molecule[iatom].label, step_molecule[iatom].x,
                              step_molecule[iatom].y,     step_molecule[iatom].z);
                                    
                     }
                 }

               sprintf(step_output,"%s%d", poscar_output, iframe);

               open_file( &fp_vasp_output, step_output, "w");

               if ( iloop < num_inter+1)
                 {
/***************************************/
/*** Check sum of position vectors *****/
/***************************************/
                   if (need_interpolate)
                    {
                       if (need_shift)
                          {
                             for ( iatom=0; iatom < 3; iatom++) sum_pos[iatom] = 0.0;

                             for (iatom=0; iatom < num_atoms; iatom++)
                               {
                                  sum_pos[0] += step_molecule[iatom].x; 
                                  sum_pos[1] += step_molecule[iatom].y; 
                                  sum_pos[2] += step_molecule[iatom].z; 
                               }
               
                             printf("Sum of components for step %d : %10.6f %10.6f %10.6f\n",
                                    iloop, sum_pos[0], sum_pos[1], sum_pos[2]);

                             if (iloop == 1)
                               {
                                  printf("This is first image recording reference centre\n");
                                  sum_pos1[0] = sum_pos[0];
                                  sum_pos1[1] = sum_pos[1];
                                  sum_pos1[2] = sum_pos[2];
                               }
                             else
                               {
                                  shift_centre[0] = (sum_pos1[0]-sum_pos[0])/num_atoms;
                                  shift_centre[1] = (sum_pos1[1]-sum_pos[1])/num_atoms;
                                  shift_centre[2] = (sum_pos1[2]-sum_pos[2])/num_atoms;

                                  for ( iatom=0; iatom < 3; iatom++) sum_check[iatom] = 0.0;

                                  for (iatom=0; iatom < num_atoms; iatom++)
                                    {
                                       step_molecule[iatom].x += shift_centre[0]; 
                                       step_molecule[iatom].y += shift_centre[1]; 
                                       step_molecule[iatom].z += shift_centre[2]; 
 
                                       sum_check[0] += step_molecule[iatom].x; 
                                       sum_check[1] += step_molecule[iatom].y; 
                                       sum_check[2] += step_molecule[iatom].z; 
                                    }
                                   printf("Sum of components for step %d after shift : %10.6f %10.6f %10.6f\n",
                                                      iloop, sum_check[0], sum_check[1], sum_check[2]);
                               }
                           }
                        }

                   write_poscar( fp_vasp_output, &step_molecule[0],  
                                 &fract_coords[0], &types[0], num_types, 
                                 &latt_vec[0], &scale_factor, num_atoms,
                                 &title_line[0], &c_title_line[0], pbc, 
                                 is_fract, &fix_flags[0]);
                 }
               else
                 {
                   for (imol=0; imol<=num_chosen_atoms; imol++)
                     {
                       iatom = chosen_indices[imol];
                       end_molecule[iatom].x = end_molecule[iatom].x + end_cofm[0];
                       end_molecule[iatom].y = end_molecule[iatom].y + end_cofm[1];
                       end_molecule[iatom].z = end_molecule[iatom].z + end_cofm[2];
                     }

                   if (need_interpolate)
                    {
/***************************************/
/*** Check sum of position vectors *****/
/***************************************/
                      for ( iatom=0; iatom < 3; iatom++) sum_pos[iatom] = 0.0;

                      for (iatom=0; iatom < num_atoms; iatom++)
                        {
                           sum_pos[0] += step_molecule[iatom].x; 
                           sum_pos[1] += step_molecule[iatom].y; 
                           sum_pos[2] += step_molecule[iatom].z; 
                        }
               
                      printf("Sum of components for step %d : %10.6f %10.6f %10.6f\n",
                              iloop, sum_pos[0], sum_pos[1], sum_pos[2]);

                      if (iloop == 1)
                        {
                           printf("This is first image recording reference centre\n");
                           sum_pos1[0] = sum_pos[0];
                           sum_pos1[1] = sum_pos[1];
                           sum_pos1[2] = sum_pos[2];
                        }
                      else
                        {
                           shift_centre[0] = (sum_pos1[0]-sum_pos[0])/num_atoms;
                           shift_centre[1] = (sum_pos1[1]-sum_pos[1])/num_atoms;
                           shift_centre[2] = (sum_pos1[2]-sum_pos[2])/num_atoms;

                           for ( iatom=0; iatom < 3; iatom++) sum_check[iatom] = 0.0;

                           for (iatom=0; iatom < num_atoms; iatom++)
                             {
                                step_molecule[iatom].x += shift_centre[0]; 
                                step_molecule[iatom].y += shift_centre[1]; 
                                step_molecule[iatom].z += shift_centre[2]; 
 
                                sum_check[0] += step_molecule[iatom].x; 
                                sum_check[1] += step_molecule[iatom].y; 
                                sum_check[2] += step_molecule[iatom].z; 
                             }

                         printf("Sum of components for step %d after shift : %10.6f %10.6f %10.6f\n",
                                           iloop, sum_check[0], sum_check[1], sum_check[2]);

                        }
                     }

                   write_poscar( fp_vasp_output, &end_molecule[0],  
                                 &fract_coords[0], &types[0], num_types, 
                                 &latt_vec[0], &scale_factor, num_atoms,
                                 &title_line[0], &c_title_line[0], pbc, 
                                 is_fract, &fix_flags[0]);

               fclose(fp_vasp_output);
                }

/***** write arc file if requested *************************/

               if (need_arc)
                 {
                   if ( iloop == 1 && !done_start_frame)
                     {
                        sprintf(c_title_line,"Inter_vasp generated arcfile record of interpolation.");
                        start_frame=TRUE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], &step_molecule[0], 
                                   &mol_number[0], pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
                                   &magmom[0], num_magmom);

                        done_start_frame = TRUE;
                     }
                   else if ( iloop < num_inter+1)
                     {
                        start_frame=FALSE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], &step_molecule[0], 
                                   &mol_number[0], pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
                                   &magmom[0], num_magmom);
                     }
                   else
                     {
                        start_frame=FALSE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], &end_molecule[0], 
                                   &mol_number[0], pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
                                   &magmom[0], num_magmom);
                     }
                      
                 }

/***** write pdb file if requested *************************/

               if (need_pdb)
                 {
                   if ( iloop < num_inter+1)
                     {
                       write_pdb(fp_pdb_output, &step_molecule[0], 
                             &abc[0], num_atoms, &super[0], &latt_vec[0]);
                     }
                   else
                     {
                       write_pdb(fp_pdb_output, &end_molecule[0], 
                             &abc[0], num_atoms, &super[0], &latt_vec[0]);
                     }
                      
                 }

            if (need_car)
              {
/***** Strip .car from file to allow insertion of number ***/

               iletter=0;
               strcpy(step_output, car_output);
               while (strncmp(&step_output[iletter],".",1) != 0) iletter++;
               step_output[iletter]= '\0';

               sprintf(step_output,"%s%d.car", step_output, iframe);

               open_file( &fp_car_output, step_output, "w");

               start_frame=TRUE;

               if ( iloop < num_inter+1)
                 {
                   write_car( fp_car_output, &header_line[0], &title_line[0], 
                              &c_title_line[0], &date_line[0], &step_molecule[0], 
                              &mol_number[0], pbc, &abc[0], 
                              num_mol_members[0], scale_factor, start_frame,
                              &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0],  
                              &magmom[0], num_magmom);
                 }
               else
                 {
                   write_car( fp_car_output, &header_line[0], &title_line[0], 
                              &c_title_line[0], &date_line[0], &end_molecule[0], 
                              &mol_number[0], pbc, &abc[0], 
                              num_mol_members[0], scale_factor, start_frame,
                              &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0],  
                              &magmom[0], num_magmom);
                 }

               fclose(fp_car_output);
             }
           }
      if (need_arc) fclose(fp_arc_output);
   }
/****************************************/
/** Write out vibrational modes in ******/
/** Requested format               ******/
/****************************************/
    if (need_arc && need_freq)
      {
          if (which_mode<0)
              printf("Will write %d modes\n", num_modes);
          else
              printf("Have %d modes but just writing mode %d\n", num_modes, which_mode);

          for ( imode = 0; imode < num_modes; imode++)
            {
              if (which_mode < 0 || which_mode == imode+1)
                {
              sprintf(frame_output,"%s%d.arc",arc_output,imode+1);

              sprintf(c_title_line,"frame 1 Mode %d frequency %10.6f cm-1",
                                             imode+1, eigenvalues[imode]);

              open_file( &fp_arc_output, frame_output, "w");

              if (need_pdb)
               {
                 sprintf(frame_output_pdb,"%s%d.pdb",arc_output,imode+1);
                 printf("Will write pdbfile of mode to %s\n",frame_output_pdb);
                 open_file( &fp_pdb_output, frame_output_pdb, "w");

                 fprintf( fp_pdb_output,"COMENT %s\n", c_title_line);
               }

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

                  for (iatom=0; iatom < num_atoms; iatom++)
                    {
                      step_molecule[iatom] = molecule[iatom];
                      step_molecule[iatom].x = step_molecule[iatom].x 
                                            + step * eigenvecs[imode].dx[iatom];
                      step_molecule[iatom].y = step_molecule[iatom].y 
                                            + step * eigenvecs[imode].dy[iatom];
                      step_molecule[iatom].z = step_molecule[iatom].z 
                                            + step * eigenvecs[imode].dz[iatom];
                    }

                  if (need_poscar)
                    {
                      sprintf(frame_output,"%sMode_%d_step%d",
                                      poscar_output, imode+1, iloop);

                      open_file( &fp_vasp_output, frame_output, "w");
      			printf("Trying to write poscar file : %s\n", poscar_output);

                      sort_by_elem( &step_molecule[0], num_mol_members[0], 
                                    &types[0], num_types);

                      write_poscar( fp_vasp_output, &step_molecule[0],  
                                    &fract_coords[0], &types[0], 
                                    num_types, &latt_vec[0], 
                                    &scale_factor, num_mol_members[0],
                                    &title_line[0], &c_title_line[0], 
                                    pbc, is_fract, &fix_flags[0]);

                      fclose(fp_vasp_output);
                    }

		if (need_car)
 		  {
			sprintf(frame_output,"%sMode_%d_step%d.car",
                                      car_output, imode+1, iloop);
      			printf("Trying to write car file : %s\n", car_output);
      			open_file( &fp_car_output, frame_output, "w");

   			write_car( fp_car_output, &header_line[0], &title_line[0],
                 		&c_title_line[0], &date_line[0], &step_molecule[0],
                 		&mol_number[0], pbc, &abc[0],
                 		num_mol_members[0], scale_factor, TRUE,
                                &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
                                &magmom[0], num_magmom);

      fclose(fp_car_output);
   }


/***** write pdb file if requested *************************/

                  if (need_pdb)
                    {
                      write_pdb(fp_pdb_output, &step_molecule[0], 
                                &abc[0], num_atoms, &super[0], &latt_vec[0]);
                    }

                  write_car( fp_arc_output, &header_line[0], 
                             &title_line[0], &c_title_line[0],
                             &date_line[0], &step_molecule[0], 
                             &mol_number[0], pbc, &abc[0], 
                             num_mol_members[0], scale_factor, start_frame,
                             &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
                             &magmom[0], num_magmom);

/*********************************************/
/** Remove repeats of frequency information **/
/*********************************************/

                  sprintf(c_title_line,"frame %d ", iloop+2); 

                  start_frame = FALSE;
                }
              fclose(fp_arc_output);
              if (need_pdb) fclose(fp_pdb_output);
            }
         }
       }

/*********************************************/
/*** Give force data *************************/
/*********************************************/
     if (need_force)
       {
          printf("Forces read for %d atoms from OUTCAR file (eV/Angstrom)\n",
                               num_atoms);
          rms_x=0.0;
          rms_y=0.0;
          rms_z=0.0;
          rms_t=0.0;
          for (iatom=0; iatom < num_atoms; iatom++)
            {
               printf("%4s : %10.6f %10.6f %10.6f\n", molecule[iatom].label,
                                                      forces.dx[iatom],
                                                      forces.dy[iatom],
                                                      forces.dz[iatom]);
              
               rms_x += forces.dx[iatom] * forces.dx[iatom];
               rms_y += forces.dy[iatom] * forces.dy[iatom];
               rms_z += forces.dz[iatom] * forces.dz[iatom];
            }
          rms_x = sqrt( rms_x / num_atoms);
          rms_y = sqrt( rms_y / num_atoms);
          rms_z = sqrt( rms_z / num_atoms);

          rms_t = sqrt( rms_x * rms_x + rms_y * rms_y + rms_z * rms_z );

          printf("\n RMS : %10.6f %10.6f %10.6f\n", rms_x, rms_y, rms_z);
          printf("\n RMS magnitude : %10.6f\n\n", rms_t);

          if (have_band)
             {
                printf("This OUTCAR is one of a set from an elastic band calculation\n");

                printf("Chain spring forces: \n");
                rms_x=0.0;
                rms_y=0.0;
                rms_z=0.0;
                rms_t=0.0;
                for (iatom=0; iatom < num_atoms; iatom++)
                  {
                     printf("%4s : %10.6f %10.6f %10.6f\n", molecule[iatom].label,
                                                            chain_forces.dx[iatom],
                                                            chain_forces.dy[iatom],
                                                            chain_forces.dz[iatom]);
              
                     rms_x += chain_forces.dx[iatom] * chain_forces.dx[iatom];
                     rms_y += chain_forces.dy[iatom] * chain_forces.dy[iatom];
                     rms_z += chain_forces.dz[iatom] * chain_forces.dz[iatom];
                  }
                rms_x = sqrt( rms_x / num_atoms);
                rms_y = sqrt( rms_y / num_atoms);
                rms_z = sqrt( rms_z / num_atoms);

                rms_t = sqrt( rms_x * rms_x + rms_y * rms_y + rms_z * rms_z );

                printf("\n Tangent RMS : %10.6f %10.6f %10.6f\n", rms_x, rms_y, rms_z);
                printf("\n Tangent RMS magnitude : %10.6f\n\n", rms_t);

                printf("Actual forces in chain direction ( without band forces )\n");
                rms_x=0.0;
                rms_y=0.0;
                rms_z=0.0;
                rms_t=0.0;
                for (iatom=0; iatom < num_atoms; iatom++)
                  {
                     forces.dx[iatom]=forces.dx[iatom]-chain_forces.dx[iatom];
                     forces.dx[iatom]=forces.dy[iatom]-chain_forces.dy[iatom];
                     forces.dx[iatom]=forces.dz[iatom]-chain_forces.dz[iatom];

                     printf("%4s : %10.6f %10.6f %10.6f\n", molecule[iatom].label,
                                                            forces.dx[iatom],
                                                            forces.dy[iatom],
                                                            forces.dz[iatom]);
              
                     rms_x += forces.dx[iatom] * forces.dx[iatom];
                     rms_y += forces.dy[iatom] * forces.dy[iatom];
                     rms_z += forces.dz[iatom] * forces.dz[iatom];
                  }
                rms_x = sqrt( rms_x / num_atoms);
                rms_y = sqrt( rms_y / num_atoms);
                rms_z = sqrt( rms_z / num_atoms);

                rms_t = sqrt( rms_x * rms_x + rms_y * rms_y + rms_z * rms_z );

                printf("\n Chain RMS : %10.6f %10.6f %10.6f\n", rms_x, rms_y, rms_z);
                printf("\n Chain RMS magnitude : %10.6f\n\n", rms_t);

             }
       }
   }
 else
   {
     printf("Bad read!!! returned value: %d\n",good_read);
   }

return 0;
}
