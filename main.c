/******************************************************************************/
/*                                                                            */
/* Inter_vasp is a general purpose program for setting up VASP jobs           */
/* and for converting the results for analysis in Materials Studio.           */
/*                                                                            */
/*                                      Begun June 03, Dave Willock           */
/*                                      Updated Sept 2014                     */
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
#include "debug.h"
#include "own_maths.h"
#include "reader.h"
#undef MAIN

void open_file(FILE **p_file, char *p_filename, char *p_status);

void set_defaults( is_list *p_is, have_list *p_have, need_list *p_need);

void set_up_groups( atom *p_molecule, atom *p_end_molecule, int num_atoms,
                    group_lists *p_groups, int num_groups,
                    double *p_orig_bond1, double *p_orig_bond2,
                    double *p_delta_bond1, double *p_delta_bond2,
                    have_list *p_have, need_list *p_need);

int read_input( FILE *input_fp, char *p_title, char *p_master_input,
                is_list *p_is, have_list *p_have, need_list *p_need, 
                group_lists *p_groups,  labels *p_shell_species,
                char_list *p_image_files, char_list *p_outcar_files, 
                char_list *p_dos_files, perms *p_permute,
                char_list *p_elems_to_match,
                int *p_num_outcars,       int *p_num_to_set,    
                int *p_num_inter,         int *p_num_groups, 
                int *p_linear,            int *p_num_steps, 
                int *p_mode,              int *p_num_shell_species, 
                int *p_compare_modes,     int *p_num_images,     
                int *p_num_atoms_pdos,    int *p_spd, 
                int *p_num_dos_files,     int *p_read_restart,
                int *p_part_dos_list,     int *p_end_min_image,  
                int *p_miller,            int *p_num_miller,   
                int *p_monit_type1,       int *p_monit_type2, 
                int *p_monit_at1,         int *p_monit_at2,   
                int *p_num_per_formula,   int *p_super, 
                int  *p_num_after_switch, int *p_num_elem_to_match, 
                double *p_dos_smear,    double *p_dos_weights, 
                double *p_amplitude,    double *p_switch,
                double *p_temperature,  double *p_min_weight,
                double *p_expansion,    double *p_oshift,
                double *p_vgap,
                char *p_variable_label, char *p_end_input,
                char *p_potcar_input,   char *p_outcar_input,
                char *p_incar_input,    char *p_car_output,
                char *p_arc_output,     char *p_poscar_output, 
                char *p_onetep_output,  char *p_energy_output, 
                char *p_mol_cnt_lab,    char *p_pdb_output,
                char *p_gulp_output,    char *p_mode_compare_input, 
                char *p_doscar_input,   char *p_dos_output,
                char *p_mdtraj_input,   char *p_restart_file,
                char *p_mon_str1,       char *p_mon_str2, 
                char *p_report_file,    char *p_potcar_output,
                char *p_cif_output );

int read_car( FILE *fp, int *p_header_line, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_abc, int *p_been_before, group_lists *p_groups,
              int have_grp, int *p_num_groups,
              int find_fixed, coord_flags *p_fix_flags, int just_count);

int read_cif( FILE *fp, int *p_title_line, 
              atom *p_molecule, int *p_pbc, int *p_num_atoms, 
              double *p_abc, int *p_been_before, int *p_set_labels,
              is_list *p_is, int just_count);

int read_pdb( FILE *fp, int *p_header_line, int *p_title_line, 
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms, 
              double *p_abc, int *p_been_before, int just_count);

int read_gulp( FILE *fp, int *p_header_line, int *p_title_line,
              atom *p_molecule, atom *p_shells, int *p_date_line, int *p_pbc,
              int *p_num_atoms, int *p_num_shells,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_abc, int *p_been_before, int *p_top_bit,
              int *p_num_top_chars, int *p_bottom_bit, int *p_num_bottom_chars,
              int *p_space_group, charge_list *p_spec_charges, int *p_num_species);

int read_poscar( FILE *fp, int *p_title_line,
                 atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
                 int *p_have_labels, int *p_num_of_mols, int *p_num_mol_members, 
                 int *p_mol_number, double *p_latt_vec, double *p_recip_latt_vec,  
                 double *p_abc, int *p_been_before, int *p_is_fract,
		 int *p_is_cart, atom_number *p_types, int *p_num_types, 
                 double *p_scale_factor, coord_flags *p_fix_flags, int just_count);

int read_fdf( FILE *fp, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_latt_vec, double *p_recip_latt_vec,  
              double *p_abc, int *p_been_before, int *p_is_fract,
              int *p_is_cart, int *p_ion_number, int *p_num_types, 
              double *p_scale_factor, int *p_num_free, coord_flags *p_fix_flags);

int read_onetep( FILE *fp, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_latt_vec, double *p_recip_latt_vec,  
              double *p_abc, int *p_been_before, int *p_is_fract,
              int *p_is_cart, int *p_ion_number, int *p_num_types, 
              double *p_scale_factor, int *p_num_free, coord_flags *p_fix_flags,
              char *line_ptrs[], int *p_line_lengths, int *p_num_lines, 
              int just_count );

int read_outcar( FILE *fp, atom *p_molecule, double *p_latt_vec,
                 double *p_recip_latt_vec, double *p_abc, double *p_eigenvals,
                 e_vec *p_eigenvecs, e_vec *p_forces, e_vec *p_chain_forces,
                 int *p_num_atoms, atom_number *p_types, int *p_num_types,
                 int *p_num_modes, int need_freq, int need_force,
                 int need_energy, int *p_have_band, double *p_energy,
                 int need_fermi, double *p_fermi, int *p_num_energies,
                 int just_count );

int read_punch( FILE *fp, atom *p_molecule, atom *p_shells, int *p_num_shells,
                coord_flags *p_fix_flags, int just_count );

int read_incar( FILE *fp, double *p_magmom, int *p_num_magmom );

int read_potcar( FILE *fp, atom *p_molecule,
                 atom_number *p_types, int *p_num_types);

int read_xdatcar(FILE *fp, int *p_num_frames, atom *p_molecule, int num_atoms,
                 coord_flags *p_fix_flags, double *p_latt_vec, double *p_recip_latt_vec,
                 double *p_abc, int *p_is_fract, int *p_is_cart, labels *p_atom_names, 
                 int *p_atom_num, int *p_num_labels );

void generate_neighbours( atom *p_molecule, int num_atoms, 
                          atom_number *p_types, int *p_num_types,
                          int use_pbc,  double *p_recip_latt_vec, double *p_latt_vec,
                          charge_list *p_spec_charges, int set_labels);

void sort_by_elem( atom *p_molecule, int num_atoms, atom_number *p_types, 
                   int num_types );

void sort_by_miller( atom *p_molecule, int num_atoms, atom_number *p_types, int num_types,
                     double *p_latt_vec, double *p_recip_latt, int *p_miller );

void cart_latt_vecs( double *p_abc, double *p_latt_vec, 
                                         double *p_recip_latt_vec);

int compare_strings( char *p_ichar1, char *p_ichar2 );

void string_to_int(char *p_ichar2, int *p_ichar1, int max_position );

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number, int use_mols,
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
                   int is_fract, coord_flags *p_fix_flags, int need_zsort);

void write_cif( FILE *fp, atom *p_molecule, double *p_fract_coords,
                atom_number *p_types, int num_types,
                double *p_latt_vec, double *p_abc, double *p_scale_factor, int num_atoms,
                int *p_title_line, char *p_c_title_line, int pbc, int is_fract,
                coord_flags *p_fix_flags, int need_zsort);

void write_onetep( FILE *fp, atom *p_molecule, double *p_fract_coords,
                   atom_number *p_types, int num_types,
                   double *p_latt_vec, double *p_scale_factor, int num_atoms,
                   int *p_title_line, char *p_c_title_line, int pbc, 
                   int is_fract, int is_onetep, char *file_line_ptrs[],
                   int num_lines);

void rotate_with_flags(atom *p_molecule, double *p_axis, double *p_origin,
                       double theta, int *p_flag_list, int num_atoms);

void inter_atom_vector(atom *p_atom, atom *p_atom2, double *p_vec);

void vec_cross(double *p_A, double *p_B, double *p_cross);

double vec_dot(double *p_A, double *p_B);

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

void move_molecule_with_flags(atom *p_molecule, int *p_flag_list, 
                                           int num_atoms, double *move_vec);

void move_selected_atoms(atom *p_molecule, int num_chosen,
                         int *p_chosen_indices, double *move_vec);

void move_atom(atom *p_atom, double scale, double *move_vec);

void interpolate_all_but_chosen(atom *p_molecule, atom *p_step_molecule, int num_atoms,
                                int *p_chosen, int num_chosen, double scale, double *p_vec );

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
               int *p_super, double *p_recip_latt_vec, double *p_latt_vec);

void min_image( double *x, double *y, double *z,
                           double *p_recip_latt_vec, double *p_latt_vec);

double atom_separation_squared(atom *p_A, atom *p_B, int pbc,
                           double *p_recip_latt_vec, double *p_latt_vec);

/******************************************************************/
/******************************************************************/
/******************************************************************/

void read_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                  int num_atoms_pdos, dos *p_pdos, int *p_spd, int num_ndos_columns, 
                  int num_pdos_columns, int num_dos );

int count_doscar( FILE *fp, dos *p_ndos, int need_pdos, int *p_part_dos_list, 
                 int num_atoms_pdos, dos *p_pdos, int *p_spd, int *p_num_dos_points,
                 int *p_num_ndos_columns, int just_count );

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
                      double *p_recip_latt, double *p_abc, int *p_miller, int num_miller,
                      double *p_new_latt, double *p_new_recip_latt, 
                      double *p_new_abc, atom *p_slab_mol, int *p_num_slab_atoms);

void define_cluster(atom *p_molecule, int *p_num_atoms, double *p_latt, 
                    double *p_recip_latt, double *p_abc, double *p_miller, int num_miller,
                    atom *p_cluster_mol, int *p_num_cluster_atoms, double radius);

void define_hemi(atom *p_molecule, int *p_num_atoms, double *p_latt, 
                 double *p_recip_latt, double *p_abc, double radius,
                 double *p_origin, atom *p_cluster_mol, int *p_num_cluster_atoms);

void write_csv(FILE *fp, char *p_title_x, char *p_title_y, 
               char *p_title_z,
               double *p_x, double *p_y, double *p_z,
               int have_y, int have_z, int num);

int read_report( FILE *fp, int *p_num_frames, int just_count,
                 double *p_etot, double *p_epot, double *p_ekin, double *p_tsim, double *p_tinst);

int permutation(atom *p_master, int *p_num_atoms, double *p_latt_vec, double *p_recip_latt_vec,
                int pbc, double *p_abc, perms *p_permute, double *p_magmom, int num_magmom );

void count_oxygen_species(atom *p_molecule, int num_atoms, int *oxy_list_ptrs[],
                          atom_number *p_num_oxy_species, int num_oxy_types, int just_count);

void find_oh_bonds( atom *p_molecule, int num_atoms, int *oxy_list_ptrs[],
                    atom_number *p_num_oxy_species, int num_oxy_types, 
                    int use_pbc, double *p_latt_vec, double *p_recip_latt_vec);

void block_print_int(int *p_array, int num, int num_per_row);

/******************************************************************/
/*** Start of MAIN ************************************************/
/******************************************************************/

int main(argc, argv)
  int   argc;
  char  **argv;
{

/**** These arrays should all be replaced by Malloced pointers ****/
/**** the molecule[] case has already been done so just need   ****/
/**** to extend the approach to the other arrays.              ****/

  atom *p_molecule, *p_atom, *p_atom2, *p_at1, *p_at2;
  atom *p_end_molecule, *p_end_atom, *p_end_atom2;
  atom *p_trans_atom, *p_rec_atom;

  atom *p_step_molecule;
  atom *p_step_molecule_last;
  atom *p_step_atom, *p_step_atom2, *p_step_cnt_atom;
  atom *p_step_atom_last;

  atom slab_mol[MAXATOMS], cluster[MAXATOMS];
/* atom step_molecule[MAXATOMS]; */
  atom chosen_one[MAXATOMS];
  atom chosen_one_end[MAXATOMS];
  atom *p_shells, *p_shll;
  atom *p_end_shells;

  atom_number types[MAXATOMS], end_types[MAXATOMS];

  bonds old_bonds[MAXBONDS], new_bonds[MAXBONDS];
  int num_old_bonds, num_new_bonds;

  e_vec eigenvecs[MAXMODES], forces, chain_forces;
  e_vec new_eigenvecs[MAXMODES];

  char_list image_files[MAXIMAGES];
  char_list outcar_files[MAXOUTCARS];
  char_list dos_files[MAX_DOS_FILES];
  char_list elems_to_match[MAX_MATCH_LIST];

  double eigenvalues[MAXMODES];
  double magmom[MAXATOMS];
  double new_eigenvalues[MAXMODES];
  double dist, dmin, diffx, diffy, diffz, max_freq_shift;
  double rms_x, rms_y, rms_z, rms_t;
  double *p_energy_vasp;
  double fermi;
  double *p_d_monit, *p_this_d_monit, *p_fdum;
  double *p_etot,*p_epot,*p_ekin,*p_tsim,*p_tinst;
  double data[4];
  double vgap, oshift[3];

  int num_atoms, num_shells, good_read, num_in_group, num_groups, num_modes;
  int num_energies, num_elem_to_match;
  int num_new_atoms, num_new_types, num_new_modes;
  int num_chosen_atoms, num_free, match, imatch, imat, isign;
  int end_num_atoms, end_num_shells, start_frame, imode, jmode;
  int done_start_frame = FALSE;
  int header_line[LINESIZ], title_line[LINESIZ], date_line[LINESIZ];
  int end_header_line[LINESIZ], end_title_line[LINESIZ], end_date_line[LINESIZ];
  int pbc, num_of_mols, num_mol_members[MAXMOL], mol_number[MAXATOMS]; 
  int end_num_of_mols, end_num_mol_members[MAXMOL], end_mol_number[MAXATOMS]; 
  int map_atoms[MAXATOMS];
  int been_before, grp_cnt, grp_mem_ind, igrp, jgrp, imol;
  int find_molecules, mol_cnt_ind, mol_ind;
  int iii,iloop, jloop, iatom, jatom, iatom2; 
  int ishell, icomp, icoord, iletter, iline;
  int idos, index, itype, ispecies;
  int ineigh, neigh_index, num_types, end_num_types;
  int count_monit, ineigh2, still_there;
  int start_mol, num_in_this_mol;
  int ion_number[MAXTYPES];
  int num_steps, linear;
  int starting_late, now_close, num_steps_left;
  int iframe, num_frames, num_rep_frames;
  int flags[MAXATOMS];
  int num_images, num_outcars;
  int end_min_image;
  int compare_modes;
  int num_lines, line_lengths[MAX_LINES];
  int just_count, use_mols;
  int molecule_malloced;
  int leave, ivec[3];
  int set_labels, clash;
/*** variables needed for count_oxy stuff **/
  atom_number  num_oxy_species[4]; 
  int num_oxy_types=4;
  
/*** An array of pointers for the lists of oxygen atom indices of each type ***/
  int *oxy_list_ptrs[4];
/*** variables needed for count_oxy stuff **/

  coord_flags *p_fix_flags;

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
  double d, d2, this_d2, dx, dy, dz;
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
  char mol_cnt_lab[7];
  char elem_rec[7];
  char mon_str1[7];
  char mon_str2[7];
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
  char report_file[60];
  char dos_output[60];
  char poscar_output[60];
  char onetep_output[60];
  char step_output[60];
  char car_output[60];
  char cif_output[60];
  char arc_output[60];
  char pdb_output[60];
  char gulp_output[60];
  char mode_freqs_output[60];
  char frame_output[60];
  char frame_output_pdb[60];
  char mode_compare_input[60];
  char command[180];
  char energy_output[60];
  char potcar_output[60];
  char *file_line_ptrs[MAX_LINES];

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
  FILE *fp_cif_input;
  FILE *fp_outcar;
  FILE *fp_punch_file;
  FILE *fp_report;
  FILE *fp_car_output;
  FILE *fp_arc_output;
  FILE *fp_pdb_output;
  FILE *fp_pdb_input;
  FILE *fp_vasp_input;
  FILE *fp_siesta_input;
  FILE *fp_onetep_input;
  FILE *fp_vasp_output;
  FILE *fp_cif_output;
  FILE *fp_vasp_energy;
  FILE *fp_vasp_report;
  FILE *fp_onetep_output;
  FILE *fp_doscar_input;
  FILE *fp_dos_output;
  FILE *fp_monit;
  FILE *fp_ecsv;
  FILE *fp_input;
  FILE *fp_mode_freqs;

double *p_fract_coords, *p_end_fract_coords;
double *p_this_frac;
double *p_inter_vec, *p_vec;
double temperature,  min_weight;
double scale_factor, step, cstep, step_inc;
double axis[3], origin[3];
double dos_smear, dos_weights[MAX_DOS_FILES];
double dos_norm, hemi_rad, clus_rad;
double expansion[3], norm[3];
double new_miller[MAX_MILLER3];

int num_to_set, num_per_formula;
int top_bit[500], bottom_bit[500];
int num_top_chars, num_bottom_chars, num_species;
int space_group, state;
int monit_type1, monit_type2;
int num_mon_list1, num_mon_list2;
int *p_this_mon_list1, *p_this_mon_list2;
int num_dos_files, num_perms;
int num_dos_points, spd[5];
int part_dos_list[MAXATOMS], num_atoms_pdos;
int lowpdos, hipdos, this_pdos, nrang;
int num_shell_species;
int atom_transfered, atom_losing, atom_receiving;
int which_mode;
int num_inter;
int num, num_magmom;
int check;
int miller[MAX_MILLER3], num_miller, look_for_layers;
int imiller, h_miller, k_miller, l_miller;
int chosen_indices[MAXATOMS];
int super[3], surf_super[3];
int any_rec, matched_neigh;
int read_restart;
int num_ndos_columns, num_pdos_columns;
int num_slab_atoms, num_cluster_atoms;
int atom_num[MAXTYPES];
int num_labels;
int monit_at1, monit_at2;
int *p_mon_list1, *p_mon_list2;
int num_after_switch;

group_lists groups;

charge_list spec_charges[MAXATOMS];

labels shell_species[10];
labels ion_labels[MAXTYPES];
labels atom_names[MAXTYPES];

dos *p_ndos, *p_pdos;
dos *p_temp_ndos, *p_temp_pdos;

perms permute;

/**** Structures that bring together control variables ************/
have_list have;
is_list is;
need_list need;

/*********************************************************/
/*** Start of executable code, set some default values ***/
/*********************************************************/

printf("Starting inter_vasp version 3.1.0\n");
printf("Last Built : 28th April 2019\n");
printf("Larger array mallocing begun.\n\n");

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
use_mols = TRUE;
which_mode= -1;
amplitude = 0.5;
dos_smear = 0.2;
num_magmom=-1;

pi = acos(-1.0);
printf("Pi is %10.6f\n", pi);

monit_type1 = -1; 
monit_type2 = -1; 
num_shell_species=-1;

spd[0]= TRUE;
spd[1]= TRUE;
spd[2]= TRUE;
spd[3]= TRUE;

any_rec = FALSE; compare_modes = FALSE; read_restart= FALSE;

num_miller = -1;

num_to_set = 0;
num_per_formula = 0;
traj_switch = 1.5;
num_after_switch= 5;
num_shells=-1;
num_dos_files=0;
num_images=-1;

super[0] = 1;
super[1] = 1;
super[2] = 1;

expansion[0] = 0.10;
expansion[1] = 10.0;
expansion[2] = 1.00;

surf_super[0] = 1;
surf_super[1] = 1;
surf_super[2] = 1;

temperature = 0.0;
min_weight  = 0.0;

num_groups=-1;

permute.centre=-1;
permute.check=FALSE;
permute.mind=-1;
permute.maxd=-1;
permute.debug=FALSE;

molecule_malloced=FALSE;

find_molecules=TRUE;

/*******************************************************************/
/*********** read input file (name input by user on the ************/
/*********** command line).                             ************/
/*******************************************************************/

read_input( fp_control_file, &c_title_line[0], &master_input[0],
            &is, &have, &need, 
            &groups,           &shell_species[0], 
            &image_files[0],   &outcar_files[0],
            &dos_files[0],     &permute,
            &elems_to_match[0],
            &num_outcars,      &num_to_set, 
            &num_inter,        &num_groups,
            &linear,           &num_steps, 
            &which_mode,       &num_shell_species, 
            &compare_modes,    &num_images,
            &num_atoms_pdos,   &spd[0],
            &num_dos_files,    &read_restart, 
            &part_dos_list[0], &end_min_image, 
            &miller[0],        &num_miller,
            &monit_type1,      &monit_type2, 
            &monit_at1,        &monit_at2,
            &num_per_formula,  &super[0],
            &num_after_switch, &num_elem_to_match,
            &dos_smear,        &dos_weights[0], 
            &amplitude,        &traj_switch, 
            &temperature,      &min_weight,
            &expansion[0],     &oshift[0],
            &vgap,
            &variable_label[0],&end_input[0], 
            &potcar_input[0],  &outcar_input[0],
            &incar_input[0],   &car_output[0], 
            &arc_output[0],    &poscar_output[0], 
            &onetep_output[0], &energy_output[0],
            &mol_cnt_lab[0],   &pdb_output[0], 
            &gulp_output[0],   &mode_compare_input[0],
            &doscar_input[0],  &dos_output[0], 
            &mdtraj_input[0],  &restart_file[0],
            &mon_str1[0],      &mon_str2[0], 
            &report_file[0],   &potcar_output[0],
            &cif_output[0] );

fclose(fp_control_file);

/*******************************************************************/
/** Check what the user has asked for and try to catch problems ****/
/** such as conflicting directives.                             ****/
/*******************************************************************/

printf("c_title_line >>%s<<\n", c_title_line);
printf("master_input >>%s<<\n", master_input);

if (need.md_run)
  {
     printf("Will perform our own MD run using VASP forces\n");
  }

if (need.multi_dos)
  {
    printf("Found %d DOSCAR files will expect to amalgamate them\n", num_dos_files);

    for (iloop=0; iloop <= num_dos_files; iloop++)
      {
         printf("%s with weight %10.6f\n", dos_files[iloop].name, dos_weights[iloop]);
      }
  }

need.interpolate = (is.end_gulp || is.end_car || is.end_vasp) && num_images < 0;
need.react_coord = (is.end_gulp || is.end_car || is.end_vasp) && num_images > 0;

if (num_images > 0 && !need.react_coord)
  {
    printf("ERROR: images supplied but no end point reference structure\n");
    exit(0);
  }
else if (num_images==0)
  {
    printf("ERROR: Cannot have zero images in a react coords run\n");
    exit(0);
  }

printf("Job Title        : %s\n\n", c_title_line);
printf("Master file      : %s\n", master_input);
if (read_restart)
   {
      printf("Will take atom co-ordinates from the restart (CONTCAR) file: %s\n", restart_file);
   }
if (is.gulp)
   {
     printf("                           this is a gulp structure file\n");
   }
if (is.car)
   {
     printf("                           this is a car structure file\n");
   }
if (is.cif)
   {
     printf("                           this is a cif structure file\n");
   }
if (is.pdb)
   {
     printf("                           this is a pdb structure file\n");
   }
if (have.out)
   {
     printf("                           this is a VASP OUTCAR structure file\n");
     strcpy(outcar_input, master_input);
     printf(" The POTCAR file is : %s\n", potcar_input);
   }
if (is.vasp)
   {
     printf("                           this is a VASP POSCAR structure file\n");
     if (have.potcar) printf(" The POTCAR file is : %s\n", potcar_input);
   }
if (is.siesta)
   {
     printf("                           this is a SIESTA structure file\n");
     printf(" The fdf file is : %s\n", master_input);
   }
if (is.onetep)
   {
     printf("                           this is a ONETEP structure file\n");
     printf(" The dat file is : %s\n", master_input);
   }

if (is.end_gulp)
   {
     printf("End point file will be read from %s, this is a gulp structure file\n", end_input);
   }
if (is.end_car)
   {
     printf("End point file will be read from %s, this is a car structure file\n", end_input);
   }
if (is.end_vasp)
   {
     printf("End point file will be read from %s, this is a VASP POSCAR structure file\n", end_input);
   }


if ( need.interpolate ) 
  {
    if (end_min_image)
       {
         printf("Will use minimum image of end point atoms to corresponding start point for interpolation\n");
       }
    printf("Will provide %d interpolated structures\n", num_inter);

    if (is.end_car) printf("based on master and endpoint car file co-ordinates\n");

    if (need.shift)
      {
         printf("Will shift all structures to elliminate drift during interpolation\n");
      }
    if (need.arc)
      {
        printf("Will generate an arc file of the elastic band set of structures named : %s\n",arc_output);

        open_file( &fp_arc_output, arc_output, "w");
      }

    if (have.grp && !need.late) 
      {
             printf("Have %d groups :\n", num_groups);

                  if (groups.group_type1 == CENTRE_TYPE)
                    {
                       printf("Will interpolate GRP1 group centred on >>%s<<\n", groups.cnt_lab1);
                    }
                  else if (groups.group_type1 == ANGLE_TYPE)
                    {
                       printf("Will rotate the GRP1 atoms around the axis defined by the line\n");
                       printf("from atom >>%s<< to atom >>%s<<\n", groups.axis1_lab, groups.axis2_lab);
                    }
                  else
                    {
                       printf("ERROR: No type set for group interpolation\n");
                       exit(0);
                    }

             if (num_groups >= 1)
               {
                  if (groups.group_type2 == CENTRE_TYPE)
                    {
                       printf("Will interpolate GRP2 group centred on >>%s<<\n", groups.cnt_lab2);
                    }
                  else if (groups.group_type2 == ANGLE_TYPE)
                    {
                       printf("Will rotate the GRP2 atoms around the axis defined by the line\n");
                       printf("from atom >>%s<< to atom >>%s<<\n", groups.axis1_lab, groups.axis2_lab);
                    }
               }

        if (!is.car) 
           {
              printf("ERROR : group centre defined but structure not given as car file\n");
              exit(0);
           }
      }

    if (have.mol)
      {
        if ( need.late )
          {

        if (have.grp)
          {
            printf("Working with a molecule coming to a surface that also contains a group.............\n");
          }
        printf("Will interpolate molecule containing");
        printf(" >>%s<< to look for late transition state.\n",mol_cnt_lab);
        printf("Using %10.6f as the multiplier of the final bondlength at which the", 
                                                                          traj_switch);
         
        printf("After switch will do a further %d images\n", num_after_switch);
        printf(" transfer occurs\n");
          }
        else
          {
        printf("Will interpolate molecule containing >>%s<< as rigid body.\n"
                      ,mol_cnt_lab);
            if (need.morph)
              {
                printf("Internal geometry of molecule will be altered smoothly to end point\n");
              }
          }

        if (!is.car) 
           {
              printf("ERROR : molecule centre defined but structure not given as car file\n");
              exit(0);
           }
      }

    if (need.match)
      {
        printf("Atoms of elements: ");
        for ( iii=0; iii<= num_elem_to_match; iii++) printf("%s ", elems_to_match[iii].name );
        printf(" will be matched with their closest versions between start and end.\n");
      }

  }
else
  {
    printf("This is not an interpolation run\n");
  }

if (need.zsort)
  {
    printf("Will sort POSCAR atoms into z-order\n");
  }
if (need.oshift)
  {
    printf("Will shift all atoms by %10.6f  %10.6f  %10.6f \n", oshift[0], oshift[1], oshift[2]);
  }
if (need.vgap)
  {
    printf("Will add a vacuum gap of %10.6f Angstroms.\n", vgap);
  }

if (need.react_coord)
  {
    printf("Will provide reaction co-ordinate values for images with respect the master and end points\n");

  }

if (have.out) printf("OUTCAR file is   : %s\n", outcar_input);
if (have.incar) 
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

if (need.force)
  {
     if (have.out)
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
if (need.freq)
  {
     if (have.out)
       {
         printf("Will analyse OUTCAR file for vibrational frequencies\n");
         if (which_mode > 0)
              printf("Animation files for mode %d only requested\n", 
                            which_mode);
       }
     else if (is.siesta)
       {
         printf(" Will look for Eigenvalues and vectors in siesta '.vectors' file\n");
       }
     else 
       {
         printf("ERROR : Frequencies requested but no OUTCAR supplied\n");
       }
     if (need.arc)
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

/***************************************************************/
/*** Report on what will be done while writing car file ********/
/***************************************************************/

if (need.car)
	{
	  printf("Will produce a car file called %s\n", car_output);
          printf("Will generate super cell of %d %d %d extent\n",
                            super[0], super[1], super[2]);

          if (need.monit) 
            {
               printf("Will monitor distance between atoms %d and %d\n", monit_at1+1, monit_at2+1);
            }                                                  
	}
if ( need.monit && !(need.car || need.mdtraj) )
   {
      printf("ERROR: Currently only monitor atom distances when need car_file or need.mdtraj (XDATCAR) also set\n");
      exit(0);
   }
      
if (need.expansion)
	{
	  printf("Will perform a unit cell expansion\n");
          printf("Using expansion parameter of %4.2f\n",
                            expansion[0]);
        }
if (need.pdb)
	{
	  printf("Will produce a pdb file named >>%s<<\n", pdb_output);
          printf("pdb files will show super cell of %d %d %d extent\n",
                            super[0], super[1], super[2]);
	}
if (need.gulp)
	{
	  printf("Will produce a gulp file named >>%s<<\n", gulp_output);
          printf("gulp files will show super cell of %d %d %d extent\n",
                            super[0], super[1], super[2]);
          if (need.shells)
            {
               printf("Will include shells for species: ");
               for (iloop=0; iloop<=num_shell_species; iloop++)
                 {
                   printf("%s ", shell_species[iloop].label);
                 }
               printf("\n");
            }
	}
if (need.onetep)
	{
	  printf("Will produce a ONETEP file called %s\n", onetep_output);
	}
if (need.poscar)
	{
	  printf("Will produce a POSCAR file called %s\n", poscar_output);
	}
if (need.mdtraj)
        {
          printf("Will convert the trajectory file %s to format ", mdtraj_input);
          if (need.arc)
            {
               printf(" arc file.\n");
               printf("Producing file with stem : %s\n", arc_output); 
            }
          else if (need.pdb)
            {
               printf(" pdb file.\n");
               printf("Producing file with stem : %s\n", pdb_output); 
            }
          else
            {
              printf("\nERROR : No format type defined for trajectory output\n");
              exit(0);
            }
          if (need.monit)
            {
              printf("Monitoring option for MD trajectory analysis also enabled\n");
              printf("DEBUG>> monit_type ints are %d and %d\n", monit_type1, monit_type2);

/*** monit_type will have been set so that 1 means a particular atom number, 2 an element and 3 a label **/

              if (monit_type1 == 1)
                  printf("Monitor type 1 is atom number and atom %d has been requested\n", monit_at1+1);
              else if (monit_type1 == 2)
                  printf("Monitor type 1 is an element %s has been requested\n", mon_str1);
              else if (monit_type1 == 3)
                  printf("Monitor type 1 is an atom label %s has been requested\n", mon_str1);

              if (monit_type2 == 1)
                  printf("Monitor type 2 is atom number and atom %d has been requested\n", monit_at2+1);
              else if (monit_type2 == 2)
                  printf("Monitor type 2 is an element %s has been requested\n", mon_str2);
              else if (monit_type2 == 3)
                  printf("Monitor type 2 is an atom label %s has been requested\n", mon_str2);
            }
        }


/*******************************************************************/
/** MD analysis options.                                        ****/
/*******************************************************************/

if (have.report)
  {
     printf("Have a REPORT file to analyse....\n");
     open_file( &fp_report, report_file, "r");

     just_count=TRUE;
     num_rep_frames=-1;
     read_report( fp_report, &num_rep_frames, just_count,
                  p_etot, p_epot,p_ekin,p_tsim,p_tinst);

     fclose(fp_report);
     printf("Have had a first look...mallocing\n");
     p_etot=(double*)malloc(num_rep_frames*sizeof(double));
     p_epot=(double*)malloc(num_rep_frames*sizeof(double));
     p_ekin=(double*)malloc(num_rep_frames*sizeof(double));
     p_tsim=(double*)malloc(num_rep_frames*sizeof(double));
     p_tinst=(double*)malloc(num_rep_frames*sizeof(double));
     printf("In main with %d frames from report\n", num_rep_frames);

     open_file( &fp_report, report_file, "r");

     just_count=FALSE;
     read_report( fp_report, &num_rep_frames, just_count,
                  p_etot, p_epot,p_ekin,p_tsim,p_tinst);
     fclose(fp_report);
     printf("Have had a second look............\n");

     sprintf(filename, "report_sum.csv");
     open_file( &fp_ecsv, filename, "w");

     sprintf(title_x, "etot");
     sprintf(title_y, "epot");
     sprintf(title_z, "ekin");

     write_csv(fp_ecsv, title_x, title_y, title_z,
               p_etot, p_epot, p_ekin,
               TRUE, TRUE, num_rep_frames-1);

     fclose(fp_ecsv);
     sprintf(filename, "report_sum_temp.csv");
     open_file( &fp_ecsv, filename, "w");
     sprintf(title_x, "tsim");
     sprintf(title_y, "tinst");
     sprintf(title_z, " ");
     write_csv(fp_ecsv, title_x, title_y, title_z,
               p_tsim, p_tinst, p_fdum,
               TRUE,FALSE, num_rep_frames-1);
    fclose(fp_ecsv);
    exit(0);
  }


/*******************************************************************/
/****** Try and drive VASP MD using our own integration routines ***/
/****** needs development work                                   ***/
/*******************************************************************/

if (need.md_run)
  {
/**** Needs changing for malloced version ***/
     printf("ERROR: Malloced version cannot currently cope with MD OUTCARs\n");
     exit(0);
//     timestep = 0.001;
//     set_temp = 300;
//     open_file( &fp_outcar, master_input, "r");

/**** Read in current forces and energy from OUTCAR ****/

//     need.fermi=FALSE;
//     good_read = read_outcar( fp_outcar, &molecule[0], &latt_vec[0],
//                              &recip_latt_vec[0], &abc[0], &eigenvalues[0],
//                              &eigenvecs[0], &forces, &chain_forces,
//                              &num_atoms, &types[0], &num_types, &num_modes, need_freq,
//                              need_force, need_energy, &(have.band), p_energy_vasp, need_fermi,
//                              &fermi, &num_energies, just_count);

//     pbc = TRUE;   /* This should be set by read_outcar */

//     printf("Back from read_outcar with %d atoms\n\n", num_atoms);

//     for (iatom= 0; iatom <= num_atoms; iatom++)
//       {
//         printf("%d %s %10.6f  %10.6f  %10.6f  %10.6f force: %10.6f  %10.6f  %10.6f\n",
//                                 iatom+1, molecule[iatom].label,
//                                 molecule[iatom].mass,
//                                 molecule[iatom].x, 
//                                 molecule[iatom].y, 
//                                 molecule[iatom].z, 
//                                 forces.dx[iatom],
//                                 forces.dy[iatom],
//                                 forces.dz[iatom]);
//       }

/**** Write arc file ****/

//     open_file( &fp_arc_output, "md_movie.arc", "w");

//     sprintf(c_title_line,
//             "Inter_vasp generated arc file.");

//     scale_factor = 1.0;

//     start_frame=TRUE;

//     write_car( fp_arc_output, &header_line[0],
//                &title_line[0], &c_title_line[0],
//                &date_line[0], &molecule[0],
//                &mol_number[0], pbc, &abc[0],
//                num_atoms+1, scale_factor, start_frame,
//                &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0],
//                &magmom[0], num_magmom);

//     start_frame=FALSE;

/**** Set all velocities to the mean for temperature ( in 2D ) ****/
/**** using KE = RT and 1/2 m v^2                              ****/

//       ke = BOLTZ * set_temp;

//     for (iatom= 0; iatom <= num_atoms; iatom++)
//       {
//       speed = sqrt( 2.0 * ke / molecule[iatom].mass);

//       if ( iatom == 0 )
//               printf("Initial speed set to %10.6f\n", speed);

/****  ( 1.0 - 2.0*drand48() ) will be distributed +/- 1 **********/

//       vec[0] =  1.0 - 2.0*drand48();
//       vec[1] =  1.0 - 2.0*drand48();
//       vec[2] =  1.0 - 2.0*drand48();

//       unit_vector(&vec[0]);

//       molecule[iatom].vx = vec[0] * speed;
//       molecule[iatom].vy = vec[1] * speed;
//       molecule[iatom].vz = vec[2] * speed;
//
//     }

/**** Loop over number of MD runs ****/

//    for (iatom=0; iatom <= num_atoms; iatom++)
//       {
//
///*** a(t) = f(t) / m *****/
 
/*** Include conversion between eV and internal units for forces ***/

//          molecule[iatom].ax = forces.dx[iatom]/(INTUNITS_TO_EV*molecule[iatom].mass);
//          molecule[iatom].ay = forces.dy[iatom]/(INTUNITS_TO_EV*molecule[iatom].mass);
//          molecule[iatom].az = forces.dz[iatom]/(INTUNITS_TO_EV*molecule[iatom].mass);

/*** v(t+dt/2) = v(t-dt/2) + a(t) ***/
/*** for first step put in half iteration to offset velocity from positions ****/
//          if (iloop==0)
//            {
//              molecule[iatom].vx += molecule[iatom].ax * 0.5*timestep;
//              molecule[iatom].vy += molecule[iatom].ay * 0.5*timestep;
//              molecule[iatom].vz += molecule[iatom].az * 0.5*timestep;
//            }
//          else
//            {
//              molecule[iatom].vx += molecule[iatom].ax * timestep;
//              molecule[iatom].vy += molecule[iatom].ay * timestep;
//              molecule[iatom].vz += molecule[iatom].az * timestep;
//
//            }

/*** r(t+dt) = r(t) + v(t+dt/2) t *******/
//          molecule[iatom].x += molecule[iatom].vx * timestep;
//          molecule[iatom].y += molecule[iatom].vy * timestep;
//          molecule[iatom].z += molecule[iatom].vz * timestep;


//       }
//     printf("\n");
//
//     for (iatom= 0; iatom <= num_atoms; iatom++)
//       {
//         printf("%d %s %10.6f  %10.6f  %10.6f  %10.6f acc: %10.6f  %10.6f  %10.6f\n",
//                                 iatom+1, molecule[iatom].label,
//                                 molecule[iatom].mass,
//                                 molecule[iatom].x,
//                                 molecule[iatom].y,
//                                 molecule[iatom].z,
//                                 molecule[iatom].ax,
//                                 molecule[iatom].ay,
//                                 molecule[iatom].az);
//       }


/**** Write the new co-ordinates to a POSCAR file ******/

//    open_file( &fp_vasp_output, "POSCAR", "w");
//
//    scale_factor = 1.0;
//    write_poscar( fp_vasp_output, &molecule[0],  &fract_coords[0],
//                  &types[0], num_types, &latt_vec[0], &scale_factor, num_atoms,
//                  &title_line[0], &c_title_line[0], pbc, FALSE, &fix_flags[0], 
//                  need.zsort);

//    fclose(fp_vasp_output);

/**** Run the vasp job to get next set of forces ***/

//    fflush(stdout);

//    sprintf(command, "./vasp_script");

//    system(command);

/**** Write arc file ****/

//    write_car( fp_arc_output, &header_line[0],
//               &title_line[0], &c_title_line[0],
//               &date_line[0], &molecule[0],
//               &mol_number[0], pbc, &abc[0],
//               num_atoms+1, scale_factor, start_frame,
//               &super[0], &latt_vec[0], &recip_latt_vec[0], &fix_flags[0], 
//               &magmom[0], num_magmom);

/**** End loop over number of MD runs here ****/

//    fclose(fp_arc_output); 
//    exit(0);
  }

/*** Report variables for a permutation run ****/
if (need.permute)
  {
    printf("This is a run to permute atoms between sites:\n");
    printf("Expecting marker element : %s\n", permute.elem);
    printf("Expecting number to leave: %d\n", permute.num);
    if ( !permute.mind )
       printf("Expecting min_distance   : %10.6f\n", permute.mind);
    if ( !permute.maxd )
       printf("Expecting max_distance   : %10.6f\n", permute.maxd);
    if ( permute.centre != -1 )
      {
        printf("This will consider the distances of atoms to selected centre with index : %d\n",
                                                                        permute.centre-1 );
      }
    if ( permute.check )
      {
        printf("This is a run to check the number of structures that would be generated.\n");
      }
    else
      {
        printf("This run will generate the structures.\n");
      }
  }

/******************************************************************/
/*********** read the master structure file            ************/
/*********** The just_count option allows the routines ************/
/*********** to count the number of atoms present so   ************/
/*********** that memory allocation can be used to     ************/
/*********** dimension the "p_molecule" array.         ************/
/******************************************************************/

printf("\nReading master file.....\n\n");

if (is.car)
  {
   printf("Master is a car file....\n");
   open_file( &fp_car_input, master_input, "r");


/**** First of all just count the number of atoms in the car file ****/

    if (!molecule_malloced) 
      {
        just_count = TRUE;
      }
    else
      {
        just_count = FALSE;
      }

    been_before= FALSE;
    good_read=  read_car( fp_car_input, &header_line[0], &title_line[0],
                          p_molecule, &date_line[0], &pbc, &num_atoms,
                          &num_of_mols, &num_mol_members[0], &mol_number[0],
                          &abc[0], &been_before, &groups, have.grp, 
                          &num_groups, TRUE, p_fix_flags, just_count);

    printf("Been to car file and counted %d atoms\n", num_atoms);

    fclose(fp_car_input);

/**** Reserve memory space for the molecule array ****/
/**** Note that num_atoms is set to the actual      **/
/**** number of atoms.                              **/

    if (just_count)
      {
        printf("mallocing molecule for %d atoms, molecule and fix_flags\n", num_atoms);
        p_molecule =(atom*)malloc(num_atoms*sizeof(atom));
        p_fix_flags=(coord_flags*)malloc(num_atoms*sizeof(coord_flags));
        molecule_malloced=TRUE;

/**** Read in the actual atom data now that space has been reserved **/
/**** For master input check for fixing flags in group labels ****/

        open_file( &fp_car_input, master_input, "r");
        just_count = FALSE;
        been_before= FALSE;
        good_read=  read_car( fp_car_input, &header_line[0], &title_line[0],
                              p_molecule, &date_line[0], &pbc, &num_atoms,
                              &num_of_mols, &num_mol_members[0], &mol_number[0],
                              &abc[0], &been_before, &groups, have.grp, 
                              &num_groups, TRUE, p_fix_flags, just_count);

        fclose(fp_car_input);

        printf("Back from read_car for the second time with %d atoms..\n\n", num_atoms);
     }

/*    if (read_restart) */
/*      { */
/*        printf("ERROR: Restart directive given but master file is a car format file\n"); */
/*        printf("       Restart is only intended for copying fixed atom flags over \n"); */
/*        printf("       from POTCAR type files to CONTCAR type structures. \n"); */
/*        exit(0); */
/*      } */

    scale_factor=1.0;

    is.cart = TRUE;
    is.fract= FALSE;
    set_labels = FALSE;
  }
else if (is.pdb)
  {
   printf("Master is a pdb file....\n");
   open_file( &fp_pdb_input, master_input, "r");

/**** First of all just count the number of atoms in the car file ****/

    if (!molecule_malloced) 
      {
        just_count = TRUE;
      }
    else
      {
        just_count = FALSE;
      }

    been_before= FALSE;

    good_read=  read_pdb( fp_pdb_input, &header_line[0], &title_line[0],
                          p_molecule, &date_line[0], &pbc, &num_atoms,
                          &abc[0], &been_before, just_count);

    printf("Been to pdb file and counted %d atoms\n", num_atoms);

    fclose(fp_pdb_input);

/**** Reserve memory space for the molecule array ****/
/**** Note that num_atoms is set to the actual      **/
/**** number of atoms.                              **/

    if (just_count)
      {
        printf("mallocing molecule for %d atoms, molecule and fix_flags\n", num_atoms);
        p_molecule =(atom*)malloc(num_atoms*sizeof(atom));
        p_fix_flags=(coord_flags*)malloc(num_atoms*sizeof(coord_flags));
        molecule_malloced=TRUE;

/**** Read in the actual atom data now that space has been reserved **/
/**** For master input check for fixing flags in group labels ****/

        open_file( &fp_pdb_input, master_input, "r");
        just_count = FALSE;
        been_before= FALSE;
        good_read=  read_pdb( fp_pdb_input, &header_line[0], &title_line[0],
                              p_molecule, &date_line[0], &pbc, &num_atoms,
                              &abc[0], &been_before, just_count);

        fclose(fp_pdb_input);

        printf("Back from read_pdb for the second time with %d atoms..\n\n", num_atoms);
     }

/*    if (read_restart) */
/*      { */
/*        printf("ERROR: Restart directive given but master file is a car format file\n"); */
/*        printf("       Restart is only intended for copying fixed atom flags over \n"); */
/*        printf("       from POTCAR type files to CONTCAR type structures. \n"); */
/*        exit(0); */
/*      } */

    scale_factor=1.0;
    num_mol_members[0]=num_atoms;

    is.cart = TRUE;
    is.fract= FALSE;
    set_labels = TRUE;
  }
else if (is.cif)
  {
   printf("Master is a cif file.. not sure what to do......\n");
   open_file( &fp_cif_input, master_input, "r");


/**** First of all just count the number of atoms in the cif file ****/

    if (!molecule_malloced) 
      {
        just_count = TRUE;
      }
    else
      {
        just_count = FALSE;
      }

    good_read = read_cif( fp_cif_input, &title_line[0], 
                          p_molecule, &pbc, &num_atoms, 
                          &abc, &been_before, &set_labels,
                          &is, just_count);

    num_mol_members[0]=num_atoms;

    printf("read_cif found %d atoms...\n", num_atoms);

    fclose(fp_cif_input);

/**** Reserve memory space for the molecule array ****/
/**** Note that num_atoms is set to the actual      **/
/**** number of atoms.                              **/

    if (just_count)
      {
        printf("mallocing molecule for %d atoms, molecule and fix_flags\n", num_atoms);
        p_molecule =(atom*)malloc(num_atoms*sizeof(atom));
        p_fix_flags=(coord_flags*)malloc(num_atoms*sizeof(coord_flags));
        molecule_malloced=TRUE;

/**** Read in the actual atom data now that space has been reserved **/
/**** For master input check for fixing flags in group labels ****/

        open_file( &fp_cif_input, master_input, "r");
        just_count = FALSE;
        been_before= FALSE;
        good_read = read_cif( fp_cif_input, &title_line[0], 
                              p_molecule, &pbc, &num_atoms, 
                              &abc, &been_before, &set_labels,
                              &is, just_count);

        fclose(fp_cif_input);

        printf("Back from read_cif for the second time with %d atoms..\n\n", num_atoms);
     }

/*    if (read_restart) */
/*      { */
/*        printf("ERROR: Restart directive given but master file is a car format file\n"); */
/*        printf("       Restart is only intended for copying fixed atom flags over \n"); */
/*        printf("       from POTCAR type files to CONTCAR type structures. \n"); */
/*        exit(0); */
/*      } */

    if (is.cart)  printf("This cif file has cartessian co-ordinates defined\n");
    if (is.fract) printf("This cif file has fractional co-ordinates defined\n");
    if (is.cart && is.fract)
      {
         printf("ERROR: Both cartessian and fractional co-ordinates defined in cif file, cannot mix!\n");
         exit(0);
      }
    if (!is.cart && !is.fract)
      {
         printf("ERROR: Neither cartessian or fractional co-ordinates defined in cif file, must have one or other!\n");
         exit(0);
      }

    scale_factor=1.0;
  }
else if (is.gulp)
  {
   printf("Master is a gulp file....\n");
//
// Note need just_count version of read_gulp 
//
   printf("ERROR: mallocing version cannot deal with gulp files...\n");
   exit(0);
//   open_file( &fp_gulp_input, master_input, "r");

//   read_gulp( fp_gulp_input, &header_line[0], &title_line[0],
//              &molecule[0], &shells[0], &date_line[0], &pbc,
//              &num_atoms, &num_shells, &num_of_mols, 
//              &num_mol_members[0], &mol_number[0],
//              &abc[0], &been_before, &top_bit[0],
//              &num_top_chars, &bottom_bit[0], &num_bottom_chars,
//              &space_group, &spec_charges[0], &num_species);

//    if (read_restart)
//      {
//        printf("ERROR: Restart directive given but master file is a gulp format file\n");
//        printf("       Restart is only intended for copying fixed atom flags over \n");
//        printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
//        exit(0);
//      }

//   scale_factor=1.0;

//   printf("For GULP inputs need to set is.cart or is.fract\n");

//    fclose(fp_gulp_input);
    set_labels = FALSE;
  }
else if (is.vasp)
  {
    printf("Master is a POSCAR or CONTCAR file....\n");
    open_file( &fp_vasp_input, master_input, "r");

/*** Use been_before flag to pick up fixed atom flags from master ****/
    been_before=TRUE;
    just_count =TRUE;
    
    good_read = read_poscar( fp_vasp_input, &title_line[0],
                             p_molecule, &date_line[0], &pbc, &num_atoms,
                             &(have.labels), &num_of_mols, &num_mol_members[0], &mol_number[0],
                             &latt_vec[0], &recip_latt_vec[0], &abc[0], 
                             &been_before, &(is.fract), &(is.cart),
                             &types[0], &num_types, &scale_factor, p_fix_flags,
                             just_count);

    num_mol_members[0]=num_atoms;

    printf("read_poscar found %d atoms and %d types...\n", num_atoms, num_types);
    pbc=TRUE;
    fclose(fp_vasp_input);

/**** Reserve memory space for the molecule array ****/
/**** Note that num_atoms is set to actual number   **/
/**** of atoms.                                     **/

    printf("mallocing molecule for %d atoms and fix_flags\n", num_atoms);
    p_molecule=(atom*)malloc(num_atoms*sizeof(atom));
    p_fix_flags=(coord_flags*)malloc(num_atoms*sizeof(coord_flags));

/*** Now read in the atoms **/

    open_file( &fp_vasp_input, master_input, "r");

/*** Use been_before flag to pick up fixed atom flags from master ****/
    been_before=TRUE;
    just_count =FALSE;
    
    good_read = read_poscar( fp_vasp_input, &title_line[0],
                             p_molecule, &date_line[0], &pbc, &num_atoms,
                             &(have.labels), &num_of_mols, &num_mol_members[0], &mol_number[0],
                             &latt_vec[0], &recip_latt_vec[0], &abc[0], 
                             &been_before, &(is.fract), &(is.cart),
                             &types[0], &num_types, &scale_factor, p_fix_flags,
                             just_count);

    set_labels = FALSE;
/***** Get atom label information from the POTCAR file for the job *****/
/***** Unless these have already been gotten from the new POSCAR   *****/
/***** Added feature May 2016. Dave Willock                        *****/

    if (!have.labels)
      {
        if (have.potcar)
          {
             printf("No labels in POSCAR so checking with POTCAR file %s\n", potcar_input);
             open_file( &fp_vasp_input, potcar_input, "r");

             good_read = read_potcar( fp_vasp_input, p_molecule, 
                                      &types[0], &num_types); 

             fclose(fp_vasp_input);
          }
        else
          {
             printf("ERROR : No labels provided in POSCAR file and no POTCAR file defined\n");
             printf("ERROR : Please supply either a POSCAR file with labels line or the POTCAR file.\n");
             exit(0);
          }
      }
//    if (read_restart)
//      {
//        printf("Reading co-ordinates from the restart file : %s\n", restart_file);
//        printf("While taking flags for fixed atoms from the master file: %s\n", master_input);
//
//        open_file( &fp_vasp_input, restart_file, "r");

/*** Use been_before flag to ignore fixed atom flags from restart file ****/
//        been_before=FALSE;
    
//        good_read = read_poscar( fp_vasp_input, &title_line[0],
//                                 &molecule[0], &date_line[0], &pbc, &num_atoms,
//                                 &num_of_mols, &num_mol_members[0], &mol_number[0],
//                                 &latt_vec[0], &recip_latt_vec[0], &abc[0], 
//                                 &been_before, &(is.fract), &(is.cart),
//                                 &ion_number[0], &num_types, &scale_factor, &fix_flags[0]);

//        if (num_atoms != num_mol_members[0])
//          {
//             printf("ERROR: master file and restart file contain different numbers of atoms\n");
//             exit(0);
//          }
//      }
  }
/****************************************************************************/
/*** Process the OUTCAR file to make convergence movie DJW & EJ May 04 ******/
/****************************************************************************/
    else if (have.out)
      {
    printf("Master is an OUTCAR file....\n");
   printf("DEBUG: mallocing version trying to deal with OUTCAR file...\n");

      open_file( &fp_outcar, master_input, "r");

      if (!molecule_malloced) 
        {
          just_count = TRUE;
        }
      else
        {
          just_count = FALSE;
        }
      need.fermi=FALSE;
      if (need.dos) need.fermi=TRUE;
      
      if (need.energy)
        {
           num_energies=-1;
        }

      good_read = read_outcar( fp_outcar, p_molecule, &latt_vec[0],
                               &recip_latt_vec[0], &abc[0], &eigenvalues[0],
                               &eigenvecs[0], &forces, &chain_forces,
                               &num_atoms, &types[0], &num_types, &num_modes, need.freq,
                               need.force, need.energy, &(have.band), p_energy_vasp, need.fermi,
                               &fermi, &num_energies, just_count);

      pbc=TRUE; scale_factor=1.0;
      num_mol_members[0]=num_atoms;

      printf("DEBUG>> Back in main from read_outcar closing file pointer\n");

      fclose(fp_outcar);

      printf("OUTCAR states there are %d atoms...will use to malloc\n", num_atoms);

/**** Reserve memory space for the molecule array ****/
/**** Note that num_atoms is set to the actual      **/
/**** number of atoms.                              **/

      if (just_count)
        {
          printf("mallocing molecule for %d atoms\n", num_atoms);
          p_molecule=(atom*)malloc(num_atoms*sizeof(atom));
          printf("mallocing step molecule for %d atoms\n", end_num_atoms);
          p_step_molecule=(atom*)malloc(num_atoms*sizeof(atom));
          printf("mallocing fix_flags for %d atoms\n", end_num_atoms);
          p_fix_flags=(coord_flags*)malloc(num_atoms*sizeof(coord_flags));
          molecule_malloced=TRUE;

          if (need.energy)
            {
              printf("\nmallocing for energies list with %d members from OUTCAR...\n", num_energies+1);
              p_energy_vasp=(double*)malloc((num_energies+1)*sizeof(double));
            }
          else
            {
              p_energy_vasp=(double*)malloc(sizeof(double));
            }

/**** Read in the actual atom data now that space has been reserved **/
/**** For master input check for fixing flags in group labels ****/

          open_file( &fp_outcar, master_input, "r");
          just_count = FALSE;
          been_before= FALSE;
          good_read = read_outcar( fp_outcar, p_molecule, &latt_vec[0],
                                   &recip_latt_vec[0], &abc[0], &eigenvalues[0],
                                   &eigenvecs[0], &forces, &chain_forces,
                                   &num_atoms, &types[0], &num_types, &num_modes, need.freq,
                                   need.force, need.energy, &(have.band), p_energy_vasp, need.fermi,
                                   &fermi, &num_energies, just_count);

         printf("Back from read_outcar DAVE for the second time with %d atoms. num_mol_members[0]=%d.\n", num_atoms, num_mol_members[0]);

         fclose(fp_outcar);


         if (need.energy)
           {
             printf("Openning file %s to put the energies in...\n", energy_output);
             open_file( &fp_vasp_energy, energy_output, "w");

             sprintf(title_x, "Frame num");
             sprintf(title_y, "Pot. Energy");
             sprintf(title_z, " ");

             write_csv(fp_vasp_energy, title_x, title_y, title_z,
                       p_energy_vasp, p_fdum, p_fdum, FALSE,
                       FALSE, num_energies);

             fclose(fp_vasp_energy);
             exit(0);
           }
      }

     if (compare_modes) 
       {
          printf("Malloced version cannot cope with comparing modes at the moment as end_molecule has to be done...\n");
          exit(0);
       }
//     if (compare_modes) 
//       {
//         open_file( &fp_outcar, mode_compare_input, "r");
//
//         need_fermi=FALSE;
//         good_read = read_outcar( fp_outcar, &end_molecule[0], &new_latt_vec[0], 
//                                  &new_recip_latt_vec[0], &new_abc[0], 
//                                  &new_eigenvalues[0], &new_eigenvecs[0], &forces, 
//                                  &chain_forces, &num_new_atoms, &end_types[0], 
//                                  &num_new_types, &num_new_modes, need.freq,
//                                  need.force, need.energy, &(have.band), p_energy_vasp, need.fermi,
//                                  &fermi);
//
//         fclose(fp_outcar);

//         printf("Read reference modes from file : %s\n", mode_compare_input);

/*** simple checks ***/
//         if ( num_atoms != num_new_atoms)
//           {
//             printf("ERROR: Comparing modes but number of atoms in two structures do not match.\n");
//             exit(0);
//           }
//         else
//           {
//             printf("Comparing atoms by co-ordinates:\n");
//
//             for (iatom=0; iatom<=num_atoms; iatom++)
//               {
//                 dmin=-1.0;
//                 map_atoms[iatom]=-1;
//                 for (jatom=0; jatom<=num_atoms; jatom++)
//                   {
//                     dx = molecule[iatom].x - end_molecule[jatom].x;
//                     dy = molecule[iatom].y - end_molecule[jatom].y;
//                     dz = molecule[iatom].z - end_molecule[jatom].z;
//
//                     min_image( &dx, &dy, &dz,
//                                &recip_latt_vec[0], &latt_vec[0]);
//
//                     dist = sqrt(dx*dx + dy*dy + dz*dz);
//
//*** map_atoms matches atom indicies from master to reference i.e. map_atoms[iatom] is the ***/
/**** index for the atom in the reference list that matches iatom in the master list.       ***/
          
//                     if (dmin < 0.0 || dist < dmin)
//                       {
//                          dmin=dist;
//                          map_atoms[iatom] = jatom;
//                       }
//                   }
//                 if (map_atoms[iatom] >= 0) 
//                   {
//                      printf("Closest atom to master atom %d in reference is %d they are %10.6f apart\n",
//                                     iatom, map_atoms[iatom], dmin);
//                   }
//                 else
//                   {
//                      printf("ERROR: Problem mapping master atoms to reference.\n");
//                      exit(0);
//                   }
//               }
//           }
//         if ( num_modes != num_new_modes)
//           {
//             printf("ERROR: Comparing modes but number of modes in two structures do not match.\n");
//             exit(0);
//           }
//
//         printf("Can compare eigenvecs as number of atoms and modes match.\n");
//
//         printf("Comparing %d modes for %d atoms.\n", num_modes, num_atoms);

/*** Check reference modes are orthonormal ****/

//         printf("Checking that reference modes form an orthonormal set.\n");
//         for (imode=0; imode<num_modes; imode++)
//           {
//              tot_dot=0.0;
//              for ( jmode=0; jmode < num_modes; jmode++)
//                {
//                  dot=0.0;
//                  for (iatom= 0; iatom <= num_atoms; iatom++)
//                     { 
//                       dot += new_eigenvecs[imode].dx[iatom]*new_eigenvecs[jmode].dx[iatom];
//                       dot += new_eigenvecs[imode].dy[iatom]*new_eigenvecs[jmode].dy[iatom];
//                       dot += new_eigenvecs[imode].dz[iatom]*new_eigenvecs[jmode].dz[iatom];
//                     }
//                  if (imode == jmode )
//                    {
//                      printf("Mode %d with self                    : %10.6f\n", 
//                                                                       imode,dot);
//                    }  
//                  else
//                    {
//                      tot_dot+=fabs(dot);
//                    }
//                }
//              printf("Abs. total of %d mode with all others: %10.6f\n\n", imode, tot_dot);
//           }
                
/*** Open file for csv output ***/
//         open_file( &fp_outcar, "mode_match.csv", "w");
//         fprintf(fp_outcar,"%s, ,%s\n",master_input, mode_compare_input);
//         fprintf(fp_outcar,"Mode Num, WaveNo (cm-1), Mode Num, WaveNo (cm-1), diff, Composition\n"); 
/*** Work out vector distance between eigenvectors of two structures ***/
//         for (imode=0; imode<num_modes; imode++)
//           {
//              dmin=-1.0;
//              for ( jmode=0; jmode < num_modes; jmode++)
//                {
/*** Allow for eigenvectors out of phase ***/
//                   for (isign=-1; isign <=1; isign+=2)
//                     {
//                       dist = 0.0;
//                       for (iatom= 0; iatom <= num_atoms; iatom++)
//                         { 
//                           diffx=eigenvecs[imode].dx[iatom]+isign*new_eigenvecs[jmode].dx[map_atoms[iatom]];
//                           diffy=eigenvecs[imode].dy[iatom]+isign*new_eigenvecs[jmode].dy[map_atoms[iatom]];
//                           diffz=eigenvecs[imode].dz[iatom]+isign*new_eigenvecs[jmode].dz[map_atoms[iatom]];
//
//                           dist += diffx*diffx + diffy*diffy + diffz*diffz;
//                         }

/*         printf("Dist mode %d to %d is %10.6f\n",imode+1, jmode+1, sqrt(dist));   */
/*** test if this is the closest ***/
                        
//                       if ( dmin < 0.0 || dist < dmin )
//                         { 
//                           dmin = dist;
//                           match = jmode;
//                         }
//                    }
/*                    printf("\n"); */
//                  dot_list[jmode]=0.0;
//                  for (iatom= 0; iatom <= num_atoms; iatom++)
//                    {
//                      dot_list[jmode]+=eigenvecs[imode].dx[iatom]*new_eigenvecs[jmode].dx[map_atoms[iatom]];
//                      dot_list[jmode]+=eigenvecs[imode].dy[iatom]*new_eigenvecs[jmode].dy[map_atoms[iatom]];
//                      dot_list[jmode]+=eigenvecs[imode].dz[iatom]*new_eigenvecs[jmode].dz[map_atoms[iatom]];
//                    }
//                }
//             printf("Mode %d ( %7.1f cm-1) in %s matches %d ( %7.1f cm-1) in %s diff %5.3f : comp ",
//                                         imode+1, eigenvalues[imode], master_input,
//                                         match+1, new_eigenvalues[match], mode_compare_input,
//                                         sqrt(dmin)); 

//             fprintf(fp_outcar,"%d,%7.1f,%d,%7.1f,%5.3f", imode+1, eigenvalues[imode], 
//                                                          match+1, new_eigenvalues[match],
//                                                          sqrt(dmin)); 

//             for ( jmode=0; jmode < num_modes; jmode++)
//                {
//                  if (fabs(dot_list[jmode]) > 0.1) 
//                    {
//                      printf("%5.3f of %d ", dot_list[jmode], jmode+1);
//                      fprintf(fp_outcar,",%5.3f,%d", dot_list[jmode], jmode+1);
//                    }
//                }
//             printf("\n");
//             fprintf(fp_outcar,"\n");
//           }
//         fclose(fp_outcar);
//         exit(0);
//       }

//  if (read_restart)
//    {
//      printf("ERROR: Restart directive given but master file is a OUTCAR format file\n");
//      printf("       Restart is only intended for copying fixed atom flags over \n");
//      printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
//      exit(0);
//    }
//
//      num_of_mols = 1;
//      num_atoms++; 
//      num_mol_members[0]=num_atoms;
//      scale_factor=1.0;
//      is.cart = TRUE;
//      is.fract= FALSE;
//      pbc=TRUE;
  
        set_labels = FALSE;
      }

else if (is.siesta)
   {
    printf("Master is a SIESTA file....\n");
//
// Note need just_count version of read for siesta
//
   printf("ERROR: mallocing version cannot deal with siesta files...\n");
   exit(0);
//    printf("Will try to get structure from siesta file\n");
//    open_file( &fp_siesta_input, master_input, "r");

//  if (read_restart)
//    {
//      printf("ERROR: Restart directive given but master file is a siesta fdf format file\n");
//      printf("       Restart is only intended for copying fixed atom flags over \n");
//      printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
//      exit(0);
//    }

//    num_free=-1;
    
//    good_read = read_fdf( fp_siesta_input, &title_line[0],
//                          &molecule[0], &date_line[0], &pbc, &num_atoms,
//                          &num_of_mols, &num_mol_members[0], &mol_number[0],
//                          &latt_vec[0], &recip_latt_vec[0], &abc[0], 
//                          &been_before, &(is.fract), &(is.cart),
//                          &ion_number[0], &num_types, &scale_factor, &num_free,
//                          &fix_flags[0]);

//    scale_factor=1.0;

 /**** read_fdf returns the highest index of an atom in the list ***/
//    num_atoms++;
//    num_mol_members[0]=num_atoms;
//    pbc=TRUE;
//    fclose(fp_siesta_input);

//    cut_to_dot(master_input);
//    sprintf(master_input, "%s.vectors", master_input);

//    if (need.freq)
//      {
//        printf("Now getting eigenvalues and vectors from %s file\n",master_input);
//        open_file( &fp_siesta_input, master_input, "r");
//
//        read_siesta_vectors( fp_siesta_input, &eigenvalues[0], 
//                             &eigenvecs[0], &num_atoms, num_free,
//                             &num_modes);
//        fclose(fp_siesta_input);

//        cut_to_dot(master_input);
//        sprintf(master_input, "%s.fdf", master_input);
//      }
    set_labels = FALSE;
   }
else if (is.onetep)
   {
    printf("Master is a ONETEP file....\n");
//
// Note need just_count version of read_onetep
//
     printf("ONETEP version not malloced yet...\n");
     exit(0);

//      printf("Will try to get structure from onetep file\n");
//      open_file( &fp_onetep_input, master_input, "r");

//    if (read_restart)
//      {
//        printf("ERROR: Restart directive given but master file is a siesta fdf format file\n");
//        printf("       Restart is only intended for copying fixed atom flags over \n");
//        printf("       from POTCAR type files to CONTCAR tyoe structures. \n");
//        exit(0);
//      }

//      num_free=-1;
    
//      just_count=TRUE;
        num_lines = -1;
//      good_read = read_onetep( fp_onetep_input, &title_line[0],
//                            &molecule[0], &date_line[0], &pbc, &num_atoms,
//                            &num_of_mols, &num_mol_members[0], &mol_number[0],
//                            &latt_vec[0], &recip_latt_vec[0], &abc[0], 
//                            &been_before, &(is.fract), &(is.cart),
//                            &ion_number[0], &num_types, &scale_factor, &num_free,
//                            &fix_flags[0], file_line_ptrs, &line_lengths[0], 
//                            &num_lines, just_count );
/*** Malloc arrays ***/
//      fclose(fp_onetep_input); 
//      for (iline = 0; iline <= num_lines; iline++)
//        {
//         file_line_ptrs[iline] = (char*)malloc((line_lengths[iline]+10)*sizeof(char)); 
//         printf("mallocing line %d to %d\n",iline,line_lengths[iline]+10);
//        }

//      just_count=FALSE;
//      num_lines = -1;
//      open_file( &fp_onetep_input, master_input, "r");
//      good_read = read_onetep( fp_onetep_input, &title_line[0],
//                            &molecule[0], &date_line[0], &pbc, &num_atoms,
//                            &num_of_mols, &num_mol_members[0], &mol_number[0],
//                            &latt_vec[0], &recip_latt_vec[0], &abc[0], 
//                            &been_before, &(is.fract), &(is.cart),
//                            &ion_number[0], &num_types, &scale_factor, &num_free,
//                            &fix_flags[0], file_line_ptrs, &line_lengths[0], 
//                            &num_lines, just_count );

//      fclose(fp_onetep_input); 
//      for (iline = 0; iline <= num_lines; iline++)
//        {
//          printf("read onetep line in as >>%s%",file_line_ptrs[iline]);
//        }
//      scale_factor=1.0;

 /**** read_fdf returns the highest index of an atom in the list ***/
//      printf("%d atoms read in from onetep file\n",num_atoms);
//      num_mol_members[0]=num_atoms;
//      pbc=TRUE;
    set_labels = FALSE;
   }
else if (is.punch)
   {
    printf("Master is a ChemShell punch file....\n");
    printf("Cannot cope, but working in it......\n");

    pbc = FALSE;
    use_mols = FALSE;
    just_count = TRUE;
    open_file( &fp_punch_file, master_input, "r");

    num_atoms = read_punch(fp_punch_file, p_molecule, p_shells, &num_shells, p_fix_flags,
                           just_count);

    fclose(fp_punch_file);

/*** Size arrays ***/ 

    p_molecule=(atom*)malloc(num_atoms*sizeof(atom));
    p_fix_flags=(coord_flags*)malloc(num_atoms*sizeof(coord_flags));
    p_shells=(atom*)malloc(num_shells*sizeof(atom));

    just_count = FALSE;
    open_file( &fp_punch_file, master_input, "r");

    num_atoms = read_punch(fp_punch_file, p_molecule, p_shells, &num_shells, p_fix_flags,
                           just_count);
    find_molecules=FALSE;
    good_read = 0;
    num_of_mols = 0;
    num_mol_members[0]=num_atoms;

    fclose(fp_punch_file);

    printf("Read %d atoms and %d shells:\n", num_atoms, num_shells);

    p_atom=p_molecule;
    for (iatom=0; iatom<num_atoms; iatom++)
       {
         printf("%s (%s) %10.6f %10.6f %10.6f\n", p_atom->label,
                                                  p_atom->elem,
                                                  p_atom->x,
                                                  p_atom->y,
                                                  p_atom->z);

         p_shll = p_shells;
         for (ishell=0; ishell<num_shells; ishell++)
            {
              if ( iatom == p_shll->neighb[0] )
                {
                   printf("Shell: %s %10.6f %10.6f %10.6f\n", p_shll->label,
                                                              p_atom->x,
                                                              p_atom->y,
                                                              p_atom->z);
                   break;
                } 
              p_shll++;
            }

         p_atom++;
       }
//    exit(0);
    set_labels = FALSE;
   }


else
  {
    printf("ERROR: Do not understand format of master file\n");
    exit(0);
  }

/*************************************************************************/
/**** Process master file ************************************************/
/*************************************************************************/
/*****************************************************************************/
/*** Convert fractional co-ordinates to cartesian  ***************************/
/*****************************************************************************/

printf("Finished reading structural data\n\n");

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

p_atom=p_molecule;
for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {
      printf("%s (%s) %10.6f %10.6f %10.6f\n",p_atom->label,
                                              p_atom->elem,
                                              p_atom->x,
                                              p_atom->y,
                                              p_atom->z);
      p_atom++;
   }


/**** For car files generate cartesian lattice vectors from abc alpha beta gamma ****/
/**** do same for SIESTA fdf files, added Nov 06, Dave Willock.                   ****/

     if (pbc && (is.car || is.siesta || is.cif || is.pdb )) 
       {
/*** If required make vacuum gap by increasing c-vector ****/
         if (need.vgap)
           {
              printf("Adding a vacuum gap of %10.6f Angstroms\n", vgap);
              abc[2] += vgap;
           }
         cart_latt_vecs( &abc[0], &latt_vec[0], &recip_latt_vec[0]);
       }
     else if (need.vgap)
       {
         latt_vec[8] += vgap;
       }


/******* Sort out cartessian and fractional versions of co-ordinates ***/

printf("Found:\n %d atoms \n %d molecules \n",num_atoms,num_of_mols);

if (num_of_mols > MAXMOL)
  {
    printf("ERROR: Too many molecules....maximum parameterised as %d\n", MAXMOL);
    exit(0);
  }

printf("Mallocing space for fractional co-ordinate list for the %d atoms.\n", num_atoms);
p_fract_coords=(double*)malloc(3*num_atoms*sizeof(double));

if (is.fract)
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

/**** malloc space for the fractional co-ordinates ****/


     p_atom=p_molecule;
     p_this_frac=p_fract_coords;
     for (iloop=0; iloop < num_atoms; iloop++)
         {
            *p_this_frac= p_atom->x; p_this_frac++;
            *p_this_frac= p_atom->y; p_this_frac++;
            *p_this_frac= p_atom->z; p_this_frac++;

            printf("Sending molecule frac to fract_to_cart: %10.6f %10.6f %10.6f\n", 
                                  p_atom->x, p_atom->y, p_atom->z);

            fract_to_cart( &(p_atom->x), &(p_atom->y), &(p_atom->z),
                           *(p_this_frac-3),*(p_this_frac-2),*(p_this_frac-1),
                           &latt_vec[0] );

            printf("Now get cartesian: %10.6f %10.6f %10.6f\n\n", 
                                    p_atom->x, p_atom->y, p_atom->z);
            p_atom++;
         }
   }
else
   {

/*** Generate fractional version of co-ordinates too ***/
      printf("Generating fractional co-ordinates from cartessian:\n");
     
      printf("\nreciprocal lattice vectors:\n"); 
      printf("%10.6f %10.6f %10.6f\n", recip_latt_vec[0], recip_latt_vec[1], recip_latt_vec[2]);
      printf("%10.6f %10.6f %10.6f\n", recip_latt_vec[3], recip_latt_vec[4], recip_latt_vec[5]);
      printf("%10.6f %10.6f %10.6f\n\n", recip_latt_vec[6], recip_latt_vec[7], recip_latt_vec[8]);

      printf("Fractional version of atom list:\n");

      p_atom=p_molecule;
      p_this_frac=p_fract_coords;
      for (iatom= 0; iatom < num_atoms; iatom++)
         {
/**** Print out fractional co-ordinates version of structure ****/
           cart_to_fract( p_atom->x, p_atom->y, p_atom->z,
                          &fa, &fb, &fc,
                          &recip_latt_vec[0] );

           printf("%s %10.6f  %10.6f  %10.6f\n", p_atom->label, fa, fb, fc); 

           *p_this_frac= fa; p_this_frac++;
           *p_this_frac= fb; p_this_frac++;
           *p_this_frac= fc; p_this_frac++;

           p_atom++;
        }
   }
/*** If required make origin shift, move atoms by vector supplied ****/
if (need.oshift)
  {
     printf("Applying shift vector: %10.6f %10.6f %10.6f\n", oshift[0], oshift[1], oshift[2]);
     p_atom=p_molecule;
     for (iatom= 0; iatom < num_mol_members[0]; iatom++)
        {
           p_atom->x = p_atom->x + oshift[0];
           p_atom->y = p_atom->y + oshift[1];
           p_atom->z = p_atom->z + oshift[2];
           p_atom++;
        }
  }

printf("Structure as currently held: %s\n",master_input);
p_atom=p_molecule;
for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {
      printf("%s %10.6f %10.6f %10.6f\n",p_atom->label,
                                         p_atom->x,
                                         p_atom->y,
                                         p_atom->z);
      p_atom++;
   }

printf("So far have %d molecules\n",num_of_mols);
/**** Now sort out Neighbour information ****/
     start_mol = 0;
     for (iloop = 0; iloop <= num_of_mols; iloop++)
       {
          printf("Off to generate neighbours with start_mol= %d, num_mol_members= %d\n",
                                                 start_mol,num_mol_members[iloop]);

          generate_neighbours( p_molecule+start_mol, num_mol_members[iloop]-1,
                               &types[0], &num_types,
                               pbc, &recip_latt_vec[0],  &latt_vec[0],
                               &spec_charges[0], set_labels);

          printf("Molecule %d has %d members starting at %d\n",iloop,num_mol_members[iloop], start_mol);
          printf("Molecule %d has %d types of atoms.\n",iloop,num_types);

          for (iatom= start_mol; iatom < start_mol+num_mol_members[iloop]; iatom++)
             {
/**** If no group number make one up ****/
/**** Dave December 2005             ****/

       p_atom= p_molecule+iatom;

       if (strncmp(&(p_atom->group_no[0])," ",1) == 0) sprintf(&(p_atom->group_no[0]),"X"); 
       if (p_atom->group_no[0] == '\0') sprintf(&(p_atom->group_no[0]),"X"); 

/*** Assign Mass ***/
          p_atom->mass= atomic_mass_list(p_atom->elem);

                printf("%s (elem= %s, mass=%10.2f, g_no >>%s<<) with %d neighbours : ", 
                                p_atom->label, 
                                p_atom->elem, 
                                p_atom->mass,
                                p_atom->group_no, 
                                p_atom->num_neigh); 

                for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                  {
                     neigh_index= start_mol+p_atom->neighb[ineigh];
                     printf("%s ",(p_molecule+neigh_index)->label);
                  }
                printf("\n");
             }
/*** Carry out test for clash between atoms ***/
          clash=FALSE;
          for (iatom= start_mol; iatom < start_mol+num_mol_members[iloop]; iatom++)
             {
               p_atom= p_molecule+iatom;
               for (jatom= iatom+1; jatom < start_mol+num_mol_members[iloop]; jatom++)
                  {
                     p_atom2= p_molecule+jatom;
                     d = sqrt(atom_separation_squared(p_atom, p_atom2, 
                                                      TRUE, &recip_latt_vec[0], 
                                                      &latt_vec[0]));
                     if ( d < CLASH_TOL ) 
                       {
                         clash=TRUE;
                         printf("ERROR: atoms %s (%d) and %s (%d) have a separation of %10.6f.\n",
                                        p_atom->label, iatom, p_atom2->label, jatom, d);
                       }
                  }
             }
          if (clash)
             {
                printf("\n\nERROR: File contains atom-atom close contacts as listed above.\n");
                printf("ERROR: Structures generated from this file are likely to cause\n");
                printf("ERROR: VASP or other DFT codes to crash, so cannot continue...\n");
                exit(0);
             }

          start_mol = start_mol + num_mol_members[iloop];
       }

/*****************************************************/
/*** Put molecule numbers on each atom according to **/
/*** connectivity found.                            **/
/*****************************************************/

printf("Back from generate_neighbours with %d types\n", num_types);

if (find_molecules)
  {
     num_of_mols = find_mol( p_molecule, num_mol_members[0]-1 );
     printf("Back from find_mol\n");

     p_atom=p_molecule;
     for (iatom= 0; iatom < num_mol_members[0]; iatom++)
        {

          printf("atom %d %s (elem= %s, mass=%10.2f) is in molecule %d with %d neighbours : ", 
                            iatom,
                            p_atom->label, 
                            p_atom->elem, 
                            p_atom->mass,
                            p_atom->mol,
                            p_atom->num_neigh); 

          for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
             {
               neigh_index= p_atom->neighb[ineigh];
               printf("%s (%d)",(p_molecule+neigh_index)->label, neigh_index);
             }
          printf("\n");
          p_atom++;
        }
   }

/*****************************************************/
/*** Deal with the user requests *********************/
/*** At this point we have the co-ordinates read in  */
/*** from whatever source and neighbour and molecule */
/*** information assigned.                           */
/***                                                 */
/*** p_molecule : Atom co-ordinates and other info   */
/*** num_mol_members[0] : number of atoms in list    */
/*** latt_vec[] : a nine membered array containing   */
/***              the real space lattice vectors     */
/***   [0], [1], [2] : ax ay az                      */
/***   [3], [4], [5] : bx by bz                      */
/***   [6], [7], [8] : cx cy cz                      */
/*****************************************************/
/****** Report on DOS expectations *******/

if (need.dos)
  {
     if (is.vasp_dos)
      {
        printf("Will write density of states files based on DOSCAR file: %s\n", doscar_input);
      }
     else if (is.siesta_dos)
      {
        printf("Will write density of states files based on SIESTA .DOS file: %s\n", doscar_input);
        printf("ERROR: Malloced version cannot currently deal with density of states calculations\n");
        printf("ERROR: from SIESTA files...\n");
        exit(0);
      }
     else
      {
        printf("ERROR: need_dos requested but no DOS data supplied\n");
        exit(0);
      }
     printf("Output csv format files for DOS will be given stem     : %s_ndos\n", dos_output);
     printf("DOS smearing to use                                    : %10.6f\n", dos_smear);

     if (need.part_dos)
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

         if (!spd[0] && !spd[1] && !spd[2] && !spd[3] )  
          {
             printf("ERROR: Additional flag given for PDOS orbitals does not include s p d or f\n");
             exit(0);
          }
        printf("Will include ");
        if (spd[0]) printf("s ");
        if (spd[1]) printf("p ");
        if (spd[2]) printf("d ");
        if (spd[3]) printf("f ");
        printf("orbitals in PDOS contributions\n");
      }

     if (need.multi_dos)
       {
          printf("Multi-DOS to sum\n");
          printf("ERROR: Multi-DOS not available in malloced code...\n");
          exit(0);

          for (idos=0; idos <= num_dos_points; idos++)
            {
//               ndos[idos].up_dos = 0.0;
//               ndos[idos].down_dos = 0.0;
//               ndos[idos].up_totdos = 0.0;
//               ndos[idos].down_totdos = 0.0;
            }

          if (need.part_dos)
            {
              for (idos=0; idos <= num_dos_points; idos++)
                {
//                   pdos[idos].up_dos = 0.0;
//                   pdos[idos].down_dos = 0.0;
//                   pdos[idos].up_totdos = 0.0;
//                   pdos[idos].down_totdos = 0.0;
                }
            }

          dos_norm=0.0;
          for (iii=0; iii <= num_dos_files; iii++)
            {
//              open_file( &fp_doscar_input, dos_files[iii].name, "r");

//              read_doscar( fp_doscar_input, &temp_ndos[0], need_part_dos, 
//                           &part_dos_list[0], num_atoms_pdos, &temp_pdos[0],
//                           &spd[0], num_columns, num_dos_points );

               if (iii==0)
                 {
//                   for (idos=0; idos <= num_dos_points; idos++)
//                     {
//                        ndos[idos].energy=temp_ndos[idos].energy;
//                     }
                 }
               else
                 {
/***** check DOSCAR files were aligned on energy axis ****/
//                   for (idos=0; idos < num_dos_points; idos++)
//                     {
//                        if ( fabs(ndos[idos].energy-temp_ndos[idos].energy) > 0.0001)
//                          {
//                             printf("ERROR: Energy axes in DOSCAR files not aligned cannot combined them.\n");
//                             printf("       Use EMIN and EMAX in INCAR file to define the range you want\n");
//                             printf("       use the same values in all single k-point calculations.\n");
//                             exit(0);
//                          }
//                     }
                 }

//               for (idos=0; idos < num_dos_points; idos++)
//                 {
//                    ndos[idos].up_dos += dos_weights[iii]*temp_ndos[idos].up_dos;
//                    ndos[idos].down_dos += dos_weights[iii]*temp_ndos[idos].down_dos;
//                    ndos[idos].up_totdos += dos_weights[iii]*temp_ndos[idos].up_totdos;
//                    ndos[idos].down_totdos += dos_weights[iii]*temp_ndos[idos].down_totdos;
//
//                    dos_norm+= dos_weights[iii]; 
//                 }

/* Sum part_dos similarly */
  
//               if (need.part_dos)
//                 {
//                    for (idos=0; idos < num_dos_points; idos++)
//                      {
//                         pdos[idos].up_dos += dos_weights[iii]*temp_pdos[idos].up_dos;
//                         pdos[idos].down_dos += dos_weights[iii]*temp_pdos[idos].down_dos;
//                         pdos[idos].up_totdos += dos_weights[iii]*temp_pdos[idos].up_totdos;
//                         pdos[idos].down_totdos += dos_weights[iii]*temp_pdos[idos].down_totdos;
//                      }
//                 }

//               printf("Found %d dos points in file %d: %s\n", num_dos_points, iii,
//                                                              dos_files[iii].name );
//               fclose(fp_doscar_input);
           }
//         for (idos=0; idos < num_dos_points; idos++)
//           {
//              ndos[idos].up_dos = ndos[idos].up_dos / dos_norm;
//              ndos[idos].down_dos = ndos[idos].down_dos  / dos_norm;
//              ndos[idos].up_totdos = ndos[idos].up_totdos / dos_norm;
//              ndos[idos].down_totdos = ndos[idos].down_totdos / dos_norm;
//           }
  
//         if (need.part_dos)
//           {
//              for (idos=0; idos < num_dos_points; idos++)
//                {
//                   pdos[idos].up_dos = pdos[idos].up_dos / dos_norm;
//                   pdos[idos].down_dos = pdos[idos].down_dos  / dos_norm;
//                   pdos[idos].up_totdos = pdos[idos].up_totdos / dos_norm;
//                   pdos[idos].down_totdos = pdos[idos].down_totdos / dos_norm;
//                }
//           }
       }
     else
       {
          printf("Reading single DOSCAR file.\n");

          open_file( &fp_doscar_input, doscar_input, "r");
               
          just_count=TRUE;
          num_pdos_columns= count_doscar( fp_doscar_input, p_temp_ndos, need.part_dos, 
                                          &part_dos_list[0], num_atoms_pdos, p_temp_pdos,
                                          &spd[0], &num_dos_points, &num_ndos_columns, just_count );

          printf("main: number of columns for ndos  = %d\n", num_ndos_columns);
          printf("main: number of orbitals for pdos = %d\n", num_pdos_columns);
          printf("main: number of DOSCAR data points= %d\n", num_dos_points);
          fclose(fp_doscar_input);

// Mallocing for the p_ndos and p_pdos arrays

          p_ndos=(dos*)malloc(num_dos_points*sizeof(dos));
          p_pdos=(dos*)malloc(num_dos_points*sizeof(dos));

         
          open_file( &fp_doscar_input, doscar_input, "r"); 

          read_doscar( fp_doscar_input, p_ndos, need.part_dos, 
                       &part_dos_list[0], num_atoms_pdos, p_pdos,
                       &spd[0], num_ndos_columns, num_pdos_columns, 
                       num_dos_points );

          fclose(fp_doscar_input);
      }

     smear_dos( p_ndos, num_dos_points, dos_smear );

     sprintf(filename, "%s_ndos.csv",dos_output);
     open_file( &fp_dos_output, filename, "w");

     sprintf(title_x, "Energy (eV)");
     sprintf(title_y, "DOS");
     sprintf(title_z, "TDOS");
     have.tot = TRUE;
     if (num_ndos_columns==2) is.restricted=TRUE;
     else if (num_ndos_columns==4) is.restricted=FALSE;

     write_doscsv( fp_dos_output, title_x, title_y, 
                  title_z, p_ndos, have.tot, num_dos_points-1, is.restricted, fermi);

     fclose(fp_dos_output);

     if (need.part_dos)
       {
         smear_dos( p_pdos, num_dos_points, dos_smear );

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
         if (spd[3]) sprintf(title_y, "%sf", title_y); 
         
         sprintf(title_z, "TPDOS");
         have.tot = TRUE;

         write_doscsv( fp_dos_output, title_x, title_y, 
                      title_z, p_pdos, have.tot, num_dos_points-1, is.restricted, fermi);

       }
     exit(0);
  }

/****** Catch request for pdos with no dos file name defined ******/

if (need.part_dos && !need.dos)
  {
     printf("ERROR: Partial DOS requested without 'need dos_file' directive.\n");
     exit(0);
  }


if (need.hbond)
  {
/*****************************************************/
/*** Temporary code to identify oxygen species *******/
/*** to select add directive :                 *******/
/*** hbond                                     *******/
/*** to input file                             *******/
/*****************************************************/

/*** num_oxy_types will hold the number of expected types of oxygen ***/
/***                                                                ***/
/*** num_oxy_species will hold the number of each oxygen species found ***/
/*** These are structures of type atom_number so can also hold label for */
/*** oxygen type ***/
/*** num_oxy_species[0] is for O lattice ***/
/*** num_oxy_species[1] is for OH        ***/
/*** num_oxy_species[2] is for water O   ***/
/*** num_oxy_species[3] is for O alcohol ***/

     sprintf(&(num_oxy_species[0].atom_type),"lattice");
     sprintf(&(num_oxy_species[1].atom_type),"hydroxyl");
     sprintf(&(num_oxy_species[2].atom_type),"water");
     sprintf(&(num_oxy_species[3].atom_type),"alcohol");

/**** count the oxygens of each type ****/
     just_count=TRUE;
     count_oxygen_species(p_molecule, num_mol_members[0], oxy_list_ptrs, 
                          &num_oxy_species[0], num_oxy_types, just_count);

     printf("Back from count_oxygen_species have found:\n");

/**** Malloc arrays for oxygen indices ****/
     for (ispecies=0; ispecies < num_oxy_types; ispecies++)
       {
         printf("Oxygen of type %10s: %d\n", num_oxy_species[ispecies].atom_type, num_oxy_species[ispecies].num);
         oxy_list_ptrs[ispecies]=(int*)malloc((num_oxy_species[ispecies].num+1)*sizeof(int));
       }
        
/**** Now get the indices of the oxygens too **/
     just_count=FALSE;
     count_oxygen_species(p_molecule, num_mol_members[0], oxy_list_ptrs, 
                          &num_oxy_species[0], num_oxy_types, just_count);

     printf("Back from count_oxygen_species second time have found:\n");

     for (ispecies=0; ispecies < num_oxy_types; ispecies++)
       {
         printf("Oxygen of type %10s: %d\n", num_oxy_species[ispecies].atom_type, num_oxy_species[ispecies].num);

         if (num_oxy_species[ispecies].num > 0) 
                         block_print_int(oxy_list_ptrs[ispecies], num_oxy_species[ispecies].num, 10);
       }


     sprintf(filename, "count_Os.csv");

     open_file( &fp_monit, filename, "w");

     sprintf(title_x, "index");
     sprintf(title_y, "num");
     sprintf(title_z, " ");

     for (ispecies=0; ispecies < num_oxy_types; ispecies++)
                 data[ispecies] = (double) num_oxy_species[ispecies].num;

     write_csv(fp_monit, title_x, title_y, title_z,
               &data[0], p_fdum, p_fdum, FALSE,
               FALSE, num_oxy_types-1);

/*** Now try to find hydrogen bonds that involve these O atoms ****/
     find_oh_bonds( p_molecule, num_mol_members[0], oxy_list_ptrs, 
                    &num_oxy_species[0], num_oxy_types, 
                    TRUE,  &latt_vec[0], &recip_latt_vec[0]);

     printf("Bailing out for debug....\n\n\n");
     exit(0);
  }


/*****************************************************/
/*** End of temporary area ***************************/
/*****************************************************/


/*****************************************************/
/*** Permutation sub-routine *************************/
/*****************************************************/

if (need.permute)
  {
    printf("Off to permutation sub-routine\n");

    num_perms = permutation(p_molecule, &num_mol_members[0], &latt_vec[0], &recip_latt_vec[0],
                            pbc, &abc[0], &permute, &magmom[0], num_magmom);

/*******************************************************/
/* testing see if can write car file *******************/
/*******************************************************/

   if (need.car)
    {
      printf("Trying to write car file : %s\n", car_output);
      open_file( &fp_car_output, car_output, "w");

      start_frame = TRUE;
      write_car( fp_car_output, &header_line[0], &title_line[0], 
                 &c_title_line[0], &date_line[0], p_molecule, 
                 &mol_number[0], use_mols, pbc, &abc[0], 
                 num_mol_members[0], scale_factor, start_frame,
                 &super[0], &latt_vec[0], &recip_latt_vec[0],  
                 p_fix_flags, &magmom[0], num_magmom);

      fclose(fp_car_output);
    }
/*******************************************************/

    printf("Generated %d permutations.....\n", num_perms);
    exit(0);
  }


/*******************************************************/
/* Deal with cases for which Miller indices are used ***/
/*******************************************************/

if (have.miller)
 {
   if (num_miller == 1)
     {
        printf("Miller indices supplied as: %d %d %d\n", miller[0], miller[1], miller[2]);

        if (need.miller_sort)
          {

            printf("Call to sort_by_elem 1 before sort_by_miller : %d atoms, %d types\n",num_mol_members[0]-1, num_types);
            printf("First atom: %s %10.6f  %10.6f  %10.6f \n", p_molecule->label,
                                                               p_molecule->x,
                                                               p_molecule->y,
                                                               p_molecule->z);

            sort_by_elem( p_molecule, num_mol_members[0]-1, &types[0], num_types);

            printf("sorting atoms along Miller plane directions\n");
            printf("First atom: %s %10.6f  %10.6f  %10.6f \n", p_molecule->label,
                                                               p_molecule->x,
                                                               p_molecule->y,
                                                               p_molecule->z);

            sort_by_miller( p_molecule, num_mol_members[0]-1, &types[0], num_types, 
                            &latt_vec[0], &recip_latt_vec[0], &miller[0] );

//         }
//       else
//         {
            reorientate_cell(p_molecule, &num_mol_members[0], &latt_vec[0], &recip_latt_vec[0],
                             &abc[0], &miller[0], num_miller, &new_latt_vec[0], &new_recip_latt_vec[0],
                             &new_abc[0], &slab_mol[0], &num_slab_atoms);
          }
   
        look_for_layers = TRUE;
        if (look_for_layers)
          {
/*    NEW ROUTINE TO BE ADDED HERE */
            if (need.car)
              {
                 sprintf(step_output,"%s_%d%d%d.car",car_output,miller[0],miller[1],miller[2]);
                 printf("Trying to write car file : %s\n", step_output);
                 open_file( &fp_car_output, step_output, "w");

                 start_frame = TRUE;
                 write_car( fp_car_output, &header_line[0], &title_line[0], 
                            &c_title_line[0], &date_line[0], &slab_mol[0], 
                            &mol_number[0], pbc, use_mols, &new_abc[0], 
                            num_slab_atoms, scale_factor, start_frame,
                            &super[0], &new_latt_vec[0], &new_recip_latt_vec[0],  
                            p_fix_flags, &magmom[0], num_magmom);

                 fclose(fp_car_output);
              }
            if (need.pdb)
              {
                 sprintf(step_output,"%s_%d%d%d.pdb",pdb_output,miller[0],miller[1],miller[2]);
                 printf("Trying to write pdb file : %s\n", step_output);
                 open_file( &fp_pdb_output, step_output, "w");


                 write_pdb(fp_pdb_output, &slab_mol[0], 
                              &new_abc[0], num_slab_atoms, &super[0],
                              &new_recip_latt_vec[0], &new_latt_vec[0]);
                 fclose(fp_pdb_output);
             }
         }
      }
    else
      {
          printf("Have %d Miller index sets will try to construct shape...\n", num_miller);

        printf("Using first set of Miller indices supplied to define the surface of interest: %d %d %d\n", miller[0], miller[1], miller[2]);

// Need to make copies of Miller that are rotated to align with the slab axes
//
        reorientate_cell(p_molecule, &num_mol_members[0], &latt_vec[0], &recip_latt_vec[0],
                         &abc[0], &miller[0], num_miller, &new_latt_vec[0], &new_recip_latt_vec[0],
                         &new_abc[0], &slab_mol[0], &num_slab_atoms);

        hemi_rad = 1.83;
        origin[0]=slab_mol[0].x; origin[1]=slab_mol[0].y; origin[2]=slab_mol[0].z; 

// Shift the origin for a test 
//
        origin[0]=1.0528; origin[1]=1.0528;  origin[2]=1.0528; 

        printf("First defineing a hemi-spherical cluster of the radius %10.6f Angstroms.\n", hemi_rad);

        define_hemi(&slab_mol[0], &num_slab_atoms, &new_latt_vec[0], 
                    &new_recip_latt_vec[0], &new_abc[0], hemi_rad,
                    &origin[0], &cluster[0],&num_cluster_atoms);
        
            if (need.car)
              {
                 sprintf(step_output,"%s_hemi.car",car_output);
                 printf("Trying to write car file : %s for hemispherical cluster.\n", step_output);
                 open_file( &fp_car_output, step_output, "w");

                 printf("Created %d atoms in cluster:\n", num_cluster_atoms);
                 for (iatom=0; iatom<num_cluster_atoms; iatom++)
                    {
                       printf("%s %10.6f %10.6f %10.6f\n", cluster[iatom].label, cluster[iatom].x, 
                                                           cluster[iatom].y, cluster[iatom].z);
                    }

                 start_frame = TRUE;
                 write_car( fp_car_output, &header_line[0], &title_line[0], 
                            &c_title_line[0], &date_line[0], &cluster[0], 
                            &mol_number[0], use_mols, FALSE, &new_abc[0], 
                            num_cluster_atoms, scale_factor, start_frame,
                            &super[0], &new_latt_vec[0], &new_recip_latt_vec[0],  
                            p_fix_flags, &magmom[0], num_magmom);

                 fclose(fp_car_output);

                 sprintf(step_output,"%s_%d%d%d.car",car_output,miller[0],miller[1],miller[2]);
                 printf("Trying to write car file : %s\n", step_output);
                 open_file( &fp_car_output, step_output, "w");

                 start_frame = TRUE;
                 write_car( fp_car_output, &header_line[0], &title_line[0], 
                            &c_title_line[0], &date_line[0], &slab_mol[0], 
                            &mol_number[0], use_mols, pbc, &new_abc[0], 
                            num_slab_atoms, scale_factor, start_frame,
                            &super[0], &new_latt_vec[0], &new_recip_latt_vec[0],  
                            p_fix_flags, &magmom[0], num_magmom);

                 fclose(fp_car_output);

             }

        printf("Now DJW cutting the cluster based on the remaining miller indices list\n");
        printf("cutting from the hemi-spherical cluster...........................\n");

          for (imiller=0; imiller< num_miller; imiller++)
            {
               h_miller= 3*imiller; k_miller= h_miller+1; l_miller= k_miller+1;
               printf("Miller indices set %d supplied as: %d %d %d", 
                                imiller, miller[h_miller], miller[k_miller], miller[l_miller]);
               if (imiller==0)
                {
                   printf(" ...used for surface plane\n");
                }
               else
                {
                   printf("\n");
                }
            }

          printf("\n");
/*** Convert Miller indicies to new system orientation ***/
          for (imiller=1; imiller< num_miller; imiller++)
            {
               h_miller= 3*imiller; k_miller= h_miller+1; l_miller= k_miller+1;

               norm[0]= recip_latt_vec[0] * miller[h_miller] +  recip_latt_vec[1] * miller[k_miller] +  recip_latt_vec[2] * miller[l_miller];  
               norm[1]= recip_latt_vec[3] * miller[h_miller] +  recip_latt_vec[4] * miller[k_miller] +  recip_latt_vec[5] * miller[l_miller];  
               norm[2]= recip_latt_vec[6] * miller[h_miller] +  recip_latt_vec[7] * miller[k_miller] +  recip_latt_vec[8] * miller[l_miller];  

               new_miller[h_miller-3]= new_latt_vec[0] * norm[0] + new_latt_vec[1] * norm[1] + new_latt_vec[2] * norm[2];
               new_miller[k_miller-3]= new_latt_vec[3] * norm[0] + new_latt_vec[4] * norm[1] + new_latt_vec[5] * norm[2];
               new_miller[l_miller-3]= new_latt_vec[6] * norm[0] + new_latt_vec[7] * norm[1] + new_latt_vec[8] * norm[2];

               printf("Miller set %d give normal vector in old orientation of %10.6f %10.6f %10.6f\n",
                                                                                      imiller,norm[0], norm[1], norm[2]);

               printf("Miller indices set %d reorientated to : %3.1f %3.1f %3.1f\n", 
                                imiller, new_miller[h_miller-3], new_miller[k_miller-3], new_miller[l_miller-3]);

            }

        clus_rad=hemi_rad/2.0;
        define_cluster(&cluster[0], &num_cluster_atoms, &new_latt_vec[0], &new_recip_latt_vec[0],
                       &new_abc[0], &new_miller[0], num_miller-1, &slab_mol[0], &num_slab_atoms,
                       clus_rad);
         
        
            if (need.car)
              {
                 sprintf(step_output,"%s_hkl_clus.car",car_output);
                 printf("Trying to write car file : %s for hkl plane terminated cluster.\n", step_output);
                 open_file( &fp_car_output, step_output, "w");

                 start_frame = TRUE;
                 write_car( fp_car_output, &header_line[0], &title_line[0], 
                            &c_title_line[0], &date_line[0], &slab_mol[0], 
                            &mol_number[0], use_mols, FALSE, &new_abc[0], 
                            num_slab_atoms, scale_factor, start_frame,
                            &super[0], &new_latt_vec[0], &new_recip_latt_vec[0],  
                            p_fix_flags, &magmom[0], num_magmom);

                 fclose(fp_car_output);

             }
      }
  }

/*****************************************************/
/***** Need to extrapolate for onetep cell expansion**/
/*****************************************************/
if (need.expansion)
   {
    printf("Grabbing fractional coordinates\n");

    p_atom=p_molecule;
    for (iatom= 0; iatom < num_mol_members[0]; iatom++)
       {
         /**** Print out fractional co-ordinates version of structure ****/
         cart_to_fract( p_atom->x, p_atom->y, p_atom->z,
                        &fa, &fb, &fc,
                        &recip_latt_vec[0] );
                          
         p_atom->x = fa;
         p_atom->y = fb;
         p_atom->z = fc;

         p_atom++;
       }

    printf("Expanding unit cell by %3.1f\n",expansion[0]);
    printf("original lattice vectors a b c:\n");
    printf("     %10.6f %10.6f %10.6f \n",abc[0],abc[1],abc[2]);
    abc[0] = abc[0]*expansion[0];
    abc[1] = abc[1]*expansion[0];
    abc[2] = abc[2]*expansion[0];
    printf("new lattice vectors      a b c:\n");
    printf("     %10.6f %10.6f %10.6f \n",abc[0],abc[1],abc[2]);
    cart_latt_vecs( &abc[0], &latt_vec[0], &recip_latt_vec[0]);
    printf("Converting back to cartesian\n");

    p_atom=p_molecule;
    for (iatom= 0; iatom < num_mol_members[0]; iatom++)
       {
         /**** Print out cartesian co-ordinates version of structure *****/
         fract_to_cart( &fa, &fb, &fc,
                        p_atom->x,p_atom->y,p_atom->z,
                        &latt_vec[0] );
         p_atom->x = fa;
         p_atom->y = fb;
         p_atom->z = fc;
         printf("new expanded coordinates : %10.6f, %10.6f, %10.6f \n",
                        p_atom->x,p_atom->y,p_atom->z);

         p_atom++;
       }

    printf("done expansion\n");
   }


/*****************************************************/
/*** If requested write out start points *************/
/*****************************************************/
/*** If requested write out start points *************/

if (need.poscar)
  {
    open_file( &fp_vasp_output, poscar_output, "w");

/*** If needed do an element sort for poscar *********/
/*****************************************************/

    if (!need.miller_sort)
      {
         printf("Call to sort_by_elem 1 : %d atoms, %d types\n",num_mol_members[0]-1, num_types);
         sort_by_elem( p_molecule, num_mol_members[0]-1, &types[0], num_types);
      }
    else
      {
         printf("These co-ordinates already grouped by type for Miller sort....\n");
      }


printf("Structure as held prior to poscar writing: %s\n",master_input);
     start_mol = 0;
     for (iloop = 0; iloop <= num_of_mols; iloop++)
       {
          printf("Molecule %d has %d members starting at %d\n",iloop,num_mol_members[iloop], start_mol);
          printf("Molecule %d has %d types of atoms.\n",iloop,num_types);

          for (iatom= start_mol; iatom < start_mol+num_mol_members[iloop]; iatom++)
             {
/**** If no group number make one up ****/
/**** Dave December 2005             ****/


                 p_atom= p_molecule+iatom;

                 printf("atom %d: %s %10.6f %10.6f %10.6f (elem= %s, mass=%10.2f, g_no >>%s<<) with %d neighbours : ",
                                iatom,
                                p_atom->label,
                                p_atom->x,
                                p_atom->y,
                                p_atom->z,
                                p_atom->elem,
                                p_atom->mass,
                                p_atom->group_no,
                                p_atom->num_neigh);

                for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                  {
                     neigh_index= start_mol+p_atom->neighb[ineigh];
                     printf("%s ",(p_molecule+neigh_index)->label);
                  }
                printf("\n");
             }
          start_mol = start_mol + num_mol_members[iloop];
       }
    printf("Writing new POSCAR format file: %s\n", poscar_output);

    printf("First atom %s  %10.6f   %10.6f   %10.6f \n",
                  p_molecule->label, p_molecule->x,
                  p_molecule->y    , p_molecule->z);

    printf("num_types: %d num_atoms: %d num_mol_members[0]: %d\n", num_types, num_atoms, num_mol_members[0]);

    write_poscar( fp_vasp_output, p_molecule, p_fract_coords,
                  &types[0], num_types, &latt_vec[0], &scale_factor, num_mol_members[0],
                  &title_line[0], &c_title_line[0], pbc, need.poscar_frac, p_fix_flags,
                  need.zsort);

    fclose(fp_vasp_output);
  }


/*****************************************************/

if (need.onetep)
  {
    printf("writing ONETEP file >>%s<<\n", onetep_output);
    open_file ( &fp_onetep_output, onetep_output, "w");
    printf("found %d atom types\n",num_types); 
    
    write_onetep( fp_onetep_output, p_molecule, p_fract_coords,
                  &types[0], num_types, &latt_vec[0], &scale_factor, num_atoms,
                  &title_line[0], &c_title_line[0], pbc, is.fract, is.onetep,    
                  file_line_ptrs, num_lines);
  }

printf("Structure as held prior to poscar writing: %s\n",master_input);
     start_mol = 0;
     for (iloop = 0; iloop <= num_of_mols; iloop++)
       {
          printf("Molecule %d has %d members starting at %d\n",iloop,num_mol_members[iloop], start_mol);
          printf("Molecule %d has %d types of atoms.\n",iloop,num_types);

          for (iatom= start_mol; iatom < start_mol+num_mol_members[iloop]; iatom++)
             {
/**** If no group number make one up ****/
/**** Dave December 2005             ****/

                 p_atom= p_molecule+iatom;

                 printf("atom %d: %s %10.6f %10.6f %10.6f (elem= %s, mass=%10.6f, g_no >>%s<<) with %d neighbours : ", 
                                iatom,
                                p_atom->label, 
                                p_atom->x,
                                p_atom->y,
                                p_atom->z, 
                                p_atom->elem, 
                                p_atom->mass,
                                p_atom->group_no, 
                                p_atom->num_neigh); 

                for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                  {
                     neigh_index= start_mol+p_atom->neighb[ineigh];
                     printf("%s ",(p_molecule+neigh_index)->label);
                  }
                printf("\n");
             }
          start_mol = start_mol + num_mol_members[iloop];
       }

if (need.potcar)
  {
    printf("Asked to compile a POTCAR file.>>%s<<...\n", potcar_output);

    sprintf(command, "cat ");

    printf("Have types:\n");
    for (itype=0; itype <= num_types; itype++)
      {
         printf("%d %s\n", itype, types[itype].atom_type);

         sprintf(command, "%s POTCAR_%s ", command, types[itype].atom_type);
         
      }

    sprintf(command, "%s > %s", command, potcar_output);

    printf("So will system call with: %s\n", command);

    system(command);
  }

if (need.cif)
  {
    open_file( &fp_cif_output, cif_output, "w");

    printf("Writing new cif format file: %s\n", cif_output);

    write_cif( fp_cif_output, p_molecule, p_fract_coords,
               &types[0], num_types,
               &latt_vec[0], &abc[0], &scale_factor, num_mol_members[0], 
               &title_line[0], &c_title_line[0], pbc, TRUE,
               p_fix_flags, need.zsort);

    fclose(fp_cif_output);
  }

printf("Structure as currently held: %s\n",master_input);
p_atom= p_molecule;
for (iatom= 0; iatom < num_mol_members[0]; iatom++)
   {
      printf("%s %10.6f %10.6f %10.6f\n",p_atom->label,
                                         p_atom->x,
                                         p_atom->y,
                                         p_atom->z);
      p_atom++;
   }


if (need.car)
   {
      printf("Trying to write car file : %s\n", car_output);
      open_file( &fp_car_output, car_output, "w");

      start_frame = TRUE;
      write_car( fp_car_output, &header_line[0], &title_line[0], 
                 &c_title_line[0], &date_line[0], p_molecule, 
                 &mol_number[0], use_mols, pbc, &abc[0], 
                 num_mol_members[0], scale_factor, start_frame,
                 &super[0], &latt_vec[0], &recip_latt_vec[0],  p_fix_flags, 
                 &magmom[0], num_magmom);

      fclose(fp_car_output);

      if (need.monit)
        {
           count_monit=0;
           p_atom=p_molecule;
           for (iatom=0; iatom < num_atoms; iatom++) 
             {
               if (iatom == monit_at1) 
                 {
                   printf("Found atom %d it is %s %10.6f %10.6f %10.6f\n",
                           monit_at1+1, p_atom->label, p_atom->x, p_atom->y, p_atom->z);
                   p_at1=p_atom;
                   count_monit++;
                 }
               if (iatom == monit_at2) 
                 {
                   printf("Found atom %d it is %s %10.6f %10.6f %10.6f\n",
                           monit_at2+1, p_atom->label, p_atom->x, p_atom->y, p_atom->z);
                   p_at2=p_atom;
                   count_monit++;
                 }
               if (count_monit==2)
                 {
                   count_monit++;
                   inter_atom_vector(p_at1, p_at2, &vec[0]);

                   min_image( &vec[0], &vec[1], &vec[2], 
                              &recip_latt_vec[0], &latt_vec[0]);

                   dist = size_vector(&vec[0]);
        
                   printf("Monitor atoms distance: %10.6f\n",dist);
                 }

               p_atom++;
             }
       }
   }

/**********************************************************************/
/** Process XDATCAR file to create a trajectory from MD run ***********/
/**********************************************************************/

if (need.mdtraj)
  {
    open_file( &fp_input_frame, mdtraj_input, "r");
    
    printf("Processing trajectory from file %s\n", mdtraj_input);

    if (need.pdb)
     {
        open_file( &fp_pdb_output, pdb_output, "w");
     }
    if (need.arc)
     {
        open_file( &fp_arc_output, arc_output, "w");
     }


    num_frames = -1;

/*** On first call read_xdatcar will send back number of frames *****/
/*** from header line of XDATCAR file                           *****/

/*** No need just_count version as first call just reads atom counts **/
    read_xdatcar(fp_input_frame, &num_frames, p_step_molecule, num_atoms,
                 p_fix_flags, &latt_vec[0], &recip_latt_vec[0], &abc[0],
                 &is.fract, &is.cart, &atom_names[0],                                                    
                 &atom_num, &num_labels );

    printf("First call to read_xdatcar found %d atoms\n", num_atoms);

    p_step_molecule=(atom*)malloc(num_atoms*sizeof(atom));
    p_step_molecule_last=(atom*)malloc(num_atoms*sizeof(atom));

/** if montoring set up the monitor lists ***/
    if (need.monit)
      {
        if (monit_type1 > 1) p_mon_list1=(int*)malloc(num_atoms*sizeof(int));
        if (monit_type2 > 1) p_mon_list2=(int*)malloc(num_atoms*sizeof(int));
      }

/*** reset file pointer **/

    rewind(fp_input_frame);

/*** Note that num_frames is zero referenced ***/
/**** latt_vec is the array of lattice vectors in cartessian vector form ***/
/**** so the a-vector is latt_vec[0] to latt_vec[2] etc.                 ***/

    printf("Have %d frames to read\n", num_frames+1);
    printf("latt_vecs:\n");
    printf("%10.6f %10.6f %10.6f \n", latt_vec[0],  latt_vec[1],  latt_vec[2]);
    printf("%10.6f %10.6f %10.6f \n", latt_vec[3],  latt_vec[4],  latt_vec[5]);
    printf("%10.6f %10.6f %10.6f \n", latt_vec[6],  latt_vec[7],  latt_vec[8]);


/********************************************************************/
/** Read in frame by frame and process to required output format ****/
/********************************************************************/

    p_atom=p_molecule; p_step_atom=p_step_molecule;
    for (iatom=0; iatom < num_atoms; iatom++) 
      {
        *p_step_atom=*p_atom;
        p_atom++;p_step_atom++;
      }

    if (need.monit)
      {
/********************************************************************/
/*** set aside space for monitor arrays *****************************/
/*** monitor arrays are to hold lists of bond distances etc *********/
/********************************************************************/
     
         p_d_monit=(double*)malloc((num_frames+1)*sizeof(double));
         p_this_d_monit=p_d_monit;  


         if ( monit_type1 == 1 && monit_type2 == 1 )
           {
              printf("\nMonitoring atoms %d (%s) and %d (%s)\n\n", monit_at1+1, (p_step_molecule+monit_at1)->label,  
                                                              monit_at2+1, (p_step_molecule+monit_at2)->label); 
           }
         else if ( monit_type1 > 1 && monit_type2 == 1 )
           {
              printf("\nMonitoring with type 1 a list (%s) and type 2 an atom %d (%s)\n\n", mon_str1,  
                                                              monit_at2+1, (p_step_molecule+monit_at2)->label); 
           }
         else if ( monit_type1 == 1 && monit_type2 > 1 )
           {
              printf("\nMonitoring with type 1 an atom %d (%s) and type 2 a list (%s) \n\n", monit_at1+1, 
                                                                           (p_step_molecule+monit_at1)->label,  
                                                                           mon_str2);
           }
         else if ( monit_type1 > 1 && monit_type2 > 1 )
           {
              printf("\nMonitoring with type 1 a list (%s) and type 2 a list (%s)\n\n", mon_str1, mon_str2);  
           }
         else 
           {
              printf("ERROR: Do not understand a monitor setting with type1: %d and type2: %d\n", monit_type1, monit_type2);
              exit(0);
           }

/*** Make required lists */

         if ( monit_type1 > 1 )
           {
             p_atom=p_molecule;
             p_this_mon_list1 = p_mon_list1;
             num_mon_list1=-1;
             for (iatom=0; iatom < num_atoms; iatom++) 
               {
                 if (monit_type1 == 2)
                   {
                      if (strcmp( p_atom->elem, mon_str1) == 0 )
                       {
                         *p_this_mon_list1 = iatom;
                         num_mon_list1++;
                         p_this_mon_list1++;
                       }
                   }
                 else if (monit_type1 == 3)
                   {
                      if (strcmp( p_atom->label, mon_str1) == 0 )
                       {
                         *p_this_mon_list1 = iatom;
                         num_mon_list1++;
                         p_this_mon_list1++;
                       }
                   }
                 p_atom++;
               }
/** Check we have found some ***/
              if ( num_mon_list1 < 0 )
                {
                   printf("ERROR: No element/label found for given monitor label >>%s<< none in file!\n", mon_str1);
                   exit(0);
                }
              else
                {
                   printf("Found %d atoms of type %s to monitor\n", num_mon_list1+1, mon_str1);
                }
           }


         if ( monit_type2 > 1 )
           {
             p_atom=p_molecule;
             p_this_mon_list2 = p_mon_list2;
             num_mon_list2=-1;
             for (iatom=0; iatom < num_atoms; iatom++) 
               {
                 if (monit_type2 == 2)
                   {
                      if (strcmp( p_atom->elem, mon_str2) == 0 )
                       {
                         *p_this_mon_list2 = iatom;
                         num_mon_list2++;
                         p_this_mon_list2++;
                       }
                   }
                 else if (monit_type2 == 3)
                   {
                      if (strcmp( p_atom->label, mon_str2) == 0 )
                       {
                         *p_this_mon_list2 = iatom;
                         num_mon_list2++;
                         p_this_mon_list2++;
                       }
                   }
                 p_atom++;
               }
/** Check we have found some ***/
              if ( num_mon_list2 < 0 )
                {
                   printf("ERROR: No element/label found for given monitor label >>%s<< none in file!\n", mon_str2);
                   exit(0);
                }
              else
                {
                   printf("Found %d atoms of type %s to monitor\n", num_mon_list2+1, mon_str2);

                   printf("List is:\n");
                   p_this_mon_list2 = p_mon_list2;
                   for (iatom=0; iatom <= num_mon_list2; iatom++)
                      {
                         printf("List member %d is atom %d label %s elem %s\n", iatom, *p_this_mon_list2,
                                                                                (p_molecule+*p_this_mon_list2)->label,
                                                                                (p_molecule+*p_this_mon_list2)->elem);
                         p_this_mon_list2++;
                      } 
                }
           }
      }

/**** Loop over MD frames, i.e. each run through this loop will process one MD frame ****/
          p_atom=p_molecule; p_step_atom=p_step_molecule;
          for (iloop=0; iloop < num_atoms; iloop++)
             {
               *p_step_atom= *p_atom;
               p_atom++; p_step_atom++;
             }

    for (iframe=0; iframe <= num_frames; iframe++)
      {
        read_xdatcar(fp_input_frame, &num_frames, p_step_molecule, num_atoms,
                     p_fix_flags, &latt_vec[0], &recip_latt_vec[0], &abc[0],
                     &(is.fract), &(is.cart), &atom_names[0],                                                    
                     &atom_num, &num_labels );

        printf("Read in frame %d of %d\n", iframe+1, num_frames);

        if (is.fract)
          {
             printf("Processing the fractional co-ords\n");
             p_step_atom=p_step_molecule;

/**** Loop over atoms copying the x,y,z into temporary vector vec[] so that the ****/
/**** conversion can be carried out.                                            ****/
             for (iloop=0; iloop < num_atoms; iloop++)
                 {
                   vec[0] = p_step_atom->x;
                   vec[1] = p_step_atom->y;
                   vec[2] = p_step_atom->z;

                   fract_to_cart( &(p_step_atom->x), &(p_step_atom->y), &(p_step_atom->z),
                                  vec[0], vec[1], vec[2], &latt_vec[0] );
                   p_step_atom++;
                 }
          }

/**** Repair any translations that go across boundaries so that we can have a nice movie ***/
/**** Connie's work area Feb 2019  Begins                                                ***/

/*** Work out for each atom what its displacement has been since last frame ****/
     if (iframe > 0)
       {
          p_step_atom=p_step_molecule;
          p_step_atom_last=p_step_molecule_last;
          for (iloop=0; iloop < num_atoms; iloop++)
             {
                 inter_atom_vector(p_step_atom_last, p_step_atom, &vec[0]);

                 printf("Atom %d labeled %s moved by %10.6f  %10.6f  %10.6f \n", iloop, p_step_atom->label, vec );

                 dist = size_vector(&vec[0]);
			
                 ivec[0] = 0; ivec[1] = 0; ivec[2] = 0;
                 if ( dist > 0.1 ) 
                   {
                      printf(".................Oh! that's moved a long way....\n");

                      vec[0] = (latt_vec[0]*dx+  latt_vec[1]*dx+  latt_vec[2]*dx)/abc[0]; 
                      vec[1] = (latt_vec[3]*dy+  latt_vec[4]*dy+  latt_vec[5]*dy)/abc[1];
                      vec[2] = (latt_vec[6]*dz+  latt_vec[7]*dz+  latt_vec[8]*dz)/abc[2];

                      printf("vec of dots: %10.6f  %10.6f  %10.6f \n", vec[0], vec[1], vec[2] );

/**** Use the vector of dot products to decide on a integer array **/
                      ivec[0] = round(vec[0]/abc[0]);
                      ivec[1] = round(vec[1]/abc[1]);
                      ivec[2] = round(vec[2]/abc[2]);

                      printf("ivec of dots: %d      %d      %d\n", ivec[0], ivec[1], ivec[2] );

/*** NEXT : use the ivec array to knock off the right combination of whole lattice vectors from the co-ordinate of p_step_atom ***/
                     
                   }

/*** If displacement is bigger than tolerance shift atom by appropriate lattice vector ***/

/*ivec shows the amount of each lattice vector that is in the shift
*/

	
		/*if ( ivec > 0 );*/
		p_step_atom->x = p_step_atom->x - (latt_vec[0] * ivec[0]) - (latt_vec[3] * ivec[1]) - (latt_vec[6] * ivec[2]);
		p_step_atom->y = p_step_atom->y - (latt_vec[1] * ivec[0]) - (latt_vec[4] * ivec[1]) - (latt_vec[7] * ivec[2]);
		p_step_atom->z = p_step_atom->z - (latt_vec[2] * ivec[0]) - (latt_vec[5] * ivec[1]) - (latt_vec[8] * ivec[2]);
		
                p_step_atom++; p_step_atom_last++;
             }
       }

/*** Update step_molecule_last with the current molecule for next time around *****/

          p_step_atom=p_step_molecule;
          p_step_atom_last=p_step_molecule_last;
          for (iloop=0; iloop < num_atoms; iloop++)
             {
                 *p_step_atom_last = *p_step_atom;

                 p_step_atom++; p_step_atom_last++;
             }

/**** Connie's work area Feb 2019  Ends                                                  ***/ 

/*** Now the step molecule co-odrinates are cartessian even if the XDATCAR file was in fractional coordinates ***/

        if (need.pdb)
          {

            write_pdb(fp_pdb_output, p_step_molecule, 
                      &abc[0], num_atoms, &super[0],
                      &recip_latt_vec[0], &latt_vec[0]);
          }

/*** This is how to write an arc file *****/

        if (need.arc)
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
                            &date_line[0], p_step_molecule, 
                            &mol_number[0], use_mols, pbc, &abc[0], 
                            num_mol_members[0], scale_factor, start_frame,
                            &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                            &magmom[0], num_magmom);
              }
            else
              {
                 printf("NOT first frame\n");
                 start_frame=FALSE;
                 write_car( fp_arc_output, &header_line[0], 
                            &title_line[0], &c_title_line[0],
                            &date_line[0], p_step_molecule, 
                            &mol_number[0], use_mols, pbc, &abc[0], 
                            num_mol_members[0], scale_factor, start_frame,
                            &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags,
                            &magmom[0], num_magmom);
              }
          }

/**** Carry out monitoring tasks on this frame ****/

        if (need.monit)
          {

/*** Working point ***/
/*** monit_at1 and monit_at2 are indicies of specific atoms to monitor ****/
/*** Note: in read_input these are set so that monit_at1 < monit_at2   ****/
/*** Note: monit_at1 and monit_at2 refer to a list starting from 0     ****/

         if ( monit_type1 == 1 && monit_type2 == 1 )
           {
/*** Case 1 both items are atom numbers ***/

             p_step_atom=p_step_molecule+monit_at1; p_step_atom2=p_step_molecule+monit_at2; 
             *p_this_d_monit = sqrt(atom_separation_squared(p_step_atom, p_step_atom2, 
                                                            TRUE, &recip_latt_vec[0], 
                                                            &latt_vec[0]));

              printf("Found distance between atoms %d and %d for frame %d is %10.6f\n",
                                            monit_at1+1, monit_at2+1, iframe+1, *p_this_d_monit);

           }
/*** Case 2 item 1 is a list item 2 a single atom ***/
         else if ( monit_type1 > 1 && monit_type2 == 1 )
           {
              p_this_mon_list1= p_mon_list1;
              p_step_atom2=p_step_molecule+monit_at2; 
            
              for ( iatom = 0; iatom <= num_mon_list1; iatom++ )
                {
                  p_step_atom= p_step_molecule + *p_this_mon_list1;

/*** Check these are not bonded neighbours ***/
                  leave=FALSE;
                  for (ineigh=0; ineigh< p_step_atom->num_neigh; ineigh++)
                    {
                       neigh_index= p_step_atom->neighb[ineigh];
                       if ( neigh_index == monit_at2 ) leave=TRUE;
                    }

                  if ( !leave )
                    {
                      dist = sqrt(atom_separation_squared(p_step_atom, p_step_atom2, 
                                                          TRUE, &recip_latt_vec[0], 
                                                          &latt_vec[0]));

                      if (iatom == 0 || dist < *p_this_d_monit ) *p_this_d_monit = dist;
                    }

                  p_this_mon_list1++;
                }
           }
/*** Case 2 item 1 is a single atom item 2 a list ***/
         else if ( monit_type1 == 1 && monit_type2 > 1 )
           {
              p_this_mon_list2= p_mon_list2;
              p_step_atom=p_step_molecule+monit_at1; 
            
              for ( iatom = 0; iatom <= num_mon_list2; iatom++ )
                {
                  p_step_atom2= p_step_molecule + *p_this_mon_list2;
/*** Check these are not bonded neighbours ***/
                  leave=FALSE;
                  for (ineigh=0; ineigh< p_step_atom2->num_neigh; ineigh++)
                    {
                       neigh_index= p_step_atom2->neighb[ineigh];
                       if ( neigh_index == monit_at1 ) leave=TRUE;
                    }

                  if ( !leave )
                    {
                      dist = sqrt(atom_separation_squared(p_step_atom, p_step_atom2, 
                                                          TRUE, &recip_latt_vec[0], 
                                                          &latt_vec[0]));

                      if (iatom == 0 || dist < *p_this_d_monit )
                        {
                          *p_this_d_monit = dist;
                        }
                    }

                  p_this_mon_list2++;
                }
           }
/*** Case 4 both items are lists ***/
         else if ( monit_type1 > 1 && monit_type2 > 1 )
           {
              p_this_mon_list1= p_mon_list1;
              for ( iatom = 0; iatom <= num_mon_list1; iatom++ )
                {
                  p_step_atom= p_step_molecule + *p_this_mon_list1;

                  p_this_mon_list2= p_mon_list2;
                  for ( iatom2 = 0; iatom2 <= num_mon_list2; iatom2++ )
                    {
                       p_step_atom2= p_step_molecule + *p_this_mon_list2;
                       dist = sqrt(atom_separation_squared(p_step_atom, p_step_atom2, 
                                                           TRUE, &recip_latt_vec[0], 
                                                           &latt_vec[0]));

                       if (iatom == 0 || dist < *p_this_d_monit )
                         {
                           *p_this_d_monit = dist;
                         }

                      p_this_mon_list2++;
                    }

                  p_this_mon_list1++;
                }
             }
            p_this_d_monit++;
          }

/*** End of iframe loop ****/
      }

/*** Write out any monitoring csv files ***/
    if (need.monit)
      {
         if ( monit_type1 == 1 && monit_type2 == 1 )
           {
/*** Case 1 both items are atom numbers ***/
             sprintf(filename, "monit_sep_atom_%d_atom_%d.csv", monit_at1+1, monit_at2+1);
           }
/*** Case 2 item 1 is a list item 2 a single atom ***/
         else if ( monit_type1 > 1 && monit_type2 == 1 )
           {
             sprintf(filename, "monit_sep_%s_atom_%d.csv", mon_str1, monit_at2+1);
               
           }
/*** Case 2 item 1 is a single atom item 2 a list ***/
         else if ( monit_type1 == 1 && monit_type2 > 1 )
           {
             sprintf(filename, "monit_sep_atom_%d_%s.csv", monit_at1+1, mon_str2);
           }
/*** Case 4 both items are lists ***/
         else if ( monit_type1 > 1 && monit_type2 > 1 )
           {
             sprintf(filename, "monit_sep_%s_%s.csv", mon_str1, mon_str2);
           }

	 p_step_atom=p_step_molecule+monit_at1; p_step_atom2=p_step_molecule+monit_at2; 
         open_file( &fp_monit, filename, "w");

         sprintf(title_x, "Frame num");
         sprintf(title_y, "separation");
         sprintf(title_z, " ");

         write_csv(fp_monit, title_x, title_y, title_z,
                   p_d_monit, p_fdum, p_fdum, FALSE,
                   FALSE, num_frames);
      }

    free(p_step_molecule);

    if (need.pdb)
      {
         fclose(fp_pdb_output);
      }
    if (need.arc)
      {
         fclose(fp_arc_output);
      }
    fclose(fp_input_frame);
    exit(0);
  }
/**********************************************************************/
/*** convert master structure to other file formats if needed *********/
/**********************************************************************/

if (need.pdb)
   {
      open_file( &fp_pdb_output, pdb_output, "w");

      write_pdb(fp_pdb_output, p_molecule, 
                &abc[0], num_atoms, &super[0],
                &recip_latt_vec[0], &latt_vec[0]);

      if (!need.interpolate) fclose(fp_pdb_output);                
   }

if (need.gulp)
   {
      open_file( &fp_gulp_output, gulp_output, "w");

      space_group=1;
      write_gulp(fp_gulp_output, &title_line[0], 
                 p_molecule, num_atoms,
                 p_shells, num_shells, GULP_CART, &abc[0],
                 space_group, &recip_latt_vec[0], &latt_vec[0],
                 &spec_charges[0], num_types, &super[0],
                 need.shells, &shell_species[0], num_shell_species);

      if (!need.interpolate) fclose(fp_gulp_output);                
   }
/*****************************************************/
/***** Read End point file for interpolation cases ***/
/*****************************************************/
if ( need.interpolate || need.react_coord ) 
  {
if (is.end_car)
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
/*** For mallocing we should be able to create a      **/
/*** pointer for the end_molecule, which is the array **/
/*** holding the end point of the NEB. This should be **/
/*** the same length as the master molecule but should**/
/*** do a call to read_car first to check that the    **/
/*** two structures agree on the number of atoms.     **/
/*******************************************************/

    just_count=TRUE;

    good_read=  read_car( fp_car_input, &end_header_line[0], &end_title_line[0],
                          p_end_molecule, &end_date_line[0], &pbc, &end_num_atoms,
                          &end_num_of_mols, &end_num_mol_members[0], &end_mol_number[0],
                          &abc[0], &been_before, &groups, have.grp, &num_groups,
                          FALSE, p_fix_flags, just_count);

    fclose(fp_car_input);

    printf("Back from read_car for the second time with %d atoms for the end point structure..\n", end_num_atoms);

/**** Reserve memory space for the end_molecule array ****/
/**** Note that end_num_atoms is set to the actual      **/
/**** number of atoms.                                  **/
/****                                                   **/
/**** Must also need space for the inter_vec array      **/

    printf("mallocing end molecule for %d atoms\n", end_num_atoms);
    p_end_molecule=(atom*)malloc((end_num_atoms+1)*sizeof(atom));
    printf("mallocing step molecule for %d atoms\n", end_num_atoms);
    p_step_molecule=(atom*)malloc((end_num_atoms+1)*sizeof(atom));

    printf("mallocing inter-structural vector with %d components.\n", 3*end_num_atoms);
    p_inter_vec=(double*)malloc((3*end_num_atoms)*sizeof(double));

/**** Read in the actual atom data now that space has been reserved **/
/**** For master input check for fixing flags in group labels ****/

    open_file( &fp_car_input, end_input, "r");
    just_count = FALSE;
    been_before= FALSE;
    good_read=  read_car( fp_car_input, &header_line[0], &title_line[0],
                          p_end_molecule, &date_line[0], &pbc, &end_num_atoms,
                          &num_of_mols, &num_mol_members[0], &mol_number[0],
                          &abc[0], &been_before, &groups, have.grp,
                          &num_groups, TRUE, p_fix_flags, just_count);

    fclose(fp_car_input);

    printf("Back from read_car for the second time with %d atoms for the end point structure..\n", end_num_atoms);

/****************************************************************************/
/*** Make end a minimum image with master added Dec 05, Dave Willock ********/
/*** Unless asked not to.                 added Nov 09, Dave Willock ********/
/****************************************************************************/

    if (pbc)
     {
/**** match atoms in list of element supplied by distance ****/
/**** re-order end image to get this element in same order ***/
         p_at1=(atom*)malloc(sizeof(atom));

         if (need.match)
           {
             p_atom=p_molecule;
             for (iatom=0; iatom < num_mol_members[0]; iatom++)
               {
                  for ( imat=0; imat<= num_elem_to_match; imat++)
                    {
                      if (strcmp(p_atom->elem, elems_to_match[imat].name) == 0)
                        {
                          d2 = 1000.0; p_end_atom=p_end_molecule;
                          for (jatom=0; jatom < num_mol_members[0]; jatom++)
                           {
                              if (strcmp(p_end_atom->elem, p_atom->elem) == 0)
                                {
                                   inter_atom_vector(p_end_atom, p_atom, &vec[0]);

                                   min_image( &vec[0], &vec[1], &vec[2],
                                              &recip_latt_vec[0], &latt_vec[0]);

                                   this_d2 = size_vector(&vec[0]);

/*** Remember index and set d2 closer ***/
                                   if (this_d2 < d2)
                                     {
                                       imatch = jatom;
                                       d2 = this_d2;
                                     }
                                }
                              p_end_atom++; 
                           }
/*** here imatch is index of best matching atom ***/
                         if (imatch != iatom)
                           {
                              printf("Found end atom out of order during matching, start atom %d closest to end atom %d\n", iatom, imatch);

                              *p_at1 = *(p_end_molecule+iatom);
                              printf("here 1\n");
                              *(p_end_molecule+iatom) = *(p_end_molecule+imatch);
                              printf("here 2\n");
                              *(p_end_molecule+imatch) = *p_at1;
                              printf("here 3\n");
                           }
                       }
                    }
                 p_atom++;
               }
           }
         free(p_at1);

         if (end_min_image)
           {
             p_atom=p_molecule; p_end_atom=p_end_molecule;
             for (iatom=0; iatom < num_mol_members[0]; iatom++)
              {
                 dx = p_end_atom->x - p_atom->x;
                 dy = p_end_atom->y - p_atom->y;
                 dz = p_end_atom->z - p_atom->z;

                 min_image( &dx, &dy, &dz,
                           &recip_latt_vec[0], &latt_vec[0]);

                 p_end_atom->x = p_atom->x + dx;
                 p_end_atom->y = p_atom->y + dy;
                 p_end_atom->z = p_atom->z + dz;
 
                 p_atom++; p_end_atom++;
               }
            }

         p_atom=p_molecule; p_end_atom=p_end_molecule;
         for (iatom=0; iatom < num_mol_members[0]; iatom++)
           {

           d = sqrt(atom_separation_squared(p_atom, 
                                            p_end_atom,
                                            FALSE, &recip_latt_vec[0], 
                                            &latt_vec[0]));

           printf("Atoms %d : %s (Start) and %s (End) are %10.6f apart\n",
                          iatom, p_atom->label, 
                          p_end_atom->label, d);

           p_atom++; p_end_atom++;
           }
     }

    is.end_cart = TRUE; is.end_fract= FALSE;
  }
else
  {
    printf("Can only cope with end point car files ");
    printf("for interpolation or reaction coords at the moment....\n");
    exit(0);
  }

/****************************************************************************/
/** Check consistency of end point and start point **************************/
/****************************************************************************/

   printf("Checking num_atoms match...\n");
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
     printf("Assembling end point neighbours. for %d molecules...\n", end_num_of_mols);
     for (iloop = 0; iloop <= end_num_of_mols; iloop++)
       {
          printf("First off to generate neighbours for end point, start_mol = %d, members %d\n",
                                      start_mol, end_num_mol_members[iloop]);

          p_end_atom=p_end_molecule+start_mol;
          generate_neighbours( p_end_atom,
                               end_num_mol_members[iloop]-1,
                               &end_types[0], &end_num_types,
                               pbc, &recip_latt_vec[0],  &latt_vec[0],
                               &spec_charges[0], set_labels);

          printf("Molecule %d in end point structure has %d members starts at %d\n",
                                         iloop,end_num_mol_members[iloop], start_mol);
          
          for (iatom= start_mol; iatom < start_mol+end_num_mol_members[iloop]; iatom++)
             {
                printf("iatom = %d, assigning mass...\n", iatom);
                p_end_atom->mass= atomic_mass_list(p_end_atom->elem);

                printf("iatom = %d, printing details...\n", iatom);
                printf("%s (elem= %s, mass %10.2f) x: %10.6f y: %10.6f z: %10.6f, with %d neighbours : ",
                                                   p_end_atom->label, 
                                                   p_end_atom->elem, 
                                                   p_end_atom->mass,
                                                   p_end_atom->x,
                                                   p_end_atom->y,
                                                   p_end_atom->z,
                                                   p_end_atom->num_neigh);

                printf("iatom = %d, getting neighbours details.. for %d neighs.\n", iatom
                                                                                  , p_end_atom->num_neigh);
                for (ineigh=0; ineigh< p_end_atom->num_neigh; ineigh++)
                  {
                     neigh_index= start_mol+p_end_atom->neighb[ineigh];
                     printf("%s ",(p_end_molecule+start_mol+neigh_index)->label);
                  }
                printf("\n");

                printf("iatom = %d, incrementing pointer.......\n\n", iatom);
                p_end_atom++;
             }
          start_mol = start_mol + end_num_mol_members[iloop];
       }

/***************************************************************************/
/*** Look for molecules in end structure ***********************************/
/***************************************************************************/

end_num_of_mols = find_mol( p_end_molecule, end_num_mol_members[0]-1 );

p_end_atom=p_end_molecule;
for (iatom= 0; iatom < end_num_mol_members[0]; iatom++)
   {
     printf("%s (elem= %s) is in molecule %d with %d neighbours : ", 
                       p_end_atom->label, 
                       p_end_atom->elem, 
                       p_end_atom->mol,
                       p_end_atom->num_neigh); 

     for (ineigh=0; ineigh< p_end_atom->num_neigh; ineigh++)
        {
          neigh_index= p_end_atom->neighb[ineigh];
          printf("%s ",(p_end_molecule+neigh_index)->label);
        }
     printf("\n");
     p_end_atom++;
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

/********************************************************************************/
/*** Carry out react_coords run, added December 05 Dave Willock *****************/
/********************************************************************************/

 if (need.react_coord)
   {
      react_coords(p_molecule, p_end_molecule, num_atoms, &image_files[0],
                   num_images, &recip_latt_vec[0], &latt_vec[0], pbc);
      exit(0);
   }

/********************************************************************************/
/*** Carry out interpolation : output car files and POSCAR files ****************/
/********************************************************************************/

/*****************************************/
/** Sort start and end point   ***********/
/*****************************************/

         printf("Call to sort_by_elem 2 for master structure : %d atoms\n",num_mol_members[0]-1);
         sort_by_elem( p_molecule, num_mol_members[0]-1, &types[0], num_types);
         printf("Call to sort_by_elem 3 for end point structure : %d atoms\n",num_mol_members[0]-1);
         sort_by_elem( p_end_molecule, num_mol_members[0]-1, &types[0], num_types);


/*****************************************************/
/*** locate chosen molecule if present **************/
/*****************************************************/

  if (have.mol)
    {
       check= find_atom_label(p_molecule, num_mol_members[0]-1, 
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

           flag_chosen_atoms( p_molecule, p_end_molecule,
                              num_mol_members[0], &chosen_indices[0],
                              &num_chosen_atoms, mol_ind, &mol_cnt_lab[0],
                              &start_cofm[0], &end_cofm[0], &inter_cofm[0]);

           printf("Back in main have:\n");
           printf("start_cofm      : %10.6f %10.6f %10.6f\n", start_cofm[0], start_cofm[1], start_cofm[2]);
           printf("end_cofm        : %10.6f %10.6f %10.6f\n", end_cofm[0], end_cofm[1], end_cofm[2]);
           printf("\ninter cofm vec: %10.6f %10.6f %10.6f\n", inter_cofm[0], inter_cofm[1], inter_cofm[2]);

/*****************************************************/
/*** For late transition states get bonding info *****/
/*****************************************************/

           any_rec =TRUE;

           if (need.late) 
             {
               printf("Analysing bonding before and after for late transition\n");
               printf("Bonds at start:\n");
               starting_late=TRUE;
               now_close= FALSE;
               for ( iloop=0; iloop<=num_chosen_atoms; iloop++ )
                 {
                    iatom=chosen_indices[iloop];
                    p_atom=p_molecule+iatom;
                    printf("%d %s : ", iatom, p_atom->label);

                    for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                      {
                         iatom2 = p_atom->neighb[ineigh];
                         printf("%d %s ",iatom2, (p_molecule+iatom2)->label);

                      }
                    printf("\n");
                 } 

               num_old_bonds=-1;
               num_new_bonds=-1;
               printf("\nBonds at end:\n");
               for ( iloop=0; iloop<=num_chosen_atoms; iloop++ )
                 {
                    iatom=chosen_indices[iloop];
                    p_atom=p_molecule+iatom;
                    printf("%d %s : ", iatom, (p_end_molecule+iatom)->label);

                    for (ineigh=0; ineigh< (p_end_molecule+iatom)->num_neigh; ineigh++)
                      {
                         iatom2 = (p_end_molecule+iatom)->neighb[ineigh];
                         printf("%d %s ",iatom2, (p_end_molecule+iatom2)->label);
                      }

                    if ( p_atom->num_neigh > (p_end_molecule+iatom)->num_neigh )
                     {
/*****************************************************/
/*** This atom has lost bonds  ***********************/
/*** Work out which ones       ***********************/
/*****************************************************/
                        printf("This atom has lost bonds\n");
                        atom_losing= iatom;
                        for (ineigh=0; ineigh< p_atom->num_neigh; ineigh++)
                          {
                            iatom2 = p_atom->neighb[ineigh];

                            still_there=FALSE;
                            for (ineigh2=0; ineigh2< (p_end_molecule+iatom)->num_neigh; ineigh2++)
                              {
                                still_there = iatom2 == (p_end_molecule+iatom)->neighb[ineigh2]
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
                    else if ( p_atom->num_neigh < (p_end_molecule+iatom)->num_neigh )
                     {
/*****************************************************/
/*** This atom has gained a bond *********************/
/*****************************************************/
                        printf("This atom has gained bonds\n");
                        for (ineigh=0; ineigh< (p_end_molecule+iatom)->num_neigh; ineigh++)
                          {
                            iatom2 = (p_end_molecule+iatom)->neighb[ineigh];

                            still_there=FALSE;
                            for (ineigh2=0; ineigh2< p_atom->num_neigh; ineigh2++)
                              {
                                still_there = iatom2 == p_atom->neighb[ineigh2]
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
                        for (ineigh=0; ineigh< (p_end_molecule+iatom)->num_neigh; ineigh++)
                          {
                            iatom2 = (p_end_molecule+iatom)->neighb[ineigh];
                            p_atom2= p_molecule+iatom2;

                            matched_neigh = FALSE;
                            for (ineigh2=0; ineigh2< p_atom->num_neigh; ineigh2++)
                              {
                                 printf("Testing %s as partner of %s against %s as partner of %s\n",
                                               p_atom->label, (p_molecule+p_atom->neighb[ineigh2])->label,
                                               (p_end_molecule+iatom)->label, (p_end_molecule+iatom2)->label);

                                 if (iatom2 == p_atom->neighb[ineigh2])
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
                                                  p_atom->neighb[ineigh];

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
                   p_atom= p_molecule+iatom;
                   p_atom2= p_molecule+iatom2;
                   printf("%d %s to %d %s\n", iatom, p_atom->label,
                                              iatom2, p_atom2->label);
                 }

               printf("\nGained %d new bonds\n", num_new_bonds+1);
               for (iloop=0; iloop <= num_new_bonds; iloop++)
                 {
                   iatom = new_bonds[iloop].atom1;
                   iatom2= new_bonds[iloop].atom2;
                   p_atom= p_molecule+iatom;
                   p_atom2= p_molecule+iatom2;
                   printf("%d %s to %d %s\n", iatom, p_atom->label,
                                              iatom2, p_atom2->label);

                   have.transfer=FALSE;
                   for (jloop=0; jloop <= num_old_bonds; jloop++)
                     {
                        if    (iatom == old_bonds[jloop].atom1
                            || iatom == old_bonds[jloop].atom2)
                          {
                            have.transfer=TRUE;
                            atom_transfered = iatom;
                            atom_receiving  = iatom2;
                            strcpy(elem_rec, p_atom2->elem);
                          }
                        if    (iatom2 == old_bonds[jloop].atom1
                            || iatom2 == old_bonds[jloop].atom2)
                          {
                            have.transfer=TRUE;
                            atom_transfered = iatom2;
                            atom_receiving  = iatom;
                            strcpy(elem_rec, p_atom->elem);
                          }
                     }
                 }

               if (have.transfer)
                 {
                   now_close=FALSE;
                   printf("This is a transfer reaction with atom %d %s transfered from %d %s to %d %s\n",
                                 atom_transfered, (p_molecule+atom_transfered)->label,
                                 atom_losing,     (p_molecule+atom_losing)->label,
                                 atom_receiving,  (p_molecule+atom_receiving)->label);

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
                       printf("transfer vector is the direction to send the transfering atom toward the receiver\n");

                       p_trans_atom = p_molecule+atom_transfered;
                       p_rec_atom   = p_molecule+atom_receiving;
                       inter_atom_vector(p_trans_atom, p_rec_atom, &transfer_vec[0]);

/**** Get final bond length ****/

                        p_end_atom = p_end_molecule+atom_transfered;
                        p_end_atom2= p_end_molecule+atom_receiving;
                        inter_atom_vector(p_end_atom, p_end_atom2, &vec[0]);
                     }
                   else
                     {
/**** transfered atom is not from molecule ******/
                        printf("transfer is to the molecule\n");
                        p_atom = p_molecule+atom_transfered;
                        p_end_atom = p_end_molecule+atom_transfered;
  
                        inter_atom_vector(p_atom, p_end_atom, &transfer_vec[0]);

/**** Get final bond length ****/
                        if ( new_bonds[0].atom1 == atom_transfered )
                          {
                            iatom = new_bonds[0].atom2;
                          }
                        else
                          {
                            iatom = new_bonds[0].atom1;
                          }

                        p_end_atom2= p_end_molecule+iatom;
                        inter_atom_vector(p_end_atom2, p_end_atom, &vec[0]);
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
                                new_bond_length, p_end_atom->label, p_end_atom2->label);

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
         p_atom=p_molecule;
         p_end_atom=p_end_molecule; 
         p_vec= p_inter_vec;
         for (iloop=0; iloop < num_atoms; iloop++)
           {
             inter_atom_vector(p_atom, p_end_atom, p_vec);

             printf("start atom %s end atom %s: %10.6f %10.6f %10.6f \n", 
                            p_atom->label, p_end_atom->label, *p_vec, *(p_vec+1), *(p_vec+2));

             p_atom++; p_end_atom++; p_vec+=3;
           }
         printf("\n");

/***********************************************/
/** Generate intermediates and output files ****/
/***********************************************/

         step= 1.0 / ( num_inter + 1 );
         printf("Interpolation step : %10.6f\n", step);

/****************************************************/
/** Set up groups ***********************************/
/****************************************************/

        if (have.grp) 
          {           
             printf("Setting up groups.....\n");
             set_up_groups( p_molecule, p_end_molecule, num_atoms,
                            &groups,         num_groups, 
                            &orig_bond1[0],  &orig_bond2[0], 
                            &delta_bond1[0], &delta_bond2[0], 
                            &have, &need);
          }

/******************************************************************/
/******************************************************************/
/*** Carry out the interpolation steps ****************************/
/******************************************************************/
/******************************************************************/

         iframe=0;
         for (iloop=1; iloop <= num_inter+1; iloop++)
            {
               iframe++;
               printf("DEBUG>> start iloop frame %d abc now: %10.6f  %10.6f  %10.6f \n", iloop, abc[0], abc[1], abc[2]);

               if ( have.grp && need.angle )
                 {

                   printf("**********************************************\n");
                   printf("Group interpolation using Angle option step, iloop = %d, iframe = %d\n", iloop, iframe);
                   printf("**********************************************\n\n");

                   p_atom=p_molecule; p_step_atom=p_step_molecule;

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
/*                             grp_mem_ind = group_indices[igrp];  */
                             if ( grp_mem_ind == iatom ) check = TRUE;
                          }

                       if ( !check || iatom == groups.axis1_ind[0] || iatom == groups.axis2_ind[0] )
                         {
                            *p_step_atom = *p_atom;

                            move_atom(p_step_atom, step * iloop, p_inter_vec+iatom*3);
                         }
                       p_atom++; p_step_atom++;
                     }
/****************************************************************/ 
/*** Now deal with the special geometry for the grouped atoms ***/ 
/*** Rotate each one by an increment of its total theta.      ***/ 
/****************************************************************/ 
                   for (igrp=0; igrp<=num_in_group; igrp++)
                     {
/*                       grp_mem_ind = group_indices[igrp];   */
                       p_atom2=p_molecule+grp_mem_ind; p_step_atom= p_step_molecule+grp_mem_ind;

                       printf("Using angle interpolation for atom %d\n",grp_mem_ind);
                       if (grp_mem_ind != groups.axis1_ind[0] && grp_mem_ind != groups.axis2_ind[0] )
                         {
                            *p_step_atom = *p_atom2;

                       theta= iloop * p_atom2->theta/ (num_inter+1);
  
                       for (iii=0; iii < num_atoms; iii++) flags[iii]= FALSE;
                       flags[grp_mem_ind]= TRUE;

                       rotate_with_flags(p_step_molecule, &axis[0], &origin[0],
                                         theta, &flags[0], num_atoms-1);

                         }
                     }
                 }
               else if (have.grp && !need.late)
                 {

                   printf("**********************************************\n");
                   printf("Group interpolation step, iloop = %d, iframe = %d\n", iloop, iframe);
                   printf("**********************************************\n\n");

                   p_atom=p_molecule;
                   p_step_atom=p_step_molecule;

                   for (iatom=0; iatom < num_atoms; iatom++)
                     {
/************************************************************/
/**** Check if the atom is a member of a group defined   ****/
/**** Interpolate all atoms that stay the same first,    ****/ 
/**** i.e. all non-group atoms and the central one.      ****/ 
/************************************************************/
                       check = FALSE;
                       if (num_groups == 0)
                          {
//                             printf("Have only one group.....group 1 has %d members\n", groups.num_grp1);
                             is.centre = iatom == groups.centre[0];
                             for (igrp=0; igrp<=groups.num_grp1; igrp++)
                                {
                                   grp_mem_ind = groups.group1[igrp];  
                                   if ( grp_mem_ind == iatom ) check = TRUE;
                                }
                          }
                       else             
                          {
//                             printf("Have more than one group.....group 1 has %d and group 2 %d members\n", groups.num_grp1, groups.num_grp2);
                             is.centre = iatom == groups.centre[0] || iatom == groups.centre[1];
                             for (igrp=0; igrp<=groups.num_grp2; igrp++)
                                {
                                   grp_mem_ind = groups.group2[igrp];  
                                   if ( grp_mem_ind == iatom ) check = TRUE;
                                }
                          }

                       if ( !check || is.centre)
                         {
                            *p_step_atom = *p_atom;

                            move_atom(p_step_atom, step * iloop, p_inter_vec+iatom*3);
                         }

                       p_step_atom++; p_atom++;
                     }
/****************************************************************/ 
/*** Now deal with the special geometry for the grouped atoms ***/ 
/****************************************************************/ 
                   num_in_group = groups.num_grp1; 
                   if (num_groups > 0) num_in_group = groups.num_grp1 +groups.num_grp2+1; 

//                   printf("Have %d atoms in all groups\n", num_in_group);
                         
                   for (igrp=0; igrp<=num_in_group; igrp++)
                     {
/** test which group this atom belongs to ***/
                       if (igrp <= groups.num_grp1)
                         {
                           printf("This atom is in group 1\n");
                           grp_mem_ind = groups.group1[igrp];  
                           is.centre = grp_mem_ind == groups.centre[0];
                           grp_cnt= groups.centre[0];
                           this_delta = delta_bond1[igrp];
                           this_orig  = orig_bond1[igrp];
                         } 
                       else
                         {
                           printf("This atom is in group 2\n");
                           jgrp=igrp-groups.num_grp1-1;  
                           grp_mem_ind = groups.group2[jgrp];
                           is.centre = grp_mem_ind == groups.centre[1];
                           grp_cnt= groups.centre[1];
                           this_delta = delta_bond2[jgrp];
                           this_orig  = orig_bond2[jgrp];
                         }

                       p_step_cnt_atom=p_step_molecule+grp_cnt;

                       printf("Using different interpolation for atom %d\n",grp_mem_ind);
                       if (!is.centre)
                         {
                            p_step_atom= p_step_molecule+grp_mem_ind;
                            *p_step_atom= *(p_molecule+grp_mem_ind);
/***********************************/
/*** First make the normal step ****/
/***********************************/
                            move_atom(p_step_atom, step * iloop, p_inter_vec+grp_mem_ind*3);
                                                
/******************************************/
/*** Now check the vector to the centre ***/
/******************************************/
                            new_bond = this_orig + step * iloop * this_delta;

                            inter_atom_vector(p_step_cnt_atom, p_step_atom, &vec[0]); 

                            printf("Have member %d and centre %d\n", grp_mem_ind, grp_cnt);
                            printf("Centre to member vector: %10.6f %10.6f %10.6f\n", vec[0], vec[1], vec[2]);
                            printf("orig: %10.6f delta: %10.6f new: %10.6f\n", this_orig, this_delta, new_bond);

                            unit_vector(&vec[0]); 

                            p_step_atom->x = p_step_cnt_atom->x + new_bond * vec[0]; 
                            p_step_atom->y = p_step_cnt_atom->y + new_bond * vec[1]; 
                            p_step_atom->z = p_step_cnt_atom->z + new_bond * vec[2]; 

                            printf("New coords: %10.6f %10.6f %10.6f \n", p_step_atom->x,
                                                                          p_step_atom->y,
                                                                          p_step_atom->z );
                         }
                     }
                 }
               else if ( have.mol && !need.late )
                 {


                   printf("**********************************************\n");
                   printf("Processing molecule but NOT late interpolation\n");
                   printf("**********************************************\n\n");

                   cofm_step[0] = start_cofm[0] + step * iloop * inter_cofm[0];
                   cofm_step[1] = start_cofm[1] + step * iloop * inter_cofm[1];
                   cofm_step[2] = start_cofm[2] + step * iloop * inter_cofm[2];

                   interpolate_all_but_chosen(p_molecule, p_step_molecule, num_atoms,
                                              &chosen_indices[0], num_chosen_atoms, step*iloop, p_inter_vec );
/************************************************************/
/*** Now interpolate the rigid molecule atoms ***************/
/************************************************************/

                   if ( need.morph )
                     {
                       for (imol=0; imol<=num_chosen_atoms; imol++)
                         {
                           iatom = chosen_indices[imol];
                           p_atom=p_molecule+iatom; p_step_atom=p_step_molecule+iatom;
                           *p_step_atom = *p_atom;
                           p_step_atom->x = cofm_step[0] + p_atom->x
                                          + step * iloop * *(p_inter_vec+iatom*3);
                           p_step_atom->y = cofm_step[1] + p_atom->y
                                          + step * iloop * *(p_inter_vec+iatom*3+1);
                           p_step_atom->z = cofm_step[2] + p_atom->z
                                          + step * iloop * *(p_inter_vec+iatom*3+2);
                         }
                     }
                   else
                     {
                       for (imol=0; imol<=num_chosen_atoms; imol++)
                         {
                           iatom = chosen_indices[imol];
                           p_atom=p_molecule+iatom; p_step_atom=p_step_molecule+iatom;
                           *p_step_molecule = *p_atom;
                           p_step_atom->x = cofm_step[0] + p_atom->x;  
                           p_step_atom->y = cofm_step[1] + p_atom->y;  
                           p_step_atom->z = cofm_step[2] + p_atom->z;  
                         }
                     }
                 }
               else if ( have.mol && need.late )
                 {
/************************************************************/
/*** Use late transition state algorithms *******************/
/*** This is designed for a molecule delivering an      *****/
/*** atom to a surface, so the whole molecule initially *****/
/*** moves as a rigid body and the atomic interpolation *****/
/*** begins at the switch distance.                     *****/
/*** Last updated April 2019, Dave Willock ******************/
/************************************************************/

                   printf("*******************************************\n");
                   printf("Processing molecule with late interpolation\n");
                   printf("*******************************************\n\n");

                   if (have.transfer)
                     {

                       if (!now_close)
                         {
/**** make step molecule a copy of the start point ***/
                           p_step_atom=p_step_molecule; p_atom=p_molecule;
                           for (iatom= 0; iatom <= num_atoms; iatom++)
                             {
                               *p_step_atom= *p_atom;
                               p_step_atom++; p_atom++;
                             }
                           
/**** Case of rigid body interpolation ****/
                           cofm_step[0] = step * iloop * transfer_vec[0];  
                           cofm_step[1] = step * iloop * transfer_vec[1]; 
                           cofm_step[2] = step * iloop * transfer_vec[2];

                           printf("In rigid body section of motion\n");
                           printf("iloop = %d shift vector: %10.6f %10.6f %10.6f\n", iloop, cofm_step[0], cofm_step[1], cofm_step[2]);

/**** move chosen atoms according to centre of mass shift ****/
                           for (imol=0; imol<=num_chosen_atoms; imol++)
                             {
                               iatom = chosen_indices[imol]; 
                               p_atom=p_molecule+iatom; p_step_atom=p_step_molecule+iatom;

                               p_step_atom->x = p_atom->x + cofm_step[0];
                               p_step_atom->y = p_atom->y + cofm_step[1];
                               p_step_atom->z = p_atom->z + cofm_step[2];
                             }

/*** Deal with case of group on surface forming a new molecule during reaction ***/
/*** During rigid body translation of molecule just freeze the group.          ***/
                           if (have.grp)
                             {
                               printf("Also have a group defined with %d members:\n", groups.num_grp1+1);
                               for (igrp=0; igrp <= groups.num_grp1; igrp++)
                                 {
                                    grp_mem_ind= groups.group1[igrp];
                                    p_atom= p_step_molecule+grp_mem_ind;

                                    printf("%d: %d %s ", igrp, grp_mem_ind, p_atom->label);
                                    
                                    for (iii=0; iii < 3; iii++) *(p_inter_vec+3*grp_mem_ind+iii) = 0.0;
                                 }
                               printf("\n\n");
                             }

/*** Shift atoms that are not in the chosen list by linear interpolation **/
                           interpolate_all_but_chosen(p_molecule, p_step_molecule, num_atoms,
                                                      &chosen_indices[0], num_chosen_atoms, step*iloop, p_inter_vec );

/***********************************************************/
/* test for switch criterion *******************************/
/* atom_transfered and atom_receiving should have been set */
/***********************************************************/
                          
                           p_step_atom=p_step_molecule+atom_transfered; p_step_atom2=p_step_molecule+atom_receiving;
                           inter_atom_vector(p_step_atom, p_step_atom2, &vec[0]); 

                           latest_bond = size_vector(&vec[0]);

                           printf("latest bond length %10.6f transfering atom is %s receiver %s\n", latest_bond, p_step_atom->label, p_step_atom2->label );

                           if ( latest_bond < traj_switch*new_bond_length )
                             {
                               now_close = TRUE;
                               printf("criteron for switch during late transition is met....setting up for group method from here.\n");

/**** Group 1 will be the group that has been defined plus the transferred atom ***/
 
                               (groups.num_grp1)++;
                               groups.group1[groups.num_grp1] = atom_transfered;
                               strcpy(groups.cnt_lab1, p_step_atom2->label);

/**** Group 2 will be the remainder of the molecule losing the atom with a centre on the atom that has lost a neighbour *****/ 
                               num_groups=1;
                               groups.group_type2 = CENTRE_TYPE;
                               groups.num_grp2= num_chosen_atoms-1;

                               igrp=0;
                               for (imol=0; imol<=num_chosen_atoms; imol++)
                                 {
                                   if ( chosen_indices[imol] != atom_transfered ) 
                                     {
                                       groups.group2[igrp]=chosen_indices[imol];
                                       igrp++;
                                     }
                                 }
                               strcpy(groups.cnt_lab2, (p_molecule+atom_losing)->label);

                               printf("Setting up groups.....\n");
                               set_up_groups( p_step_molecule, p_end_molecule, num_atoms,
                                              &groups,         num_groups, 
                                              &orig_bond1[0],  &orig_bond2[0], 
                                              &delta_bond1[0], &delta_bond2[0], 
                                              &have, &need);

/** Reset start and the inter_structure vector ****/
                               printf("Inter-structure vector is now:\n");
                               p_atom=p_molecule;
                               p_step_atom=p_step_molecule; 
                               p_end_atom=p_end_molecule; 
                               p_vec= p_inter_vec;
                               for (iatom=0; iatom < num_atoms; iatom++)
                                 {
                                   *p_atom= *p_step_atom;
                                 
                                   inter_atom_vector(p_atom, p_end_atom, p_vec);

                                   printf("start atom %s %6.4f %6.4f %6.4f end atom %s: %6.4f %6.4f %6.4f inter_vec: %6.4f %6.4f %6.4f\n", 
                                                  p_atom->label,      p_atom->x,      p_atom->y,      p_atom->z, 
                                                  p_end_atom->label,  p_end_atom->x,  p_end_atom->y,  p_end_atom->z, 
                                                  *p_vec, *(p_vec+1), *(p_vec+2));

                                   p_atom++; p_step_atom++; p_end_atom++; p_vec+=3;
                                 }
                               printf("\n");

/*** reset iterpolation and ***/
/*** decide on new step     ***/
                               iloop = 1;
                               num_inter = num_after_switch - 1;
                               step= 1.0 / ( num_inter + 1 );

                               printf("Interpolation step after switch: %10.6f\n", step);
                               
                               need.late=FALSE;
                             }
                         }
                     }


                   if ( any_rec )
                     {
/******************************************************/
/** Any surface atom of the right element type ********/
/** will do for the switch, find closest       ********/
/******************************************************/

                        latest_bond = -1.0;
                        p_step_atom=p_step_molecule;
                        p_step_atom2=p_step_molecule+atom_transfered;

                        printf("Looping from zero for %d atoms\n", num_atoms);

                        for (iatom=0; iatom < num_atoms; iatom++) 
                          {
                            if (strcmp(p_step_atom->elem, elem_rec) == 0)
                              {
                                 inter_atom_vector(p_step_atom, p_step_atom2, &vec[0]); 
 
                                 temp_bond = size_vector(&vec[0]);

                                 printf("temp_bond is %10.6f\n", temp_bond);

                                 if (latest_bond < 0 || temp_bond < latest_bond) 
                                   {
                                     latest_bond = temp_bond;
                                     iatom2 = iatom;
                                   }
                              }
                            p_step_atom++;
                          }
                        if (latest_bond < 0)
                          {
                             printf("ERROR: Cannot find any suitable receiving element (%s) types for late transition state.\n", elem_rec);

                             exit(0);
                          }
                     }
                   else
                     {

                        printf("Doing this bit....\n");
                        p_step_atom2 = p_step_molecule+(p_step_molecule+atom_transfered)->neighb[0];

                        inter_atom_vector(p_step_atom, p_step_atom2, &vec[0]); 
 
                        latest_bond = size_vector(&vec[0]);
                    }

                   p_step_atom  = p_step_molecule+iatom2;

                 }
/************************************************************/
/*** Otherwise this is a straight linear interpolation ******/
/************************************************************/
               else
                 {

                   printf("**********************************************\n");
                   printf("Linear interpolation applied iloop = %d, iframe = %d\n", iloop, iframe);
                   printf("**********************************************\n\n");
                   p_atom=p_molecule; p_step_atom=p_step_molecule;
                   for (iatom=0; iatom < num_atoms; iatom++)
                     {
                       *p_step_atom = *p_atom;
                       printf("Atom at start : %s %10.6f %10.6f %10.6f\n",
                              p_step_atom->label, p_step_atom->x,
                              p_step_atom->y,     p_step_atom->z);
                       printf("Step length: %10.6f along vector %10.6f %10.6f %10.6f\n", step*iloop,
                                                  *(p_inter_vec+iatom*3),  
                                                  *(p_inter_vec+iatom*3+1),  
                                                  *(p_inter_vec+iatom*3+2));
                                    
                       move_atom(p_step_atom, step * iloop, p_inter_vec+iatom*3);

                       printf("Atom after move: %s %10.6f %10.6f %10.6f\n\n",
                              p_step_atom->label, p_step_atom->x,
                              p_step_atom->y,     p_step_atom->z);
                                    
                       p_atom++;p_step_atom++;
                     }
                 }

               sprintf(step_output,"%s%d", poscar_output, iframe);

               open_file( &fp_vasp_output, step_output, "w");

               if ( iloop < num_inter+1)
                 {
/***************************************/
/*** Check sum of position vectors *****/
/***************************************/
                   if (need.interpolate)
                    {
                       if (need.shift)
                          {

                             printf("Applying shift to keep centre of mass of whole system fixed...\n");

                             for ( iatom=0; iatom < 3; iatom++) sum_pos[iatom] = 0.0;

                             p_step_atom=p_step_molecule;
                             for (iatom=0; iatom < num_atoms; iatom++)
                               {
                                  sum_pos[0] += p_step_atom->x; 
                                  sum_pos[1] += p_step_atom->y; 
                                  sum_pos[2] += p_step_atom->z; 

                                  p_step_atom++;
                               }
               
                             printf("Sum of components for step %d : %10.6f %10.6f %10.6f\n",
                                    iframe, sum_pos[0], sum_pos[1], sum_pos[2]);

                             if (iframe == 1)
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

                                  p_step_atom=p_step_molecule;
                                  for (iatom=0; iatom < num_atoms; iatom++)
                                    {
//                                     p_step_atom->x += shift_centre[0]; 
//                                     p_step_atom->y += shift_centre[1]; 
//                                     p_step_atom->z += shift_centre[2]; 
 
                                       sum_check[0] += p_step_atom->x; 
                                       sum_check[1] += p_step_atom->y; 
                                       sum_check[2] += p_step_atom->z; 

                                       p_step_atom++;
                                    }
                                   printf("Sum of components for step %d after shift : %10.6f %10.6f %10.6f\n",
                                                      iloop, sum_check[0], sum_check[1], sum_check[2]);
                               }
                           }
                        }

                   write_poscar( fp_vasp_output, p_step_molecule,  
                                 p_fract_coords, &types[0], num_types, 
                                 &latt_vec[0], &scale_factor, num_atoms,
                                 &title_line[0], &c_title_line[0], pbc, 
                                 need.poscar_frac, p_fix_flags, need.zsort);
                 }
               else
                 {
//                   for (imol=0; imol<=num_chosen_atoms; imol++)
////                     {
//                       p_end_atom = p_end_molecule+chosen_indices[imol];
//                       p_end_atom->x += end_cofm[0];
//                       p_end_atom->y += end_cofm[1];
//                       p_end_atom->z += end_cofm[2];
//                     }

                   if (need.interpolate)
                    {
/***************************************/
/*** Check sum of position vectors *****/
/***************************************/
                      for ( iatom=0; iatom < 3; iatom++) sum_pos[iatom] = 0.0;

                      p_step_atom=p_step_molecule;
                      for (iatom=0; iatom < num_atoms; iatom++)
                        {
                           sum_pos[0] += p_step_atom->x; 
                           sum_pos[1] += p_step_atom->y; 
                           sum_pos[2] += p_step_atom->z; 
                           p_step_atom++;
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

                           p_step_atom=p_step_molecule;
                           for (iatom=0; iatom < num_atoms; iatom++)
                             {
//                                p_step_atom->x += shift_centre[0]; 
//                                p_step_atom->y += shift_centre[1]; 
//                                p_step_atom->z += shift_centre[2]; 
 
                                sum_check[0] += p_step_atom->x; 
                                sum_check[1] += p_step_atom->y; 
                                sum_check[2] += p_step_atom->z; 
 
                                p_step_atom++;
                             }

                         printf("Sum of components for step %d after shift : %10.6f %10.6f %10.6f\n",
                                           iloop, sum_check[0], sum_check[1], sum_check[2]);

                        }
                     }

                   write_poscar( fp_vasp_output, p_end_molecule,  
                                 p_fract_coords, &types[0], num_types, 
                                 &latt_vec[0], &scale_factor, num_atoms,
                                 &title_line[0], &c_title_line[0], pbc, 
                                 need.poscar_frac, p_fix_flags, need.zsort);

               fclose(fp_vasp_output);
                }

/***** write arc file if requested *************************/

               if (need.arc)
                 {
                   if ( iloop == 1 && !done_start_frame)
                     {
                        sprintf(c_title_line,"Inter_vasp generated arcfile record of interpolation.");
                        start_frame=TRUE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], p_molecule, 
                                   &mol_number[0], use_mols, pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                                   &magmom[0], num_magmom);

                        start_frame=FALSE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], p_step_molecule, 
                                   &mol_number[0], use_mols, pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                                   &magmom[0], num_magmom);

                        done_start_frame = TRUE;
                     }
                   else if ( iloop < num_inter+1)
                     {
                        start_frame=FALSE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], p_step_molecule, 
                                   &mol_number[0], use_mols, pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                                   &magmom[0], num_magmom);
                     }
                   else
                     {
                        start_frame=FALSE;
                        write_car( fp_arc_output, &header_line[0], 
                                   &title_line[0], &c_title_line[0],
                                   &date_line[0], p_end_molecule, 
                                   &mol_number[0], use_mols, pbc, &abc[0], 
                                   num_mol_members[0], scale_factor, start_frame,
                                   &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                                   &magmom[0], num_magmom);
                     }
                      
                 }

/***** write pdb file if requested *************************/

               if (need.pdb)
                 {
                   if ( iloop < num_inter+1)
                     {
                       write_pdb(fp_pdb_output, p_step_molecule, 
                             &abc[0], num_atoms, &super[0], &recip_latt_vec[0], &latt_vec[0]);
                     }
                   else
                     {
                       write_pdb(fp_pdb_output, p_end_molecule, 
                             &abc[0], num_atoms, &super[0], &recip_latt_vec[0], &latt_vec[0]);
                     }
                      
                 }

            if (need.car)
              {
/***** Strip .car from file to allow insertion of number ***/

               iletter=0;
               strcpy(step_output, car_output);
               iletter=0;
               while (strncmp(&step_output[iletter],".",1) != 0 && iletter < strlen(step_output)) iletter++;
               step_output[iletter]= '\0';
               sprintf(step_output,"%s%d.car", step_output, iframe);

               open_file( &fp_car_output, step_output, "w");

               start_frame=TRUE;

               if ( iloop < num_inter+1)
                 {
                   write_car( fp_car_output, &header_line[0], &title_line[0], 
                              &c_title_line[0], &date_line[0], p_step_molecule, 
                              &mol_number[0], use_mols, pbc, &abc[0], 
                              num_mol_members[0], scale_factor, start_frame,
                              &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags,  
                              &magmom[0], num_magmom);
                 }
               else
                 {
                   write_car( fp_car_output, &header_line[0], &title_line[0], 
                              &c_title_line[0], &date_line[0], p_end_molecule, 
                              &mol_number[0], use_mols, pbc, &abc[0], 
                              num_mol_members[0], scale_factor, start_frame,
                              &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags,  
                              &magmom[0], num_magmom);
                 }

               fclose(fp_car_output);
             }
           }
      if (need.arc) fclose(fp_arc_output);
   }
/****************************************/
/** Write out vibrational modes in ******/
/** Requested format               ******/
/****************************************/
    if (need.arc && need.freq)
      {
          if (which_mode<0)
              printf("Will write %d modes\n", num_modes);
          else
              printf("Have %d modes but just writing mode %d\n", num_modes, which_mode);
//
// Loop over modes but only do arc file for those required.
//
          for ( imode = 0; imode < num_modes; imode++)
            {
              if (which_mode < 0 || which_mode == imode+1)
                {
              sprintf(frame_output,"%s%d.arc",arc_output,imode+1);

              sprintf(c_title_line,"Mode %d frequency %10.6f cm-1",
                                             imode+1, eigenvalues[imode]);

              open_file( &fp_arc_output, frame_output, "w");

/*** Write header information ***/

              if (need.pdb)
               {
                 sprintf(frame_output_pdb,"%s%d.pdb",arc_output,imode+1);
                 printf("Will write pdbfile of mode to %s\n",frame_output_pdb);
                 open_file( &fp_pdb_output, frame_output_pdb, "w");

                 fprintf( fp_pdb_output,"COMENT %s\n", c_title_line);
               }
              if (need.poscar)
               {
                 sprintf(frame_output,"%sMode_%d_step%d",
                                      poscar_output, imode+1, iloop);

                 open_file( &fp_vasp_output, frame_output, "w");
      			printf("Trying to write poscar file : %s\n", poscar_output);
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

                  p_atom=p_molecule; p_step_atom=p_step_molecule;
                  for (iatom=0; iatom < num_atoms; iatom++)
                    {
                      *p_step_atom = *p_atom;
                      p_step_atom->x = p_step_atom->x 
                                            + step * eigenvecs[imode].dx[iatom];
                      p_step_atom->y = p_step_atom->y 
                                            + step * eigenvecs[imode].dy[iatom];
                      p_step_atom->z = p_step_atom->z 
                                            + step * eigenvecs[imode].dz[iatom];
                      p_atom++; p_step_atom++;
                    }

                  if (need.poscar)
                    {
                      printf("Call to sort_by_elem 4 : %d atoms\n",num_mol_members[0]-1);
                      sort_by_elem( p_step_molecule, num_mol_members[0]-1, 
                                    &types[0], num_types);

                      write_poscar( fp_vasp_output, p_step_molecule,  
                                    p_fract_coords, &types[0], 
                                    num_types, &latt_vec[0], 
                                    &scale_factor, num_mol_members[0],
                                    &title_line[0], &c_title_line[0], 
                                    pbc, need.poscar_frac, p_fix_flags, 
                                    need.zsort);

                    }

		if (need.car)
 		  {
			sprintf(frame_output,"%sMode_%d_step%d.car",
                                      car_output, imode+1, iloop);
      			printf("Trying to write car file : %s\n", car_output);
      			open_file( &fp_car_output, frame_output, "w");

   			write_car( fp_car_output, &header_line[0], &title_line[0],
                 		&c_title_line[0], &date_line[0], p_step_molecule,
                 		&mol_number[0], use_mols, pbc, &abc[0],
                 		num_mol_members[0], scale_factor, TRUE,
                                &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                                &magmom[0], num_magmom);

                       fclose(fp_car_output);
                  }


/***** write pdb file if requested *************************/

                  if (need.pdb)
                    {
                      write_pdb(fp_pdb_output, p_step_molecule, 
                                &abc[0], num_atoms, &super[0], &recip_latt_vec[0], &latt_vec[0]);
                    }

/*** add to the arc file for the mode **********************/

                  write_car( fp_arc_output, &header_line[0], 
                             &title_line[0], &c_title_line[0],
                             &date_line[0], p_step_molecule, 
                             &mol_number[0], use_mols, pbc, &abc[0], 
                             num_mol_members[0], scale_factor, start_frame,
                             &super[0], &latt_vec[0], &recip_latt_vec[0], p_fix_flags, 
                             &magmom[0], num_magmom);

                  start_frame = FALSE;
                }
              printf("About to close arc file\n");
              fclose(fp_arc_output);
              if (need.car)     fclose(fp_car_output);
              if (need.pdb)     fclose(fp_pdb_output);
              if (need.poscar)  fclose(fp_vasp_output);
            }
         }
//
// Write summary of calculated frequencies
//
        if (need.arc)
          {
             sprintf(mode_freqs_output, "%s_freq_sum.csv", arc_output);
             open_file( &fp_mode_freqs, mode_freqs_output, "w");
          }
        else
          {  
             open_file( &fp_mode_freqs, "freq_sum.csv", "w");
          }

        sprintf(title_x, "Mode num");
        sprintf(title_y, "freq. cm-1");
        sprintf(title_z, " ");

        write_csv(fp_mode_freqs, title_x, title_y, title_z,
                  &eigenvalues[0], p_fdum, p_fdum, FALSE,
                  FALSE, num_modes-1);
	
        fclose(fp_mode_freqs);
       }

/*********************************************/
/*** Give force data *************************/
/*********************************************/
     if (need.force)
       {
          printf("Forces read for %d atoms from OUTCAR file (eV/Angstrom)\n",
                               num_atoms);
          rms_x=0.0;
          rms_y=0.0;
          rms_z=0.0;
          rms_t=0.0;
          p_atom=p_molecule;
          for (iatom=0; iatom < num_atoms; iatom++)
            {
               printf("%4s : %10.6f %10.6f %10.6f\n", p_atom->label,
                                                      forces.dx[iatom],
                                                      forces.dy[iatom],
                                                      forces.dz[iatom]);
              
               rms_x += forces.dx[iatom] * forces.dx[iatom];
               rms_y += forces.dy[iatom] * forces.dy[iatom];
               rms_z += forces.dz[iatom] * forces.dz[iatom];

               p_atom++;
            }
          rms_x = sqrt( rms_x / num_atoms);
          rms_y = sqrt( rms_y / num_atoms);
          rms_z = sqrt( rms_z / num_atoms);

          rms_t = sqrt( rms_x * rms_x + rms_y * rms_y + rms_z * rms_z );

          printf("\n RMS : %10.6f %10.6f %10.6f\n", rms_x, rms_y, rms_z);
          printf("\n RMS magnitude : %10.6f\n\n", rms_t);

          if (have.band)
             {
                printf("This OUTCAR is one of a set from an elastic band calculation\n");

                printf("Chain spring forces: \n");
                rms_x=0.0;
                rms_y=0.0;
                rms_z=0.0;
                rms_t=0.0;
                p_atom=p_molecule;
                for (iatom=0; iatom < num_atoms; iatom++)
                  {
                     printf("%4s : %10.6f %10.6f %10.6f\n", p_atom->label,
                                                            chain_forces.dx[iatom],
                                                            chain_forces.dy[iatom],
                                                            chain_forces.dz[iatom]);
              
                     rms_x += chain_forces.dx[iatom] * chain_forces.dx[iatom];
                     rms_y += chain_forces.dy[iatom] * chain_forces.dy[iatom];
                     rms_z += chain_forces.dz[iatom] * chain_forces.dz[iatom];

                     p_atom++;
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
                p_atom=p_molecule;
                for (iatom=0; iatom < num_atoms; iatom++)
                  {
                     forces.dx[iatom]=forces.dx[iatom]-chain_forces.dx[iatom];
                     forces.dx[iatom]=forces.dy[iatom]-chain_forces.dy[iatom];
                     forces.dx[iatom]=forces.dz[iatom]-chain_forces.dz[iatom];

                     printf("%4s : %10.6f %10.6f %10.6f\n", p_atom->label,
                                                            forces.dx[iatom],
                                                            forces.dy[iatom],
                                                            forces.dz[iatom]);
              
                     rms_x += forces.dx[iatom] * forces.dx[iatom];
                     rms_y += forces.dy[iatom] * forces.dy[iatom];
                     rms_z += forces.dz[iatom] * forces.dz[iatom];

                     p_atom++;
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

printf("All done, freeing pointers.....\n");

free (p_molecule);
free (p_fix_flags);
printf("Freeing fract_coords...\n");
free (p_fract_coords);

if (need.dos)
  {
    free(p_ndos); free(p_pdos);
  }

if (need.interpolate)
  {
     free(p_end_molecule); free(p_inter_vec);
  }
if (need.freq) free(p_step_molecule);

return 0;
}
