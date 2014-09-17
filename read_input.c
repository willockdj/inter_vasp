/**********************************************************************************/
/**** read_input does just that, it obtains the input file names ******************/
/**** and the control parameters. Car files are read by read_car ******************/
/**** Adapted May 03 Dave Willock *************************************************/
/**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "maxima.h" 
#include "structures.h"
#include "reader.h" 
#include "global_values.h" 

char * tok_get(FILE *input_fp, int skip_lines, int lower_case);

int find_kind(char *token, int level);

double get_double(FILE *input_fp, int skip_lines, int *p_error);

int get_integer( FILE *input_fp, int skip_lines, int *p_error );

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
                int *p_num_inter, char *p_grp_cnt_lab,  char *p_grp_cnt_lab2,
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
                int *p_compare_modes, char *p_mode_compare_input, int *p_end_min_image )
 {
   int iloop, icat, skip, noskip, lower_case, leave_case;
   int end_of_input, itsanum, error, ip, iii, found_doscar;

   int   token, second_token, third_token;
   int first, got_last, srang, erang;
   char  *tok, *last_tok, *p_this_char;
   char *dash, cnum1[5], cnum2[5];

   *p_assess = FALSE;
   *p_analyse = FALSE;
   *p_have_grp = FALSE;
   *p_have_mol = FALSE;
   *p_have_out = FALSE;
   *p_have_incar = FALSE;
   *p_read_restart = FALSE;
   *p_need_md_run = FALSE;
   *p_end_min_image = TRUE;
   end_of_input= FALSE;
   found_doscar= FALSE;
   skip= TRUE;
   noskip= FALSE;
   lower_case= TRUE;
   leave_case= FALSE;
   tok = NULL;

   while (!end_of_input)
    {
      tok = tok_get( input_fp, skip, lower_case ); 
      skip=FALSE;

      if ( tok != NULL )
        {
          if ( strcmp(tok , END_OF_INPUT) == 0) end_of_input= TRUE;
          last_tok= tok;
          printf("Read tok >>%s<<\n",tok);

          token = find_kind(tok, PRIME_DIRECTIVE);

          switch (token)
            {
               case TITLE : first = TRUE;
                            while ( tok != NULL) 
                               {
                                 tok=tok_get( input_fp, skip, leave_case );
                 
                                 if ( tok != NULL )
                                   {
                                     if (first)
                                       {
                                          first = FALSE;
                                          strcpy( p_title, tok);
                                       }
                                     else
                                       {
                                          strcat( p_title, tok);
                                       }

                                     strcat( p_title, " ");
                                     printf("Title now >>%s<<\n", p_title);
                                   }
                               } 
                             break;

              case POTCAR_FILE: tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_potcar_input, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for potcar reference file\n");
                                     exit(0);
                                   }

                                break;

              case INCAR_FILE: tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_incar_input, tok);
                                     *p_have_incar=TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for incar reference file\n");
                                     exit(0);
                                   }

                                break;

              case RESTART: tok=tok_get( input_fp, skip, leave_case );
                            *p_read_restart = TRUE;

                            if (tok != NULL )
                              {
                                strcpy( p_restart_file, tok);
                              }
                            else
                              {
                                printf ("ERROR : No file name given for restart (CONTCAR) reference file\n");
                                exit(0);
                              }

                            break;

              case DOSCAR_FILE: tok=tok_get( input_fp, skip, leave_case );

                                if (found_doscar)
                                  {
                                     if (*p_is_siestados)
                                      {
                                        printf("ERROR: Multiple DOS files not supported for SIESTA\n");
                                        exit(0);
                                       }
                                    ++(*p_num_dos_files);
                   /**** This is a secondary occurance so must have multiple DOSCAR files ****/

                                    *p_need_multi_dos= TRUE; 

                                    if (*p_num_dos_files == 1)
                                      {
                                        strcpy(p_dos_files->name, p_doscar_input);
                                        printf("Copied over >>%s<<\n",p_dos_files->name);
                                        p_dos_files++;
                                      } 
 
                                     if (tok != NULL )
                                        {
                                          strcpy( p_dos_files->name, tok);
                                          printf("Read in >>%s<<\n",p_dos_files->name);
                                          p_dos_files++;
   /*** Get weight *****/ 
                                          tok=tok_get( input_fp, noskip, leave_case );

                                          if (tok != NULL )
                                             {
                                                *p_dos_weights= atof(tok);
                                                p_dos_weights++;
                                             }
                                          else
                                             {
                                                printf ("ERROR : With multiple doscar files expect weights\n");
                                                printf ("        as an integer following the file name.\n");
                                                exit(0);
                                             }
                                        }
                                     else
                                        {
                                          printf ("ERROR : No file name given for doscar reference file\n");
                                          exit(0);
                                        }
                                  }

                                else
                                  {
                                     found_doscar = TRUE;
                                     if (tok != NULL )
                                        {
                                          strcpy( p_doscar_input, tok);

                                          if (strstr(p_doscar_input, ".DOS") != NULL)
                                             {
                                               printf("You input a SIESTA DOS file string\n");
                                               *p_is_siestados= TRUE;
                                               *p_is_vaspdos= FALSE;
                                              }
                                          else
                                              {
                                               *p_is_siestados= FALSE;
                                               *p_is_vaspdos= TRUE;

                                               tok=tok_get( input_fp, noskip, leave_case );

                                               if (tok != NULL )
                                                  {
                                                     *p_dos_weights= atof(tok);
                                                     p_dos_weights++;
                                                  }
                                             }
                                        }
                                     else
                                        {
                                          printf ("ERROR : No file name given for doscar reference file\n");
                                          exit(0);
                                        }
                                  }

                                break;

              case MDTRAJ :  tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     *p_need_MDtraj = TRUE;
                                     strcpy( p_mdtraj_input, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for MD trajectory file\n");
                                     exit(0);
                                   }

                                break;

              case DOS_SMEAR : *p_dos_smear = get_double( input_fp, skip, &error );
                               break;

              case GROUP_CENTRE: tok=tok_get( input_fp, skip, leave_case );
                                 *p_have_grp = TRUE;
                                if (tok != NULL )
                                   {
                                     strcpy( p_grp_cnt_lab, tok);
                                     ++*p_num_grps;

                                     if (*p_num_grps > 1)
                                       {
                                         printf("ERROR: More than two group commands (group_centre or angle) given.\n");
                                         printf("       Program currently restricted to two groups at most.\n");
                                         exit(0);
                                       }
                                     else if (*p_num_grps == 0)
                                       {
                                         p_groups->group_type1 = CENTRE_TYPE;
                                         strcpy( p_grp_cnt_lab, tok);
                                       }
                                     else if (*p_num_grps == 1)
                                       {
                                         p_groups->group_type2 = CENTRE_TYPE;
                                         strcpy( p_grp_cnt_lab2, tok);
                                       }
                                     else
                                       {
                                         printf("ERROR: Unexpected value for num_grps in read_input.c.\n");
                                         printf("       num_grps = %d\n", *p_num_grps);
                                         exit(0);
                                       }
 
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_grp_cnt_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_grp_cnt_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
        /*************************/
        /** Check for second *****/
        /** group centre.    *****/
        /*************************/
                                     tok=tok_get( input_fp, skip, leave_case );

                                     if (tok != NULL )
                                       {
                                         ++*p_num_grps;

                                     if (*p_num_grps > 1)
                                       {
                                         printf("ERROR: More than two group comands (group_centre or angle) given.\n");
                                         printf("       Program currently restricted to at most two groups.\n");
                                         exit(0);
                                       }
                                     else if (*p_num_grps == 0)
                                       {
                                         p_groups->group_type1 = CENTRE_TYPE;
                                         strcpy( p_grp_cnt_lab, tok);
                                       }
                                     else if (*p_num_grps == 1)
                                       {
                                         p_groups->group_type2 = CENTRE_TYPE;
                                         strcpy( p_grp_cnt_lab2, tok);
                                       }
                                     else
                                       {
                                         printf("ERROR: Unexpected value for num_grps in read_input.c.\n");
                                         printf("       num_grps = %d\n", *p_num_grps);
                                         exit(0);
                                       }
 

                                         *(p_grp_cnt_lab2+4)='\0';
                                         got_last=FALSE;
                                         p_this_char = p_grp_cnt_lab2;
                                         for (iloop=0; iloop<4; iloop++)
                                            {
                                              if ( *p_this_char =='\0') got_last = TRUE;
                                              if ( got_last )  *p_this_char = ' ';
                                              p_this_char++;
                                            }
                                       }
                                   }
                                else
                                   {
                                     printf ("ERROR : No atom label given for group centre\n");
                                   }
                                break;

              case MOL_CENTRE: tok=tok_get( input_fp, skip, leave_case );
                               *p_have_mol = TRUE;
                                if (tok != NULL )
                                   {
                                     strcpy( p_mol_cnt_lab, tok);
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_mol_cnt_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_mol_cnt_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
                                     printf("Read mol centre label >>%s<<\n",p_mol_cnt_lab);
 
                                   }
                                else
                                   {
                                     printf ("ERROR : No atom label given for molecule centre\n");
                                     exit(0);
                                   }
                                break;

              case LATE_CENTRE: tok=tok_get( input_fp, skip, leave_case );
                                *p_have_mol = TRUE;
                                *p_need_late = TRUE;
                                if (tok != NULL )
                                   {
                                     strcpy( p_mol_cnt_lab, tok);
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_mol_cnt_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_mol_cnt_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
                 printf("Read mol centre label for late interpolation >>%s<<\n",
                                                                   p_mol_cnt_lab);
 
                                   }
                                else
                                   {
                                     printf ("ERROR : No atom label given for molecule centre for late interpolation\n");
                                     exit(0);
                                   }
                                break;

              case ANGLE_INTER: tok=tok_get( input_fp, skip, leave_case );
                                *p_have_grp   = TRUE;
                                *p_need_angle = TRUE;

                                if (tok != NULL )
                                   {
                                     strcpy( p_axis1_lab, tok);
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_axis1_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_axis1_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
                 printf("Read first axis definition label as >>%s<<\n",
                                                                   p_axis1_lab);
                                   }
                                else
                                   {
                                     printf ("ERROR : No atom labels given to define axis for angle interpolation\n");
                                     exit(0);
                                   }
       /****************************/
       /* need 2 atom labels for ***/
       /* axis definition.       ***/
       /****************************/
                                tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_axis2_lab, tok);
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_axis2_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_axis2_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
                 printf("Read second axis definition atom >>%s<<\n", p_axis2_lab);
                                   }
                                else
                                   {
                                     printf ("ERROR : Second atom label for definging axis for angle interpolation\n");
                                     printf ("        has not been supplied.\n");
                                     exit(0);
                                   }
                                break;

              case MASTER_FILE: tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_master_input, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for master co-ordinates\n");
                                   }

                                if (strstr(p_master_input, ".car") != NULL)
                                   {
                                     printf("You input a car file string\n");
                                     *p_is_car = TRUE;
                                     *p_is_gulp= FALSE;
                                     *p_is_vasp= FALSE;
                                     *p_is_siesta= FALSE;
                                     *p_have_out= FALSE;
                                   }
                                else if (strstr(p_master_input, ".arc") != NULL)
                                   {
                                     printf("You input an arc file string, will treat as single frame car file\n");
                                     *p_is_car = TRUE;
                                     *p_is_gulp= FALSE;
                                     *p_is_vasp= FALSE;
                                     *p_is_siesta= FALSE;
                                     *p_have_out= FALSE;
                                   }
                                else if (strstr(p_master_input, ".glp") != NULL)
                                   {
                                     printf("You input a gulp file string\n");
                                     *p_is_car = FALSE;
                                     *p_is_gulp= TRUE;
                                     *p_is_vasp= FALSE;
                                     *p_is_siesta= FALSE;
                                     *p_have_out= FALSE;
                                   }
                                else if (strstr(p_master_input, "OUT") != NULL)
                                   {
                                     printf("You input a VASP OUTCAR as master file\n");
                                     *p_is_car = FALSE;
                                     *p_is_gulp= FALSE;
                                     *p_is_vasp= FALSE;
                                     *p_is_siesta= FALSE;
                                     *p_have_out= TRUE;
                                   }
                                else if (strstr(p_master_input, "POS") != NULL || strstr(p_master_input, "CONT") != NULL)
                                   {
                                     printf("You input a VASP POSCAR file format ( POSCAR or CONTCAR ) as the main file\n");
                                     *p_is_car = FALSE;
                                     *p_is_gulp= FALSE;
                                     *p_is_vasp= TRUE;
                                     *p_is_siesta= FALSE;
                                     *p_have_out= FALSE;
                                   }
                                else if (strstr( p_master_input, ".fdf") != NULL )
                                   {
                                     printf("You input a SIESTA fdf file format as the main file\n");
                                     *p_is_car = FALSE;
                                     *p_is_gulp= FALSE;
                                     *p_is_vasp= FALSE;
                                     *p_is_siesta= TRUE;
                                     *p_have_out= FALSE;
                                   }
                                else
                                   {
              printf("File not recognised as .car (.arc), .glp,  POSCAR or CONTCAR in string >>%s<<\n",
                                            p_master_input);
                                     exit(0);
                                   }
                                break;

              case COMP_MODES: *p_compare_modes=TRUE;
                               tok=tok_get( input_fp, skip, leave_case );

                               if (tok != NULL )
                                 {
                                   strcpy( p_mode_compare_input, tok);
                                 }
                               else
                                 {
                                   printf ("ERROR : No file name given for reference modes with compare_modes key.");
                                 }
                               break;

              case  IMAGES: *p_num_images = get_integer( input_fp, skip, &error );

                            if (error)
                             {
                               printf("ERROR: images directive given without specifying");
                               printf(" number of image files to read\n");
                               exit(0);
                             }
                            printf("Expecting %d image files\n",*p_num_images);

                            if ( *p_num_images > MAXIMAGES)
                             {
                               printf("ERROR: images directive states %d images when current maximum is %d\n",
                                      *p_num_images, MAXIMAGES);
                               exit(0);
                             }
                
                            for (iloop=0; iloop < *p_num_images; iloop++)
                             {
                                tok=tok_get( input_fp, TRUE, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_image_files->name, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for image file %d\n", iloop+1);
                                     exit(0);
                                   }

                                if (strstr(p_image_files->name, ".car") == NULL)
                                   {
                                     printf("ERROR: You must input car files for images\n");
                                     printf("ERROR: image %d does not comply.\n", iloop+1);
                                     exit(0);
                                   }
                                printf("Read >>%s<< for image %d\n", p_image_files->name, iloop+1);
                                p_image_files++;
                             }
                            break;

              case  OUTCARS: *p_num_outcars = get_integer( input_fp, skip, &error );

                            if (error)
                             {
                               printf("ERROR: multiple outcars directive given without specifying");
                               printf(" number to read\n");
                               exit(0);
                             }
                            printf("Expecting %d outcar files\n",*p_num_outcars);

                            if ( *p_num_outcars > MAXOUTCARS)
                             {
                               printf("ERROR: outcars directive states %d outcars when current maximum is %d\n",
                                      *p_num_outcars, MAXOUTCARS);
                               exit(0);
                             }
                
                            for (iloop=0; iloop < *p_num_outcars; iloop++)
                             {
                                tok=tok_get( input_fp, TRUE, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_outcar_files->name, tok);
                                     *p_have_out= TRUE;
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for outcar file %d\n", iloop+1);
                                     exit(0);
                                   }

                                printf("Read >>%s<< for outcar %d\n", p_outcar_files->name, iloop+1);
                                p_outcar_files++;
                             }

                            break;
/***************************************************************/
/*** read name of endpoint file ********************************/
/***************************************************************/
              case    END_FILE: tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     strcpy( p_end_input, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for end point co-ordinates\n");
                                   }

                                if (strstr(p_end_input, ".car") != NULL)
                                   {
                                     printf("You input a car file string for the end point\n");
                                     *p_end_is_car = TRUE;
                                     *p_end_is_gulp= FALSE;
                                     *p_end_is_vasp= FALSE;
                                   }
                                else if (strstr(p_end_input, ".glp") != NULL)
                                   {
                                     printf("You input a gulp file string for the end point\n");
                                     *p_end_is_car = FALSE;
                                     *p_end_is_gulp= TRUE;
                                     *p_end_is_vasp= FALSE;
                                   }
                                else if (strstr(p_end_input, "POS") != NULL)
                                   {
                          printf("You input a VASP POSCAR file as the end point file\n");
                                     *p_end_is_car = FALSE;
                                     *p_end_is_gulp= FALSE;
                                     *p_end_is_vasp= TRUE;
                                   }
                                else
                                   {
            printf("End point file not recognised as .car, .glp or POSCAR in string >>%s<<\n",
                                                                 p_end_input);
                                     exit(0);
                                   }
                                break;

              case MINIMAGE_END: { tok=tok_get( input_fp, skip, lower_case );
                                   printf("Found >>%s<< for min_image_end directive\n", tok);

                                    if (tok != NULL)
                                      {
                                        if (strcmp(tok,"yes") == 0)  *p_end_min_image= TRUE;
                                        if (strcmp(tok,"no") == 0)  *p_end_min_image= FALSE;
                                      }
 
                                    else
                                      {
                           printf("ERROR: No or unrecognised additional word with end_min_image directive.\n");
                           printf("       yes or no?\n");
                           exit(0);
                                      }
                                 }
                                break;

              case MD_RUN: printf("Recognised MD_RUN directive  \n");
                           *p_need_md_run = TRUE;
                           break;

              case INTERPOLATE: printf("Recognised interpolation flag\n");
                         *p_num_inter = get_integer( input_fp, skip, &error );
                          printf("Read request for %d interpolated structures\n", *p_num_inter);
                          break;

              case SUPER: printf("Recognised super cell flag\n");
                         *p_super = get_integer( input_fp, skip, &error );
                         *(p_super+1) = get_integer( input_fp, skip, &error );
                         *(p_super+2) = get_integer( input_fp, skip, &error );
                         
                         if (error)
                           {
                              printf("ERROR : super_cell directive given without sufficient values.\n");
                              printf("\n        use syntax : super_cell i j k\n\n");
                              printf("where i j k are integers representing the extent of the\n");
                              printf("super cell require in the a b and c directions.\n");
                              exit(0);
                           }
                          break;

              case NUM_STEPS: printf("Recognised steps flag\n");
                         *p_num_steps = get_integer( input_fp, skip, &error );
                          break;

              case SWITCH: printf("Recognised Switch flag\n");
                         *p_switch = get_double( input_fp, skip, &error );
                         *p_mag_after_switch = get_double( input_fp, skip, &error );
                          break;

              case AMPLITUDE: printf("Recognised Amplitude flag\n");
                         *p_amplitude = get_double( input_fp, skip, &error );
                          break;

              case SHIFT: *p_need_shift=TRUE;
                               break;

              case MORPH: *p_need_morph=TRUE;
                          break;

              case LINEAR: *p_linear=TRUE;
                          break;

	      case NEEDED: tok=tok_get( input_fp, skip, lower_case );
			   printf("second token >>%s<<\n", tok);
                         token = find_kind(tok, SECONDARY_DIRECTIVE);
                         switch (token)
                            {
			    case CAR_FILE_NEEDED : *p_need_car = TRUE;
                                   if (!(tok=tok_get( input_fp, skip, leave_case )))
                                       {
                                          printf("ERROR: No file name given for requested car file output\n");
                                          exit(0);
                                       }
                                   strcpy(p_car_output,tok);
                                   break;

			    case ARC_FILE_NEEDED : *p_need_arc = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_arc_output,tok);
                                   break;

			    case POSCAR_FILE_NEEDED : *p_need_poscar = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_poscar_output,tok);
				   break;

			    case PDB_FILE_NEEDED : *p_need_pdb = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_pdb_output,tok);
				   break;

                            case GULP_FILE_NEEDED : *p_need_gulp = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_gulp_output,tok);
				   break;

                            case DOS_FILE_NEEDED : *p_need_dos = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_dos_output,tok);
				   break;

                            case PART_DOS_NEEDED : *p_need_partdos = TRUE;
                                            *p_num_atoms_pdos=-1;
                                            printf("Reading atom indices for pdos:\n");
                                            while ( tok=tok_get( input_fp, noskip, leave_case ))
                                               {
                                     if (strchr(tok,'-'))
                                       {
/**** Required string filling with \0 to ensure correct termination on ZHENG machine ****/
/**** Dave Willock July 06 **************************************************************/
                                          for (iii=0; iii<5; iii++) cnum1[iii]='\0';
                                          dash=strchr(tok,'-');
                                          printf("Found dash\n");
                                          printf("%d\n", dash-tok);
                                          strncpy(cnum1,tok,dash-tok);
                                          printf("First bit >>%s<<\n",cnum1);
                                          if (atoi(cnum1))
                                            {
                                              srang=atoi(cnum1);
                                            }
                                          else
                                            {
                                              printf("ERROR: Lower range end given for PDOS is not a number\n");
                                              exit(0);
                                            }

                                          dash++;
                                          for (iii=0; iii<5; iii++) cnum1[iii]='\0';
                                          strcpy(cnum2,dash);
                                          printf("Sec.  bit >>%s<<\n",cnum2);
                                          if (atoi(cnum2))
                                            {
                                              erang=atoi(cnum2);
                                            }
                                          else
                                            {
                                              printf("ERROR: Upper range end given for PDOS is not a number\n");
                                              exit(0);
                                            }

                                          if (erang < srang)
                                           {
                                             printf("ERROR: Range given for PDOS has lower bound greater than upper\n");
                                             exit(0);
                                           }

                                          for (iloop=srang; iloop <=erang; iloop++)
                                           {
                                             (*p_num_atoms_pdos)++;
                                             *p_part_dos_list= iloop;
                                             printf("Atom %d added to list\n", *p_part_dos_list );
                                             p_part_dos_list++;
                                           }
                                           
                                       }
                                             else if (atoi(tok) > 0)
                                                   {
                                                     (*p_num_atoms_pdos)++;
                                                      *p_part_dos_list= atoi(tok);
                                                       printf("Atom %d added to list\n", *p_part_dos_list );
                                                       p_part_dos_list++;
                                                   }
                                                 else
                                                   {
                                                     *p_spd = FALSE;
                                                     *(p_spd+1) = FALSE;
                                                     *(p_spd+2) = FALSE;

                      /*** tok is not a number so test if it is the spd flag ***/
                                                     if (strchr(tok,'s')) *p_spd = TRUE; 
                                                     if (strchr(tok,'p')) *(p_spd+1) = TRUE; 
                                                     if (strchr(tok,'d')) *(p_spd+2) = TRUE; 
                                                   }
                                               }
                                   break;

                            case FORCES : *p_need_force= TRUE;
                                   break;

                           case ENERGY : *p_need_energy= TRUE;
                                   if (!(tok=tok_get( input_fp, skip, leave_case )))
                                       {
                                          printf("ERROR: No file name given for requested energy file output\n");
                                          exit(0);
                                       }
                                   strcpy(p_energy_output,tok);
                                   break;


                            case FREQUENCY : *p_need_freq= TRUE;
                                   break;

                            case SHELLS : *p_need_shells= TRUE;
                                          *p_num_shell_species=-1;
                                          printf("Reading Shells:\n");
                                          while ( tok=tok_get( input_fp, noskip, leave_case ))
                                             {
                                               (*p_num_shell_species)++;
                                               strcpy(p_shell_species->label, tok);
                                               printf("%s\n",p_shell_species->label);
                                               p_shell_species++;
                                             }
                                   break;
                            }
			 break;

              case MILLER: printf("Trying to read Miller indices\n");
                           *p_have_miller=TRUE;
                           *p_miller= get_integer( input_fp, noskip, &error );
                           printf("Got Miller: %d\n", *p_miller);
                           p_miller++;
                           *p_miller= get_integer( input_fp, noskip, &error );
                           printf("Got Miller: %d\n", *p_miller);
                           p_miller++;
                           *p_miller= get_integer( input_fp, noskip, &error );
                           printf("Got Miller: %d\n", *p_miller);
                           break;

              case MODE: *p_mode= get_integer( input_fp, skip, &error );
                         break;

              case VARIABLE_SITES: tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_variable_label,tok);
                                   break;

              case MIN_WEIGHT: *p_min_weight= get_double( input_fp, skip, &error );
                                break;
        
              case ASSESS: *p_assess= TRUE; break;

              case ANALYSE: *p_analyse= TRUE; break;


              case NUM_PER_FORMULA_UNIT: *p_num_per_formula= get_integer( input_fp, skip, &error );
                                         break;

            }
        }
      skip= TRUE;
    }
 return 0;
 }
