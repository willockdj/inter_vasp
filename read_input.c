/**********************************************************************************/
/**** read_input does just that, it obtains the input file names ******************/
/**** and the control parameters. Car files are read by read_car ******************/
/**** Adapted May 03 Dave Willock *************************************************/
/**********************************************************************************/

#include <stdio.h>
#include <sys/file.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
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
                char *p_cif_output )
 {
   int iloop, icat, skip, noskip, lower_case, leave_case;
   int end_of_input, itsanum, error, iii, found_doscar;
   int itmp;
   long int list_start, list_end;
   int  iset, token, second_token, third_token;
   int first, got_last, srang, erang;
   char  *tok, *last_tok, *p_this_char;
   char *mon_tok1, *mon_tok2, *mon_tok3, *mon_tok4;
   char *dash, cnum1[5], cnum2[5];

   char_list *p_this_match;

   *p_end_min_image = FALSE;
   end_of_input= FALSE;
   found_doscar= FALSE;
   skip= TRUE;
   noskip= FALSE;
   lower_case= TRUE;
   leave_case= FALSE;

   *p_num_miller=0;
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
                                p_have->potcar=TRUE;

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
                                     p_have->incar=TRUE;
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
                                     if (p_is->siesta_dos)
                                      {
                                        printf("ERROR: Multiple DOS files not supported for SIESTA\n");
                                        exit(0);
                                       }
                                    ++(*p_num_dos_files);
                   /**** This is a secondary occurance so must have multiple DOSCAR files ****/

                                    p_need->multi_dos= TRUE; 

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
                                               p_is->siesta_dos= TRUE;
                                               p_is->vasp_dos= FALSE;
                                              }
                                          else
                                              {
                                               p_is->siesta_dos= FALSE;
                                               p_is->vasp_dos= TRUE;

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
                                     p_need->mdtraj = TRUE;
                                     strcpy( p_mdtraj_input, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for MD trajectory file\n");
                                     exit(0);
                                   }

                                break;

              case REPORT_FILE : tok=tok_get( input_fp, skip, leave_case );

                                if (tok != NULL )
                                   {
                                     p_have->report = TRUE;
                                     strcpy( p_report_file, tok);
                                   }
                                else
                                   {
                                     printf ("ERROR : No file name given for MD REPORT file but report directive in input\n");
                                     exit(0);
                                   }
   
                                printf("Have report file : %s\n", p_report_file);
                                break;

              case DOS_SMEAR : *p_dos_smear = get_double( input_fp, skip, &error );
                               break;

              case GROUP_CENTRE: tok=tok_get( input_fp, skip, leave_case );
                                 p_have->grp = TRUE;
                                if (tok != NULL )
                                   {
                                     ++*p_num_groups;

                                     if (*p_num_groups > 1)
                                       {
                                         printf("ERROR: More than two group commands (group_centre or angle) given.\n");
                                         printf("       Program currently restricted to two groups at most.\n");
                                         exit(0);
                                       }
                                     else if (*p_num_groups == 0)
                                       {
                                         p_groups->group_type1 = CENTRE_TYPE;
                                         strcpy( p_groups->cnt_lab1, tok);
                                       }
                                     else if (*p_num_groups == 1)
                                       {
                                         p_groups->group_type2 = CENTRE_TYPE;
                                         strcpy( p_groups->cnt_lab2, tok);

                                         *(p_groups->cnt_lab2+4)='\0';
                                         got_last=FALSE;
                                         p_this_char = p_groups->cnt_lab2;
                                         for (iloop=0; iloop<4; iloop++)
                                            {
                                              if ( *p_this_char =='\0') got_last = TRUE;
                                              if ( got_last )  *p_this_char = ' ';
                                              p_this_char++;
                                            }

                                       }
                                     else
                                       {
                                         printf("ERROR: Unexpected value for num_groups in read_input.c.\n");
                                         printf("       num_groups = %d\n", *p_num_groups);
                                         exit(0);
                                       }
 
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_groups->cnt_lab1+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_groups->cnt_lab1;
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
                                         ++*p_num_groups;

                                     if (*p_num_groups > 1)
                                       {
                                         printf("ERROR: More than two group comands (group_centre or angle) given.\n");
                                         printf("       Program currently restricted to at most two groups.\n");
                                         exit(0);
                                       }
                                     else if (*p_num_groups == 0)
                                       {
                                         p_groups->group_type1 = CENTRE_TYPE;
                                         strcpy( p_groups->cnt_lab1, tok);
                                       }
                                     else if (*p_num_groups == 1)
                                       {
                                         p_groups->group_type2 = CENTRE_TYPE;
                                         strcpy( p_groups->cnt_lab2, tok);
                                       }
                                     else
                                       {
                                         printf("ERROR: Unexpected value for num_groups in read_input.c.\n");
                                         printf("       num_groups = %d\n", *p_num_groups);
                                         exit(0);
                                       }
 

                                         *(p_groups->cnt_lab2+4)='\0';
                                         got_last=FALSE;
                                         p_this_char = p_groups->cnt_lab2;
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
                               p_have->mol = TRUE;
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
                                p_have->mol = TRUE;
                                p_need->late = TRUE;
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
                                p_have->grp   = TRUE;
                                p_need->angle = TRUE;

                                if (tok != NULL )
                                   {
                                     strcpy( p_groups->axis1_lab, tok);
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_groups->axis1_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_groups->axis1_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
                 printf("Read first axis definition label as >>%s<<\n",
                                                                   p_groups->axis1_lab);
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
                                     strcpy( p_groups->axis2_lab, tok);
        /*************************/
        /* force 4 letter label **/
        /*************************/
                                     *(p_groups->axis2_lab+4)='\0';
                                     got_last=FALSE;
                                     p_this_char = p_groups->axis2_lab;
                                     for (iloop=0; iloop<4; iloop++)
                                        {
                                          if ( *p_this_char =='\0') got_last = TRUE;
                                          if ( got_last )  *p_this_char = ' ';
                                          p_this_char++;
                                        }
                 printf("Read second axis definition atom >>%s<<\n", p_groups->axis2_lab);
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
                                     p_is->car = TRUE;
                                   }
                                else if (strstr(p_master_input, ".cif") != NULL)
                                   {
                                     printf("You input a cif file string\n");
                                     p_is->cif = TRUE;
                                   }
                                else if (strstr(p_master_input, ".arc") != NULL)
                                   {
                                     printf("You input an arc file string, will treat as single frame car file\n");
                                     p_is->car = TRUE;
                                   }
                                else if (strstr(p_master_input, ".glp") != NULL)
                                   {
                                     printf("You input a gulp file string\n");
                                     p_is->gulp= TRUE;
                                   }
                                else if (strstr(p_master_input, "OUT") != NULL)
                                   {
                                     printf("You input a VASP OUTCAR as master file\n");
                                     p_have->out= TRUE;
                                   }
                                else if (   strstr(p_master_input, "POS") != NULL || strstr(p_master_input, "CONT") != NULL
                                         || strstr(p_master_input, ".vasp") != NULL )
                                   {
                                     printf("You input a VASP POSCAR file format ( POSCAR or CONTCAR ) as the main file\n");
                                     p_is->vasp= TRUE;
                                   }
                                else if (strstr( p_master_input, ".fdf") != NULL )
                                   {
                                     printf("You input a SIESTA fdf file format as the main file\n");
                                     p_is->siesta= TRUE;
                                   }
                                else if (strstr( p_master_input, ".dat") != NULL )
                                   {
                                     printf("You input a ONETEP dat file format as the main file\n");
                                     p_is->onetep= TRUE;
                                   }
                                else if (strstr( p_master_input, ".pun") != NULL )
                                   {
                                     printf("You input a ChemShell pun file format as the main file\n");
                                     p_is->punch = TRUE;
                                   }
                                else if (strstr( p_master_input, ".pdb") != NULL )
                                   {
                                     printf("You input a protein data bank pdb file format as the main file\n");
                                     p_is->pdb = TRUE;
                                   }
                                else
                                   {
              printf("File not recognised as .car (.arc), .glp,.pun,.cif, POSCAR or CONTCAR in string >>%s<<\n",
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
                                     p_have->out= TRUE;
                                     printf("Setting have_out TRUE from OUTCARS read\n");
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
                                     p_is->end_car = TRUE;
                                   }
                                else if (strstr(p_end_input, ".glp") != NULL)
                                   {
                                     printf("You input a gulp file string for the end point\n");
                                     p_is->end_gulp= TRUE;
                                   }
                                else if (strstr(p_end_input, "POS") != NULL)
                                   {
                          printf("You input a VASP POSCAR file as the end point file\n");
                                     p_is->end_vasp= TRUE;
                                   }
                                else
                                   {
            printf("End point file not recognised as .car, .glp or POSCAR in string >>%s<<\n",
                                                                 p_end_input);
                                     exit(0);
                                   }
                                break;

/*** Look for match directive, asking to match up atoms in neb set ups based on inter-atom distance ***/
/*** rather than atom list position. User should specify elements that are treated in this way      ***/
              case MATCH: { printf("Found Match directive....\n");
                            p_need->match=TRUE;
                            *p_num_elem_to_match = -1;

                            p_this_match = p_elems_to_match;
                            while ( tok=tok_get( input_fp, noskip, leave_case ))
                              {
                                printf("Read element: %s\n", tok);
                                strcpy(p_this_match->name, tok);
                                ++*p_num_elem_to_match;

                                if ( 1+*p_num_elem_to_match > MAX_MATCH_LIST )
                                     printf("ERROR: Too many elements in match list, current maximum is %d\n", MAX_MATCH_LIST);
                                else p_this_match++;
                              }
                            printf("Have %d elements to match in this interpolation...\n", *p_num_elem_to_match);
                          }
                          break;

              case MINIMAGE_END: { tok=tok_get( input_fp, skip, lower_case );
                                   printf("Found >>%s<< for min_image_end directive\n", tok);

                                    iset = FALSE;
                                    if (tok != NULL)
                                      {
                                        if (strcmp(tok,"yes") == 0 || strcmp(tok,"on") == 0 )  
                                           { iset = TRUE; *p_end_min_image= TRUE; }
 
                                        if (strcmp(tok,"no") == 0 || strcmp(tok,"off") == 0 ) 
                                           { iset= TRUE; *p_end_min_image= FALSE; }
                                      }
                                    if (!iset)
                                      {
                                         printf("ERROR: No or unrecognised additional word with end_min_image directive.\n");
                                         printf("       yes or no?\n");
                                         exit(0);
                                      }
                                 }
                                break;

              case ZSORT : printf("Recognised zsort directive  \n");
                           p_need->zsort  = TRUE;
                           printf("read_input need_zsort = %d\n", p_need->zsort);
                           break;

              case SLAB  : tok=tok_get( input_fp, skip, lower_case );
                           if (!tok)
                             {
                           printf("ERROR: Slab directive given but no second token, expecting oshift or vgap\n");
                             }
                           printf("Slab: second token >>%s<<\n", tok);

                           token = find_kind(tok, SECONDARY_DIRECTIVE);
                           switch (token)
                            {
                               case VGAP: printf("Recognised vgap directive looking for value of vacuum gap to use\n");
                                          p_need->vgap = TRUE;
                                          *p_vgap = get_double( input_fp, skip, &error );
                                          if (error)
                                             {
                                               printf("ERROR: slab vacuum gap requested but no value supplied...\n");
                                               exit(0); 
                                             }
                                          printf("Will add vacuum gap of %10.6f Angstroms\n", *p_vgap );
                                          break;

                               case ORIGIN_SHIFT : printf("Recognised origin shift directive  \n");
                                                   p_need->oshift  = TRUE;
                                                   *p_oshift = get_double( input_fp, skip, &error );
                                                   *(p_oshift+1) = get_double( input_fp, skip, &error );
                                                   *(p_oshift+2) = get_double( input_fp, skip, &error );
                                                   if (error)
                                                     {
                                                       printf("ERROR: slab origin shift requested but no vector supplied...\n");
                                                       exit(0); 
                                                     }
                                                   printf("Will shift atoms by = %10.6f %10.6f %10.6f\n", *p_oshift, *(p_oshift+1), *(p_oshift+2) );
                                                   break;
                            }
                          break;

              case MD_RUN: printf("Recognised MD_RUN directive  \n");
                           p_need->md_run = TRUE;
                           break;

              case INTERPOLATE: printf("Recognised interpolation flag\n");
                         *p_num_inter = get_integer( input_fp, skip, &error );
                          printf("Read request for %d interpolated structures\n", *p_num_inter);
                          break;

              case MONITOR: printf("Recognised monitor flag\n");
                         tok=tok_get( input_fp, FALSE, lower_case );
                         printf("Next tok= %s %d\n",tok, atoi(tok));
                         p_need->monit=TRUE;

                         if (atoi(tok) == 0)
                           {
                              printf("Found string after monitor so need to code new function....\n");
/*** pick up the sequence of tokens for this command ***/
/** Have defined new toks for the monitor command char *mon_tok1, *mon_tok2, *mon_tok3, *mon_tok4; **/
/** FALSE means don't jump lines, lower_case means convert input to lower_case tok=tok_get( input_fp, FALSE, lower_case ); ***/

                              mon_tok1=tok;
                              mon_tok2=tok_get( input_fp, FALSE, leave_case );
                              mon_tok3=tok_get( input_fp, FALSE, lower_case );
                              mon_tok4=tok_get( input_fp, FALSE, leave_case );
                              
                              printf("Read the toks: %s %s %s %s\n", mon_tok1, mon_tok2, mon_tok3, mon_tok4);

/*** Check that all mon_toks are valid strings, i.e. not "NULL" ***/
                              if (mon_tok1 == NULL || mon_tok2 == NULL || mon_tok3 == NULL || mon_tok4 == NULL  )
                                {
                                   printf("ERROR: Not enough arguements with monitor directive\n");
                                   printf("       should specify type and required element, number or label\n");
                                   printf("       eg: monitor num 4 elem H \n");
                                   exit(0);
                                }
 

/*** Recognise the cases of monit_tok1 is "elem", "label" or "num" and complain if none of these ***/

   
/** *p_monit_type1, *p_monit_type2 have values -1: not set 1: expect number; 2: expect element 3: expect label; */
                              if (mon_tok1 != NULL)
                                {
                                  if (strcmp(mon_tok1,"num") == 0 ) 
                                    {
                                      *p_monit_type1 =1;  
                                      *p_monit_at1 = atoi(mon_tok2)-1;
                                    }
                                  else if (strcmp(mon_tok1,"elem") == 0 ) *p_monit_type1=2; 
                                  else if (strcmp(mon_tok1,"label") == 0 ) *p_monit_type1=3; 
                                  else
                                    {
                                       printf("ERROR: Unrecognised arguments inputted for first set to monitor.\n");
                                       printf("      Make sure they are in the format atom num, elem or label\n");
                                       printf("      eg: monitor num 4 elem H \n");
                                       exit(0);
                                    }
                                }

                              if (*p_monit_type1 > 1) strcpy(p_mon_str1, mon_tok2);
                                      
                              if (mon_tok3 != NULL)
                                   tok=tok_get( input_fp, skip, leave_case );
                                {
                                  if (strcmp(mon_tok3,"num") == 0 )
                                    {
                                       *p_monit_type2=1;  
                                       *p_monit_at2 = atoi(mon_tok4)-1;
                                    }
                                  else if (strcmp(mon_tok3,"elem") == 0 ) *p_monit_type2=2; 
                                  else if (strcmp(mon_tok3,"label") == 0 ) *p_monit_type2=3; 
                                  else 
                                    {
                                       printf("ERROR: Unrecognised arguments inputted for second set to monitor.\n");
                                       printf("      Make sure they are in the format atom num, elem or label\n");
                                       printf("      eg: monitor num 4 elem H \n");
                                       exit(0);
                                    }
                                }

                              if (*p_monit_type2 > 1) strcpy(p_mon_str2, mon_tok4);
                           }
                         else
                           {
                              *p_monit_type1=1; *p_monit_type2=1;  
                              *p_monit_at1 = atoi(tok);
                              *p_monit_at2 = get_integer( input_fp, skip, &error );
/*** Check input is sensible ****/
                              if (*p_monit_at1 < 1 ||*p_monit_at2 < 1 ) 
                                {
                                  printf("ERROR: monitor indices for atoms must be positive\n");
                                  exit(0);
                                }
                              if (*p_monit_at1 == *p_monit_at2 ) 
                                {
                                  printf("ERROR: monitor indices for atoms cannot be the same\n");
                                  exit(0);
                                }
/*** order monit atoms so that lowest index is always first **/
                              if ( *p_monit_at2 < *p_monit_at1 )
                                {
                                   itmp=*p_monit_at2;
                                   *p_monit_at2=*p_monit_at1;
                                   *p_monit_at1=itmp;
                                }
                               --*p_monit_at1;
                               --*p_monit_at2;
                               printf("Read atoms to monitor as %d and %d\n", 1+*p_monit_at1,  
                                                                              1+*p_monit_at2);
                             }
                            
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
                         *p_num_after_switch = get_integer( input_fp, skip, &error );

                          if ( *p_num_after_switch < 1 )
                            {
                              printf("ERROR: Number of steps after the switch must be a positve integer\n");
                              exit(0);
                            }
                          break;

              case AMPLITUDE: printf("Recognised Amplitude flag\n");
                         *p_amplitude = get_double( input_fp, skip, &error );
                          break;

              case SHIFT: p_need->shift=TRUE;
                               break;

              case MORPH: p_need->morph=TRUE;
                          break;

              case LINEAR: *p_linear=TRUE;
                          break;

              case HBOND : p_need->hbond=TRUE;
                          break;

              case PERMUTE : tok=tok_get( input_fp, skip, lower_case );
                             p_need->permute=TRUE;
                             printf("read_input: need.permute = %d TRUE = %d FALSE = %d\n", p_need->permute, TRUE, FALSE);
                             printf("PERMUTE: second token >>%s<<\n", tok);
                             token = find_kind(tok, SECONDARY_DIRECTIVE);
                             switch (token)
                            {
                               case ELEMENT: if (!(tok=tok_get( input_fp, skip, leave_case )))
                                                {
                                                   printf("ERROR: No element type given for requested permute element directive.\n");
                                                   exit(0);
                                                }
                                             strcpy(p_permute->elem, tok);
                                             printf("Read PERMUTE element %s\n", p_permute->elem );
                                             break;
  
                               case MAX_DIST: p_permute->maxd= get_double( input_fp, skip, &error );
                                              printf("Read PERMUTE max distance %10.6f\n", p_permute->maxd );
                                              break;
  
                               case MIN_DIST: p_permute->mind= get_double( input_fp, skip, &error );
                                              printf("Read PERMUTE min distance %10.6f\n", p_permute->mind );
                                              break;
  
                               case NUMBER  : p_permute->num = get_integer( input_fp, skip, &error );
                                              printf("Read PERMUTE number %d\n", p_permute->num );
                                              break;
  
                               case CHECK   : p_permute->check = TRUE;
                                              printf("Read PERMUTE check flag\n");
                                              break;
  
                               case CENTRE  : printf("found Centre\n");
                                              p_permute->centre = get_integer( input_fp, skip, &error );
                                              p_have->perm_centre = TRUE;
                                              printf("Read PERMUTE centre %d\n", p_permute->centre );
                                              break;
  
                               case SUBS_ELEMENT : if (!(tok=tok_get( input_fp, skip, leave_case )))
                                                {
                                                   printf("ERROR: No substitute element type given for requested permute substitute element directive.\n");
                                                   exit(0);
                                                }
                                             strcpy(p_permute->subs_elem, tok);
                                             printf("Read PERMUTE substitute element %s\n",
                                                                                p_permute->subs_elem );
                                             p_have->perm_subs_elem = TRUE;
                                             break;
  
                               case DEBUG   : p_permute->debug = TRUE;
                                              printf("Read PERMUTE debug flag\n");
                                              break;
                            }
                            break;


	      case NEEDED: tok=tok_get( input_fp, skip, lower_case );
			   printf("second token >>%s<<\n", tok);
                         token = find_kind(tok, SECONDARY_DIRECTIVE);
                         switch (token)
                            {
			    case CAR_FILE_NEEDED : p_need->car = TRUE;
                                   if (!(tok=tok_get( input_fp, skip, leave_case )))
                                       {
                                          printf("ERROR: No file name given for requested car file output\n");
                                          exit(0);
                                       }
                                   strcpy(p_car_output,tok);
                                   break;

			    case CIF_FILE_NEEDED : p_need->cif = TRUE;
                                   if (!(tok=tok_get( input_fp, skip, leave_case )))
                                       {
                                          printf("ERROR: No file name given for requested cif file output\n");
                                          exit(0);
                                       }
                                   strcpy(p_cif_output,tok);
                                   break;

			    case ARC_FILE_NEEDED : p_need->arc = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_arc_output,tok);
                                   break;

                           case EXPANSION_NEEDED: printf("Recognised expansion flag\n");
                                      p_need->expansion = TRUE;
                                      *p_expansion = get_double( input_fp, skip, &error );
                                      if (error)
                                        {
                                           printf("ERROR : unit cell expansion directive given without sufficient values.\n");
                                           printf("\n        use syntax : need expansion i\n");
                                           printf("where i is the expansion factor (fractional)\n");
                                           printf("reverting to default settings: i = 0.9\n");
                                           *p_expansion = 0.9;
                                        }
                                      else
                                        {
                                           printf("expanding to fraction: %3.1f \n",*p_expansion);
                                        }
                                      break;

			    case POSCAR_FILE_NEEDED : p_need->poscar = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_poscar_output,tok);
                                  
                                   tok=tok_get( input_fp, skip, lower_case );
                                   if  ( tok != NULL )
                                     {
                                        printf("Extra info at end of need poscar line\n");
                                        if ( strncmp( tok, "cart", 4) == 0 ) 
                                          {
                                             p_need->poscar_frac= FALSE;
                                          }
                                        else if  ( strncmp( tok, "dire", 4) == 0 ) 
                                          {
                                             p_need->poscar_frac= TRUE;
                                          }
                                        else if ( strncmp( tok, "frac", 4) == 0 ) 
                                          {
                                             p_need->poscar_frac= TRUE;
                                          }
                                        else
                                          {
                                            printf("ERROR : Additional directive with need poscar can only be\n");
                                            printf("ERROR : cart, fract or direct.                           \n");
                                            exit(0);
                                          }
                                     }
                                   else
                                     {
                                        p_need->poscar_frac= FALSE;
                                     }
				   break;

			    case POTCAR_FILE_NEEDED : p_need->potcar = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_potcar_output,tok);
				   break;

			    case ONETEP_FILE_NEEDED : p_need->onetep = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_onetep_output,tok);
				   break;

			    case PDB_FILE_NEEDED : p_need->pdb = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_pdb_output,tok);
				   break;

                            case GULP_FILE_NEEDED : p_need->gulp = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_gulp_output,tok);
				   break;

                            case DOS_FILE_NEEDED : p_need->dos = TRUE;
                                   tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_dos_output,tok);
				   break;

                            case PART_DOS_NEEDED : p_need->part_dos = TRUE;
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
                                                     *(p_spd+3) = FALSE;

                      /*** tok is not a number so test if it is the spd flag ***/
                                                     if (strchr(tok,'s')) *p_spd = TRUE; 
                                                     if (strchr(tok,'p')) *(p_spd+1) = TRUE; 
                                                     if (strchr(tok,'d')) *(p_spd+2) = TRUE; 
                                                     if (strchr(tok,'f')) *(p_spd+3) = TRUE; 
                                                   }
                                               }
                                   break;

                            case FORCES : p_need->force= TRUE;
                                   break;

                           case ENERGY : p_need->energy= TRUE;
                                   if (!(tok=tok_get( input_fp, skip, leave_case )))
                                       {
                                          printf("ERROR: No file name given for requested energy file output\n");
                                          exit(0);
                                       }
                                   strcpy(p_energy_output,tok);
                                   break;


                            case FREQUENCY : p_need->freq= TRUE;
                                   break;

                            case SHELLS : p_need->shells= TRUE;
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

              case MILLER: ++*p_num_miller;
                           printf("Trying to read Miller indices set %d\n", *p_num_miller);
                           p_have->miller=TRUE;
                           if (*p_num_miller <= MAX_MILLER)
                             {
                               *p_miller= get_integer( input_fp, noskip, &error );
                               printf("Got Miller: %d\n", *p_miller);
                               p_miller++;
                               *p_miller= get_integer( input_fp, noskip, &error );
                               printf("Got Miller: %d\n", *p_miller);
                               p_miller++;
                               *p_miller= get_integer( input_fp, noskip, &error );
                               printf("Got Miller: %d\n", *p_miller);
                               p_miller++;
                             }
                           else
                             {
                               printf("ERROR: Too many Miller planes defined, current maximum = %d\n", MAX_MILLER);
                               exit(0);
                             }
                                  
                           tok=tok_get( input_fp, skip, lower_case );
                           if  ( tok != NULL )
                             {
                                printf("Extra info at end of Miller line: >>%s<<\n", tok);
                                if ( strncmp( tok, "sort", 4) == 0 ) 
                                  {
                                     printf("Would try to sort along this direction\n");
                                     p_need->miller_sort=TRUE;
                                  }
                                else
                                  {
                                    printf("ERROR : Additional directive with need poscar can only be\n");
                                    printf("ERROR : cart, fract or direct.                           \n");
                                    exit(0);
                                  }
                             }
                           else
                             {
                        //        define default
                                  p_need->miller_sort=FALSE;
                             }
                           break;

              case MODE: *p_mode= get_integer( input_fp, skip, &error );
                         break;

              case VARIABLE_SITES: tok=tok_get( input_fp, skip, leave_case );
                                   strcpy(p_variable_label,tok);
                                   break;

              case MIN_WEIGHT: *p_min_weight= get_double( input_fp, skip, &error );
                                break;
        

              case NUM_PER_FORMULA_UNIT: *p_num_per_formula= get_integer( input_fp, skip, &error );
                                         break;

            }
        }
      skip= TRUE;
    }

/********************************************/
/**For permutation process                ***/
/**Error will be shown if centre is present */
/**but didn't specify substitute element. ***/
/********************************************/
    if ( p_have->perm_centre && !p_have->perm_subs_elem )
    {
     printf("ERROR: No substitute element type given for requested permute element directive!\n");
     exit(0);
    }
   printf("read_input returning: need.permute = %d TRUE = %d FALSE = %d\n", p_need->permute, TRUE, FALSE);
 return 0;
 }
