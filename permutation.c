/******************************************************************/
/** Routine to permute the atoms flagged within the permute *******/
/** structure.                                              *******/
/** Begun June 2013, Dave Willock                           *******/
/******************************************************************/
/** p_molecule is a pointer to the start of the atom list   *******/
/** num_atoms is the total number of atoms in the list      *******/
/** p_latt_vec and p_recip_latt are pointers to the start of*******/
/**       3x3 arrays that hold the lattice vectors and      *******/
/**       reciprocal lattice vectors as ax, ay, az etc.     *******/
/** p_abc is a pointer to the start of the array that holds *******/
/**       the lattice vectors in a,b,c, alpha,beta,gamma    *******/
/**       form.                                             *******/
/** p_permute is a pointer to the structure holding the     *******/
/**       permute information. Look in structures.h for     *******/
/**       available members.                                *******/
/******************************************************************/
/******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "maxima.h"
#include "structures.h"
#include "global_values.h"

#define MAX_PERM 50

/* prototype list for this routine */
void *calloc(size_t nitems, size_t size);

//void set_binary(int *p_distrib, int state, int num_sites);
void set_binary(int *p_distrib, unsigned long long int state, int num_sites);

void detail_state(int *p_distrib, int *p_state, int num_sites);

void fill_structure(atom *p_master, atom *p_molecule, int *p_indices, int *p_distrib, 
                    int num_indices, int *p_num_new_atoms, int num_atoms);

void open_file(FILE **p_file, char *p_filename, char *p_status);

void write_car( FILE *fp, int *p_header_line, int *p_title_line, char *p_c_title_line,
		int *p_date_line, atom *p_molecule, int *p_mol_number, int use_mols,
                int pbc, double *p_abc, int num_atoms, double scale_factor, 
                int start_frame, int *p_super, double *p_latt_vec, 
                double *p_recip_latt_vec, coord_flags *p_fix_flags, 
                double *p_magmom, int num_magmom );

void write_poscar( FILE *fp, atom *p_molecule, double *p_fract_coords,
                    atom_number *p_types, int num_types,
                   double *p_latt_vec, double *p_scale_factor, int num_atoms,
                   int *p_title_line, char *p_c_title_line, int pbc, int is_fract,
                   coord_flags *p_fix_flags);

int read_car( FILE *fp, int *p_header_line, int *p_title_line,
              atom *p_molecule, int *p_date_line, int *p_pbc, int *p_num_atoms,
              int *p_num_of_mols, int *p_num_mol_members, int *p_mol_number,
              double *p_abc, int *p_been_before, group_lists *p_groups,
              int have_grp, int *p_num_grps,
              int find_fixed, coord_flags *p_fix_flags);

void generate_neighbours( atom *p_molecule, int num_atoms, 
                          atom_number *p_types, int *p_num_types,
                          int use_pbc,  double *p_recip_latt_vec, double *p_latt_vec,
                          charge_list *p_spec_charges);

void sort_by_elem( atom *p_molecule, int num_atoms, atom_number *p_types,
                   int num_types );

void min_image( double *x, double *y, double *z, double *p_recip_latt_vec, 
	        double *p_latt_vec);

double atom_separation_squared( atom *p_A, atom *p_B, int pbc,
				double *p_recip_latt_vec, double *p_latt_vec);

/*--------------------------------------------------------------------------------------------*/

int permutation(atom *p_master, int *p_num_atoms, double *p_latt_vec, double *p_recip_latt_vec,
                int pbc, double *p_abc, perms *p_permute, double *p_magmom, int num_magmom )
  {
  int num_set, num_set_2, iii, jjj, iatom, jatom, iselected, num_perms=0, num_indices;
  int header_line[LINESIZ], title_line[LINESIZ], date_line[LINESIZ];
  int num_new_atoms, good=0;
  int num_types, iloop, seen_count, size;
  int perm_indices[MAX_PERM], distrib[MAX_PERM];
  int mol_number[MAXATOMS], super[3];
  int ineigh, neigh_index, iiii, centre;
  int have_centre, label_num, have_mind, have_maxd;
  int seen_r, Processes=1, throw, latt_test;
 // double upper, state; 
  unsigned long long int upper, state; 
  

  atom *p_this_atom, temp_molecule[MAXATOMS];
  seen_lists * seen_perm;
  double idist_to_centre, jdist_to_centre, ij_distance, fract_coords[3*MAXATOMS];
  double scale_factor;

  char filename[80], poscar_output[80];
  char command[180];
  char c_header_line[LINESIZ], c_title_line[LINESIZ], c_date_line[LINESIZ];
  FILE *fp, *fp_vasp_output;

  coord_flags fix_flags[MAXATOMS];
  atom_number types[MAXATOMS];
  charge_list spec_charges[MAXATOMS];

/*** Set if a centre, min & max distance are specified in input ***/
  if (p_permute->centre != -1) have_centre=TRUE;
  else have_centre=FALSE;
  if (p_permute->mind != -1) have_mind=TRUE;
  else have_mind=FALSE;
  if (p_permute->maxd != -1) have_maxd=TRUE;
  else have_maxd=FALSE;

/*** Count the number of flagged atoms in the structure ***/

  printf("Starting permutation routine for structure with %d atoms\n\n", *p_num_atoms);

  printf("Permutation controls:\n");
  printf("Element to be considered: %s\n", p_permute->elem);
  printf("Number to replace in place: %d\n", p_permute->num);
  if ( have_mind )
    printf("Minimum allowed distance: %10.6f\n", p_permute->mind);
  if ( have_maxd )
    printf("Maximum allowed distance: %10.6f\n", p_permute->maxd);
  centre=p_permute->centre-1;

  
  if (have_centre)
  {
   printf("This will consider the distances of atoms to selected centre: %s , index : %d\n", 
									(p_master+centre)->label, 
									centre );
   Processes=0;
   if ( p_permute->debug ) printf("Process = %d\n", Processes);
  }
  if (p_permute->check) printf("Just checking number that would be produced\n");
  else printf("Should try and make structures.\n");
  if (p_permute->debug) printf("Will print debugging statements!\n");
  printf("\n");

  super[0]=1; super[1]=1; super[2]=1;
  num_indices = -1;
  p_this_atom=p_master;
  for (iatom=0; iatom<*p_num_atoms; iatom++)
    {
      printf("%s %10.6f %10.6f %10.6f ", p_this_atom->label,
                                         p_this_atom->x,
                                         p_this_atom->y,
                                         p_this_atom->z );

      mol_number[iatom]=0;

      fix_flags[iatom].fx = FALSE; fix_flags[iatom].fy = FALSE; fix_flags[iatom].fz = FALSE; 

      if (strcmp(p_this_atom->elem, p_permute->elem) == 0)
        {
           printf(" an atom to be considered\n");
           num_indices++;
  
           if (num_indices >= MAX_PERM )
             {
               printf("ERROR: number of indices for permutation exceeds current maximum of %d\n", MAX_PERM);
               exit(0);
             }
           perm_indices[num_indices]=iatom;   
        }  
      else
           printf("\n");

      p_this_atom++;
    }

/**** Do some sanity checks ***/
  if ( p_permute->num > num_indices+1 )
    {
       printf("ERROR: Number set to replace in permutations is greater than number of atoms of the type being replaced\n");
       exit(0);
    }
  if (  have_centre && p_permute->centre > *p_num_atoms-1 )
    {
       printf("ERROR: Index of centre for removal is greater than number of atoms in system.\n");
       exit(0);
    }
  if (p_permute->mind > p_permute->maxd)
    {
       printf("ERROR: Minimum permitted distance is greater than maximum, check input.\n");
       exit(0);
    }

  latt_test= p_permute->maxd > *p_abc/2 || p_permute->maxd > *(p_abc+1)/2 || p_permute->maxd > *(p_abc+2)/2;
  if ( latt_test)
    {
       printf("ERROR: Maximum permitted distance is greater than half smallest lattice vector so min_image will cause errors.\n");
       exit(0);
    }
 

  printf("List of the %d selected atoms for permuting:\n", num_indices+1);
  for (iselected = 0; iselected<=num_indices; iselected++)
    {
      p_this_atom= p_master+perm_indices[iselected];
      printf("%d: %d %s\n",iselected, perm_indices[iselected], p_this_atom->label);
    }

  header_line[0]=-1;
  strcpy(c_header_line, "No Header given");
  date_line[0]=-1;
  strcpy(c_date_line, "No date given");
  title_line[0]=-1;
  strcpy(c_title_line, "No Title given");
  c_title_line[0]=-1;
  strcpy(c_title_line, "No Title given");

// setting upper limit for the total number of states

  if (num_indices < 64)
   {  
     upper = pow(2.0, num_indices+1);
     set_binary(&distrib[0], upper, num_indices);
     printf("distrib array for upper %llu :\n", upper);
     //printf("distrib array for upper %.0f : ", upper);
     for (iii=num_indices; iii>=0; iii--) printf(" %d ",distrib[iii]);
             printf("\n");
     if ( have_centre )
     {
// test the size of memory allocation for seen_perm
//      size = pow(2.0, (num_indices+1)/4);

      size = pow(num_indices+1, 2);
      printf("Number of indices         : %d\n", num_indices);
      printf("So possible seen list size: %d\n", size);

//
// Use a more sensible size DJW Jan 2016
//      size = upper;
      if (p_permute->debug) printf("test the array size is : %d, %lu\n", size, size*sizeof(seen_perm));

      seen_perm=(seen_lists *) calloc(size, size*sizeof(seen_perm));
     }
   }
  else
   {
     printf("ERROR: A little too ambitious to do more than 64 in permutations....\n");
     exit(0);
   }

  printf("Trying %llu numbers.....", upper);  
  seen_count=0;
  for (state = 0; state < upper; state++)
    {
      set_binary(&distrib[0], state, num_indices);

      num_set= 0;
      for ( iii=0; iii <= num_indices; iii++) num_set +=distrib[iii];

      if (num_set == p_permute->num)
        {
	  throw=FALSE;
          printf("\n\ndistrib array for state %llu :\n ", state);
          for (iii=num_indices; iii>=0; iii--) printf(" %d ",distrib[iii]);
          printf("\n");
          num_set_2=0;
/**** Work out inter-atomic distances and flag if below minimum ****/
        //printf("User Defined Minimum Distance %f : \n", p_permute->mind);  
          for (iii=0; iii<=num_indices; iii++)
            {
              if ( distrib[iii] == 1 && num_set_2 < p_permute->num ) 
                {
                   iatom=perm_indices[iii];
	           num_set_2 +=distrib[iii];
		   if (p_permute->debug) printf(">>DEBUG<< iatom = %d \n", iatom);	     
		   if ( have_centre )
		   {
		    idist_to_centre = atom_separation_squared(p_master+centre,  p_master+iatom, pbc,
								   p_recip_latt_vec,
								   p_latt_vec);
                    idist_to_centre=sqrt(idist_to_centre);
                    idist_to_centre=ceil(idist_to_centre * 1000)/1000;
                    printf("distance between atom %d %s and centre is %.2f\n",
                                                                        iatom, (p_master+iatom)->label,
                                                                        idist_to_centre);
                    
 		   }
	
                   for (jjj=iii+1; jjj<=num_indices; jjj++)
                     {
                        if ( distrib[jjj] == 1 ) 
                          {
/**** In here we have two atoms that will be included in the output structure ****/
/**** They have indicies iatom and jatom so should measure distance.......... ****/
                             jatom=perm_indices[jjj]; 
	                     num_set_2 +=distrib[jjj];
		             if (p_permute->debug) printf(">>DEBUG<< jatom = %d \n", jatom);	     
		             if ( have_centre )
		   	     {
		    	      jdist_to_centre = atom_separation_squared(p_master+centre,p_master+jatom, pbc,
								   	    p_recip_latt_vec,
								   	    p_latt_vec);
                    	      jdist_to_centre=sqrt(jdist_to_centre);
                              jdist_to_centre=ceil(jdist_to_centre * 1000)/1000;
                              printf("distance between atom %d %s and centre is %.2f\n",
                                                                        jatom, (p_master+jatom)->label,
                                                                        jdist_to_centre);
                    
 		   	     }
  	                     if ( have_maxd && ( idist_to_centre > p_permute->maxd || jdist_to_centre > p_permute->maxd ))
                             {
	                        printf("Exceeded maximum distance! Discarded!\n");
 				throw=TRUE;				
			     }
 			     else 
 			     {
			      printf("Keeping atom pair %d %s and %d %s\n", iatom, (p_master+iatom)->label,
                                                                           jatom, (p_master+jatom)->label);
			     }
			      ij_distance = atom_separation_squared(p_master+iatom,
								   p_master+jatom, pbc,
								   p_recip_latt_vec,
								   p_latt_vec);

                    	      ij_distance=sqrt(ij_distance);
                              ij_distance=ceil(ij_distance * 1000)/1000;
			     
                    	      printf("distance between atom pair %d %s and %d %s is %.2f\n",
                                                                        iatom, (p_master+iatom)->label,  
                                                                        jatom, (p_master+jatom)->label,
                                                                        ij_distance);

                             
 			  }
                     }
                }
            }

/**check if it's using centre, then copy a set of the original structure and proceed to      **/
/**modify the structure.							       	     **/
/**Processes :- 0 = have_centre ; default = normal permutation using fill_structure function **/
  	    if ( !throw )
	    {
             switch (Processes)
	     {
	     case 0: printf("seen count = %d\n", seen_count);
	             if(p_permute->debug) printf(">>Debug<< in CASE processes!\n");
	             /* keeping the generated state in the seen list  */
	     	     switch ( seen_count )
	     	     {
	      	     case 0: //seen_perm[seen_count].label_i=(p_master+iatom)->label;
              	     	     //seen_perm[seen_count].label_j=(p_master+jatom)->label;
                             
			     if (p_permute->debug) printf(">>Debug<< in CASE seen!\n");
              	      	     seen_perm[seen_count].index_i=iatom;
              	      	     seen_perm[seen_count].index_j=jatom;
	      	      	     seen_perm[seen_count].r1=idist_to_centre;
              	             seen_perm[seen_count].r2=jdist_to_centre;
              	             seen_perm[seen_count].r3=ij_distance;
	      	             seen_perm[seen_count].state=state;
	      	             seen_count++;
			     good++;
	     		     if (!p_permute->check)
	     		     {
			      if (p_permute->debug) printf(">>Debug<< copying atoms to be modified!\n");
			      num_new_atoms=-1;
	      		      for (iloop=0;iloop<*p_num_atoms-1;iloop++)
	      		      {
  			       if (p_permute->debug) printf(">>Debug<< %d\n", iloop);
			       if ( iloop >= centre )
			       {
  			        if (p_permute->debug) printf(">>Debug<< more than centre\n");
				temp_molecule[iloop]=*(p_master+(iloop+1));
			   	num_new_atoms++;
			       }
			       else
			       {
  			        if (p_permute->debug) printf(">>Debug<< less than centre\n");
			   	temp_molecule[iloop]=*(p_master+iloop);
			   	num_new_atoms++;
			       }
			      }
			      /* changing the labels of selected atoms */
              		      label_num=1;
	      		      printf("Centre %s removed!\n", (p_master+centre)->label);
			      printf("Structure now has %d atoms.\n", num_new_atoms); 
	      		      printf("Changing %s and %s to %s\n", temp_molecule[iatom-1].label,
							      temp_molecule[jatom-1].label,
								p_permute->subs_elem);
	      		      sprintf(temp_molecule[iatom-1].elem, "%s", p_permute->subs_elem);
	      		      sprintf(temp_molecule[jatom-1].elem, "%s", p_permute->subs_elem);
	      		      sprintf(temp_molecule[iatom-1].label,"%s%d", temp_molecule[iatom-1].elem, label_num);
	      		      sprintf(temp_molecule[jatom-1].label,"%s%d", temp_molecule[jatom-1].elem,label_num+1);
          
	      		      for (iloop=0; iloop<num_new_atoms; iloop++) 
	      		      {
	       		       printf("%s %10.6f %10.6f %10.6f\n", temp_molecule[iloop].label,
							         temp_molecule[iloop].x, temp_molecule[iloop].y,
							         temp_molecule[iloop].z);
      	      		      }
			     }
		             break;
	
	             default: 
			for(iloop=seen_count-1;iloop>=0;iloop--)
	      	     	{
	       	         seen_r=0;
//               	         if(p_permute->debug)
//   	       		 {
 			  printf(">>DEBUG<<\n");
	        	  printf("current : %.2f, seenr1 : %.2f\n", idist_to_centre, seen_perm[iloop].r1);
	        	  printf("current : %.2f, seenr2 : %.2f\n", jdist_to_centre, seen_perm[iloop].r2);
	        	  printf("current : %.2f, seenr3 : %.2f diff: %.2f abs: %.2f\n", ij_distance, seen_perm[iloop].r3,
                                                                                         ij_distance - seen_perm[iloop].r3,    
                                                                                    fabs(ij_distance - seen_perm[iloop].r3));
//	       		 }
	       		 if( fabs(ij_distance - seen_perm[iloop].r3) < 1.0E-3  )
	       		 {
                           printf("Failed on ij test\n");
			//compared i with r1 and j with r2, vice versa to eliminate similar permutations
	        	  if(     fabs(idist_to_centre - seen_perm[iloop].r1 ) < 1.0E-3
                              &&  fabs(jdist_to_centre - seen_perm[iloop].r2 ) < 1.0E-3)
			  {
			   seen_r=1; 
			   throw=TRUE; break;
			  }
	        	  else if (    fabs(idist_to_centre - seen_perm[iloop].r2) < 1.0E-3
                                    && fabs(jdist_to_centre - seen_perm[iloop].r1) < 1.0E-3)
			  {
			   seen_r=1; 
		           throw=TRUE; break;
			  }
			 }
		        } break;
	     	     }

	      	     switch ( seen_r )
	     	     {
	      	      case 1 :  
			if(p_permute->debug)
			{
		 	 printf(">>DEBUG<< seen_r %d\n", seen_r);
	         	 printf("current : %.2f, seenr1 : %.2f\n", idist_to_centre, seen_perm[iloop].r1);
	         	 printf("current : %.2f, seenr2 : %.2f\n", jdist_to_centre, seen_perm[iloop].r2);
	         	 printf("current : %.2f, seenr3 : %.2f\n", ij_distance, seen_perm[iloop].r3);
		        }
			printf("This state will be discarded, same as state %llu\n", seen_perm[iloop].state);
			break;
	       
	      	      case 0 : 
			printf("This is a new state, keeping in seen list and will be generated!\n");
			//temp_molecule[seen_count]=*(p_master+iatom);
	        	//seen_perm[seen_count].label_i=temp_molecule[seen_count].label;
			//temp_molecule[seen_count]=*(p_master+jatom);
                	//seen_perm[seen_count].label_j=temp_molecule[seen_count].label;
                	seen_perm[seen_count].index_i=iatom;
                	seen_perm[seen_count].index_j=jatom;
                	seen_perm[seen_count].r1=idist_to_centre;
                	seen_perm[seen_count].r2=jdist_to_centre;
                	seen_perm[seen_count].r3=ij_distance;
	        	seen_perm[seen_count].state=state;
	        	seen_count++;
			good++;
	     		if (!p_permute->check)
	     		{
			 num_new_atoms=-1;
	      		 for (iloop=0;iloop<*p_num_atoms-1;iloop++)
	      		 {
			  if ( iloop >= centre )
			  {
			   temp_molecule[iloop]=*(p_master+(iloop+1));
			   num_new_atoms++;
			  }
			  else
			  {
			   temp_molecule[iloop]=*(p_master+iloop);
			   num_new_atoms++;
			  }
			 }
			/* changing the labels of selected atoms */
              		 label_num=1;
	      		 printf("Centre %s removed!\n", (p_master+centre)->label);
			 printf("Structure now has %d atoms.\n", num_new_atoms); 
	      		 printf("Changing %s and %s to %s\n", temp_molecule[iatom-1].label,
							      temp_molecule[jatom-1].label, p_permute->subs_elem);
	      		 sprintf(temp_molecule[iatom-1].elem, "%s", p_permute->subs_elem);
	      		 sprintf(temp_molecule[jatom-1].elem, "%s", p_permute->subs_elem);
	      		 sprintf(temp_molecule[iatom-1].label,"%s%d", temp_molecule[iatom-1].elem, label_num);
	      		 sprintf(temp_molecule[jatom-1].label,"%s%d", temp_molecule[jatom-1].elem, label_num+1);
          
	      		 for (iloop=0; iloop<num_new_atoms; iloop++) 
	      		 {
	       		  printf("%s %10.6f %10.6f %10.6f\n", temp_molecule[iloop].label,
							      temp_molecule[iloop].x, temp_molecule[iloop].y,
							      temp_molecule[iloop].z);
      	      		 }
	     		}     
			break;
             	     }break;
	    
	    default : if( ij_distance >= p_permute->mind )
           	      {
	    	       printf("Good! This configuration will be generated.\n");
            	       good++;
            	       fill_structure(p_master, &temp_molecule[0], &perm_indices[0], &distrib[0], 
                        		num_indices, &num_new_atoms, *p_num_atoms);
                       if ( !p_permute->check )
		       {
		        printf("Filled %d atoms:\n", num_new_atoms);
             		for (iii=0; iii<=num_new_atoms; iii++) printf("%s %10.6f %10.6f %10.6f\n",
										     temp_molecule[iii].label,
                                                                                     temp_molecule[iii].x,
                                                                                     temp_molecule[iii].y,
                                                                                     temp_molecule[iii].z);
		       }
 	   	      }
		      else throw=TRUE;
		      break;
	     }
	    }
/** If check tag not found, then writes carfiles and POSCARs **/
/* Writes out the car files for different permutations named after their states id */
/* and moves them to their respective directories 				   */ 
	    if ( !p_permute->check && !throw )
            {
             sprintf(filename, "perm_%llu.car", state);
             if (p_permute->debug) printf("Doing OK name sorted >>%s<< 12...\n", filename);

             open_file(&fp, &filename[0], "w");
         
             sprintf(c_title_line, "Inter_vasp permutation, structure %llu", state);
/*** Need to add iconf to title **/
             write_car( fp, &header_line[0], &title_line[0], &c_title_line[0],
    	                &date_line[0], &temp_molecule[0], &mol_number[0], TRUE,
                        pbc, p_abc, num_new_atoms+1, 1.0, 
                        TRUE, &super[0], p_latt_vec, 
                        p_recip_latt_vec, &fix_flags[0], p_magmom, num_magmom);

                        if (p_permute->debug) printf("Doing OK 9..written car file..\n");
	     fclose(fp);

/*** Writes out the poscar files for each permutations ***/
            

              generate_neighbours( &temp_molecule[0], num_new_atoms,
              	                   &types[0], &num_types,
                                   pbc, p_recip_latt_vec,  p_latt_vec,
                                   &spec_charges[0]);

     	      if (p_permute->debug) 
                {
                  printf("Found:\n %d atoms \n ",num_new_atoms);

     	          printf("Back from generate_neighbours with %d types\n\n", num_types);

     	          for (iiii= 0; iiii <= num_types; iiii++)
                    {
                     printf("%d %s %d\n", iiii+1, &(types[iiii].atom_type[0]),
                                                                     types[iiii].num);
       	            }
                  printf("Molecule 1 has %d members\n",
					num_new_atoms);

                }

              for (iiii= 0; iiii <= num_new_atoms; iiii++)
                {
/**** If no group number make one up ****/
/**** Dave December 2005             ****/

                   if (strncmp(&(temp_molecule[iiii].group_no[0])," ",1) == 0) 
                     {
                       sprintf(&(temp_molecule[iiii].group_no[0]),"X"); 
                     }
       	           if (temp_molecule[iiii].group_no[0] == '\0')
	             {
	 	sprintf(&(temp_molecule[iiii].group_no[0]),"X"); 
 	       }
               printf("%s (elem= %s, g_no >>%s<<) with %d neighbours : ", 
                                temp_molecule[iiii].label, 
                                temp_molecule[iiii].elem, 
                                temp_molecule[iiii].group_no, 
                                temp_molecule[iiii].num_neigh); 

               for (ineigh=0; ineigh< temp_molecule[iiii].num_neigh; ineigh++)
               {
                neigh_index= temp_molecule[iiii].neighb[ineigh];
                printf("%s ",temp_molecule[neigh_index].label);
               }
              printf("\n");
             }
 /**** then write poscar file *******/
	     printf("start writing poscar\n"); 
	     sprintf(poscar_output, "POSCAR_%llu", state);
	             
             open_file(&fp_vasp_output, &poscar_output[0], "w");
	     printf("before sort\n"); 
	     sort_by_elem( &temp_molecule[0], num_new_atoms, &types[0], num_types); 	     
	     printf("after sort\n"); 
	     
	     printf("Converting the %s to %s \n", filename, poscar_output);
	    
             scale_factor=1.0; 
	     write_poscar( fp_vasp_output, &temp_molecule[0], &fract_coords[0],
                    	   &types[0], num_types, p_latt_vec, &scale_factor,
			   num_new_atoms+1, &title_line[0], &c_title_line[0], pbc, 
			   FALSE, &fix_flags[0]); 


             fclose(fp_vasp_output);
	     
	     sprintf(command, "mkdir perm_%llu", state);
             system(command);

             sprintf(command, "mv perm_%llu.car POSCAR_%llu perm_%llu", state, state, state);
             system(command);

             printf("--------------------------\n");
   	    }
        num_perms++;
        printf("Doing OK, num_perms with correct site count= %d number accepted so far = %d\n", num_perms, good);
        }
    }
    
  printf("Total Number Of permutations %d and Good Configurations %d :\n", num_perms, good);
  free(seen_perm);
  seen_perm = NULL;
  return good;
  }
