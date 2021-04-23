/**********************************************************************************/
/**** set_defaults puts in starting values for control variables ******************/
/**** and other parameters.                                      ******************/
/**** started Mar 19 Dave Willock *************************************************/
/**********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "maxima.h" 
#include "structures.h"
#include "reader.h" 
#include "global_values.h" 

void set_defaults( is_list *p_is, have_list *p_have, need_list *p_need)
 {
/*** have flags for saying what information the code has ***/
p_have->out      =FALSE;
p_have->incar    =FALSE;
p_have->potcar   =FALSE;
p_have->grp      =FALSE;
p_have->mol      =FALSE;
p_have->miller   =FALSE;
p_have->report   =FALSE;
p_have->band     =FALSE;
p_have->tot      =FALSE;
p_have->labels   =FALSE;
p_have->transfer =FALSE;
p_have->perm_centre    = FALSE;
p_have->perm_subs_elem = FALSE;

/*** is flags for saying what types of files have been loaded **/
p_is->gulp       = FALSE;
p_is->car        = FALSE; 
p_is->cif        = FALSE; 
p_is->vasp       = FALSE;
p_is->siesta     = FALSE;
p_is->onetep     = FALSE;
p_is->cart       = FALSE;
p_is->fract      = FALSE;
p_is->end_cart   = FALSE;
p_is->end_fract  = FALSE;
p_is->centre     = FALSE;
p_is->end_gulp   = FALSE;
p_is->end_car    = FALSE;
p_is->end_vasp   = FALSE;
p_is->siesta_dos = FALSE;
p_is->vasp_dos   = FALSE;
p_is->restricted = FALSE;

/*** need flags for saying what the user has requested **/
p_need->car        =FALSE;
p_need->cif        =FALSE;
p_need->poscar     =FALSE;
p_need->onetep     =FALSE;
p_need->arc        =FALSE;
p_need->freq       =FALSE;
p_need->force      =FALSE;
p_need->energy     =FALSE;
p_need->fermi      =FALSE;
p_need->pdb        =FALSE;
p_need->gulp       =FALSE;
p_need->shells     =FALSE;
p_need->expansion  =FALSE;
p_need->morph      =FALSE;
p_need->angle      =FALSE;
p_need->late       =FALSE;
p_need->shift      =FALSE;
p_need->dos        =FALSE;
p_need->part_dos   =FALSE;
p_need->mdtraj     =FALSE;
p_need->multi_dos  =FALSE;
p_need->md_run     =FALSE;
p_need->monit      =FALSE;
p_need->zsort      =FALSE;
p_need->miller_sort=FALSE;
p_need->permute    =FALSE;
p_need->interpolate=FALSE;
p_need->react_coord=FALSE;
p_need->poscar_frac=FALSE;
p_need->hbond      =FALSE;
p_need->match      =FALSE;
p_need->oshift     =FALSE;
p_need->vgap       =FALSE;

 return 0;
 }
