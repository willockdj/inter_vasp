/************************************/
/* THIS IS MY STUFF                 */
/************************************/

#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

/********************** read variables********************/

EXTERNAL int read_new_line;
EXTERNAL int line_no;


/********************** files pointers********************/

EXTERNAL FILE *input_fp; /*input file */
EXTERNAL FILE *output_fp;

/************ files and temporary file read variables ***************/

EXTERNAL char inputfile[FILELEN_MAX], outputfile[FILELEN_MAX];

/****************************************************************************/
/********** Variables used to keep track of potentials **********************/
/****************************************************************************/

EXTERNAL char title[80];
EXTERNAL char buffer[256];
EXTERNAL char dummy_head[256];   /* temp storage for biosym headers */


/********************************************************************/
/***** pbc    : flags that periodic boundary conditions are *********/
/*****          to be used                                  *********/
/***** abc    : holds a b c alpha beta gamma Angstroms      *********/
/*****                                     and degrees      *********/
/***** latt_vec : holds cartessian vectors for a b c        *********/
/***** recip_latt_vec : holds recip. space a* b* c*         *********/
/********************************************************************/

EXTERNAL int pbc;
EXTERNAL double abc[6];
EXTERNAL double latt_vec[9];
EXTERNAL double real_latt_sizes[3];
EXTERNAL double recip_latt_vec[9];
EXTERNAL double recip_latt_sizes[3];
EXTERNAL double cell_volume;

EXTERNAL int verbose; /* verbose output flag */

#undef EXTERNAL
