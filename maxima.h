#define MAX_LINE_LEN 1024
#define MAX_LINES    1024
#define BUFFER  1024
#define PORE_GROUP pore
#define END_CAR "end"
#define NUM_POTENTIALS 28
#define NUM_ELEMENTS 104  /* update when more found we found deuterium! */

#define MAXATOMS 1000 /* used as max size of pore - will get round to malloc */
#define MAX_NEIGHS 15   /* used as max size of neighbours in the atom structure */
#define MAXMOL 500 /* Maximum number of discrete molecules */
#define MAXTYPES 40 /* used as max types of atom */
#define MAXMODES 300 /* used as max modes in frequency reading */
#define MAXIMAGES 20 /* used as max images in a reaction co-ordinate plot */
#define MAXOUTCARS 20 /* used as max OUTCAR files a reaction co-ordinate plot */
#define MAXBONDS 5  /* maximum number of bond changes */
#define MAX_MILLER  5 /* maximum number of Miller index sets */
#define MAX_MILLER3 15 /* maximum number of Miller indices, always 3xMAX_MILLER */
#define MAX_PERI_IMAGES 100 /* maximum number of periodic images of atoms */
#define MAX_PERI_IMAGES3  300 /* 3x maxmum number of periodic images of atoms for related vectors */
#define FILELEN_MAX 128

/****************************************************************/
/***** Maximum array dimension for users of ends in searches ****/
/****************************************************************/

#define MAX_ENDS 400

/****************************************************************/
/***** Maximum list length for pair interaction lists ***********/
/****************************************************************/

#define MAX_PAIR_LIST 1000

/****************************************************************/
/***** Maximum list for DOS reading *****************************/
/****************************************************************/

#define MAX_DOS 1000
#define MAX_DOS_FILES 700

/****************************************************************/
/** Maximum list length for matching function *******************/
/****************************************************************/

#define MAX_MATCH_LIST 10
