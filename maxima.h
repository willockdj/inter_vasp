#define MAX_LINE_LEN 1024
#define BUFFER  1024
#define PORE_GROUP pore
#define END_CAR "end"
#define NUM_POTENTIALS 28
#define NUM_ELEMENTS 104  /* update when more found we found deuterium! */

#define MAXATOMS 1000 /* used as max size of pore - will get round to malloc */
#define MAXTYPES 40 /* used as max types of atom */
#define MAXMODES 200 /* used as max modes in frequency reading */
#define MAXIMAGES 60 /* used as max images in a reaction co-ordinate plot */
#define MAXOUTCARS 60 /* used as max OUTCAR files a reaction co-ordinate plot */
#define MAXBONDS 5  /* maximum number of bond changes */
#define FILELEN_MAX 128

/****************************************************************/
/***** Maximum array dimension for users of ends in searches ****/
/****************************************************************/

#define MAX_ENDS 400

/****************************************************************/
/***** Maximum list length for pair interaction lists ***********/
/****************************************************************/

#define MAX_PAIR_LIST 100

/****************************************************************/
/***** Maximum list for DOS reading *****************************/
/****************************************************************/

#define MAX_DOS 700
#define MAX_DOS_FILES 700
