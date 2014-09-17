/* structres.h */
/* Mine first then those from orig (tim) after for reference */

typedef struct
{
 int steric;
 double non_bonded;			        /* non_bonding energy */
 double vdw_rep;                                /* van der Waals repulsive */
 double vdw_disp;                               /* van der Waals dispersive */
 double charges;				/* coulombic energy */
 int acceptance;
 double minimizer_init_total;			/* minimizer non_bonding energy */
 double minimizer_end_total;			/* minimizer total energy */
 double minimizer_init_nonbond;			/* minimizer non_bonding energy */
 double minimizer_end_nonbond;			/* minimizer total energy */
} energy;

/* rudimentary statistics structure */
 typedef struct
 {
  int tries;
  int accepted;
 }stats;

typedef struct
{
double stretch;
} internal_energy;

typedef struct
{
int code; 
double energy;
} structure_list;
 
/* structure for the atoms */
typedef struct
{
  char  label[7];
  double x;
  double y;
  double z;
  double vx;
  double vy;
  double vz;
  double ax;
  double ay;
  double az;
  char pot[4];
  char group[5];
  char group_no[9];
  char elem[3];
  double  part_chge; 
  int nb_list;                  /* index of non-bonding parameters */
  double vdw;			/* van der waals radius of the atom */
  int num_neigh; 		/* number of neighbours */
  int neighb[15]; 		/*initialise -1=not bonded */
  int num_images;               /* Number of symmetry related images added Nov 98 DJW */
  int image[10];                /* image indexes initially allow only 10 images */
  int mol;                      /* molecule index for this atom */
  int neighb_stretch_list[10];   /* index for intra stretch potential */
  double theta;                 /* angle for rotation interpolation */
  double electrostatic_pot;     /* The electrostatic potential at the atom position */
  double mass;                  /* For the atomic mass of the atom */
}atom;

typedef struct
{
  char name[80];
} char_list;

typedef struct
{
  char label[7];
} labels;

typedef struct
{
  char atom1[3];
  char atom2[3];
} bond;

typedef struct
{
  char atom_type[10];
  int  num;
} atom_number;

typedef struct
{
  int start;
  int end;
} links;

typedef struct
{
double matrix[9];
double translation[3];
} symm_ops;

typedef struct
{
int start;
int end;
int num;
} list_partition;

typedef struct
{
char name[3];
int num;
} types;

typedef struct
{
char   label[7];
int    is_core;
double part_chge;
}charge_list;

typedef struct
{
 double dx[MAXATOMS];
 double dy[MAXATOMS];
 double dz[MAXATOMS];
} e_vec;

typedef struct
{
double comp[3];
} simple_vec;

typedef struct
{
int atom1;
int atom2;
} bonds;

/**** structure for density of states data *****/

typedef struct
{
  double energy;
  double up_dos;
  double down_dos;
  double up_totdos;
  double down_totdos;
} dos;

typedef struct
{
int  fx;
int  fy;
int  fz;
} coord_flags;

typedef struct
{
int group1[MAXATOMS];
int group2[MAXATOMS];
int group_type1;
int group_type2;
int num_grp1;
int num_grp2;
} group_lists;
