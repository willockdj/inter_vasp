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
  int neighb[MAX_NEIGHS]; 		/*initialise -1=not bonded */
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
char cnt_lab1[7];
char cnt_lab2[7];
int centre[2];
char axis1_lab[7];
char axis2_lab[7];
int axis1_ind[2];
int axis2_ind[2];
} group_lists;

/**** structure for controlling permutation of atoms *****/

typedef struct
{
char elem[10];
char subs_elem[10];
double mind;
double maxd;
int num;
int check;
int centre;
int debug;
} perms;

typedef struct
{
char label_i[7];
char label_j[7];
int index_i;
int index_j;
unsigned long long int state;
double r1;
double r2;
double r3;
} seen_lists;

/** Structures for holding sets of control variables ***/

typedef struct
{
int out;
int incar;
int grp;
int mol;
int miller;
int report;
int band;
int tot;
int labels;
int transfer;
int perm_centre;
int perm_subs_elem;
int potcar;
} have_list;

typedef struct
{
int gulp;
int car;
int cif;
int pdb;
int vasp;
int siesta;
int onetep;
int punch;
int cart;
int fract;
int end_cart;
int end_fract;
int centre;
int end_gulp;
int end_car;
int end_vasp;
int siesta_dos;
int vasp_dos;
int restricted;
} is_list;

typedef struct
{
int car;
int poscar;
int onetep;
int cif;
int arc;
int freq;
int force;
int energy;
int fermi;
int pdb;
int gulp;
int shells;
int expansion;
int morph;
int angle;
int late;
int shift;
int dos;
int part_dos;
int mdtraj;
int multi_dos;
int md_run;
int monit;
int zsort;
int oshift;
int vgap;  
int miller_sort;
int permute;
int interpolate;
int react_coord;
int poscar_frac;
int potcar;
int hbond;
int match;
} need_list;

