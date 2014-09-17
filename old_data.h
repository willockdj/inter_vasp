#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

typedef struct
{
  char atom[6];
  char pot[6];
  double a;
  double b;
  double sqrt_a;
  double sqrt_b;
} potential_list;

typedef struct
{
  char elem[3];
  int  number;
  double mass;
  double covalent;
  double vdw;
} periodic_table_data;    /* for vdw and other information */


typedef struct
{
  char label1[8];
  char label2[8];
  double bond_length;
}bond_dist;

typedef struct
{
  char label[3];
  double value;
}fractions;

typedef struct
{
  char type[10];
  char combination[10];
}potential_info;

/* more complex coord for Template with bits for valencies etc */
/* under development!*/
typedef struct
{
  char name[6];
  double x_coord;
  double y_coord;
  double z_coord;
  int    free_valence;
  int    include;

} complex_coord;

#ifdef MAIN

EXTERNAL int num_bond_list = 7;
EXTERNAL bond_dist bond_table[] = {
{"C", "C", 1.51},
{"C", "N", 1.46},
{"N", "N", 1.46},
{"N", "H", 1.01},
{"C","H", 1.1},
{"O","H", 1.0},
{"Si","O",1.70},
};

EXTERNAL int num_fractions = 5;
EXTERNAL fractions fraction_list[] = {
{"1/4", 0.25},
{"3/4", 0.75},
{"1/2", 0.5},
{"1/3", 1.0/3.0},
{"2/3", 2.0/3.0},
};

/*****non bonding parameters from cvff forcefield*****/
/*****as used by Docker - additional terms for   *****/
/*****zeolites                                   *****/
EXTERNAL potential_list potent[] = {
{  "H" ,  "h", 7108.47, 32.8708 , NULL, NULL }, 
{  "C" ,  "cg", 1.79034e+06, 528.482 , NULL, NULL },
{  "O" ,  "o", 272895.0, 498.879 , NULL, NULL },
{  "N" ,  "n", 2.26687e+06, 1230.56 , NULL, NULL },
{  "C" ,  "c'", 2.96875e+06, 1325.71 , NULL, NULL },
{  "C" ,  "c", 1.98105e+06, 1126.0 , NULL, NULL },
{  "H" ,  "hn", 1e-08, 0.0 , NULL, NULL },
{  "S" ,  "s", 365906.0, 250.8 , NULL, NULL },
{  "O" ,  "o*", 629358.0, 625.5 , NULL, NULL },
{  "H" ,  "h*", 1e-08, 0.0 , NULL, NULL },
{  "P" ,  "p", 6.02589e+06, 2195.6 , NULL, NULL },
{  "C" ,  "c+", 119025.0, 240.25 , NULL, NULL },
{  "Si" ,  "si", 3.14918e+06, 710.0 , NULL, NULL },
{  "F" ,  "f", 201106.0, 235.2 , NULL, NULL },
{  "Cl" ,  "cl", 1.05917e+06, 541.0 , NULL, NULL },
{  "Br" ,  "br", 3.57203e+06, 1195.0 , NULL, NULL },
{  "Na" ,  "Na", 14000.0, 300.0 , NULL, NULL },
{  "Cl" ,  "Cl", 6238.29, 51.6719 , NULL, NULL },
{  "O" ,  "o", 272895.0, 498.879 , NULL, NULL },
{  "C" ,  "c1", 1.98105e+06, 1126.0 , NULL, NULL },
{  "C" ,  "c2", 1.98105e+06, 1126.0 , NULL, NULL },
{  "C" ,  "c3", 1.98105e+06, 1126.0 , NULL, NULL },
{  "N" ,  "n4", 2.26687e+06, 1325.71 , NULL, NULL },
{  "Si" ,  "sz", 3.14918e+06, 710.0 , NULL, NULL },
{  "O" ,  "oss", 272895.0, 498.879 , NULL, NULL },
{  "C" ,  "cp", 2.96875e+06, 1325.71 , NULL, NULL },
{  "Ca" ,  "ca", 2.96875e+06, 1325.71 , NULL, NULL },
{  "Al" ,  "az", 3.14918e+06, 710.0 , NULL, NULL },
};

/*	Element, Z, mass, covalent radius, van der Waals radius*/

EXTERNAL periodic_table_data period_table[] = {
{ "Ac", 89,  227.028,  1.88,  2.81  },
{ "Ag", 47,  107.868,  1.59,  2.37  },
{ "Al", 13,  26.982,  1.35,  2.01  },
{ "Am", 95,  243.000,  1.51,  2.25  },
{ "Ar", 18,  39.948,  0.00,  0.00  },
{ "As", 33,  74.922,  1.21,  2.20  },
{ "At", 85,  210.000,  0.00,  0.00  },
{ "Au", 79,  196.967,  1.50,  1.85  },
{ "B", 5,  10.811,  0.83,  1.70  },
{ "Ba", 56,  137.327,  1.34,  2.00  },
{ "Be", 4,  9.012,  0.35,  0.52  },
{ "Bi", 83,  208.980,  1.54,  2.30  },
{ "Bk", 97,  247.000,  0.00,  0.00  },
{ "Br", 35,  79.904,  1.21,  1.95  },
{ "C", 6,  12.011,  0.68,  1.75  },
{ "Ca", 20,  40.078,  0.99,  1.48  },
{ "Cd", 48,  112.411,  1.69,  2.52  },
{ "Ce", 58,  140.115,  1.83,  2.73  },
{ "Cf", 98,  251.000,  0.00,  0.00  },
{ "Cl", 17,  35.453,  0.99,  1.77  },
{ "Cm", 96,  247.000,  0.00,  0.00  },
{ "Co", 27,  58.933,  1.33,  1.99  },
{ "Cr", 24,  51.996,  1.35,  2.01  },
{ "Cs", 55,  132.905,  1.67,  2.49  },
{ "Cu", 29,  63.546,  1.52,  1.54  },
{ "Dy", 66,  162.500,  1.75,  2.61  },
{ "Er", 68,  167.260,  1.73,  2.58  },
{ "Es", 99,  252.000,  0.00,  0.00  },
{ "Eu", 63,  151.965,  1.99,  2.97  },
{ "F", 9,  18.998,  0.64,  1.30  },
{ "Fe", 26,  55.847,  1.34,  2.00  },
{ "Fm", 100,  257.000,  0.00,  0.00  },
{ "Fr", 87,  223.000,  0.00,  0.00  },
{ "Ga", 31,  69.723,  1.22,  1.82  },
{ "Gd", 64,  157.250,  1.79,  2.67  },
{ "Ge", 32,  72.610,  1.17,  1.75  },
{ "H", 1,  1.008,  0.23,  1.17  },
{ "D", 1,  2.000,  0.23,  1.17  },
{ "He", 2,  4.003,  0.00,  0.00  },
{ "Hf", 72,  178.490,  1.57,  2.34  },
{ "Hg", 80,  200.590,  1.70,  1.90  },
{ "Ho", 67,  164.930,  1.74,  2.60  },
{ "I", 53,  126.905,  1.40,  2.10  },
{ "In", 49,  114.820,  1.63,  2.43  },
{ "Ir", 77,  192.220,  1.32,  1.97  },
{ "K", 19,  39.098,  1.33,  1.99  },
{ "Kr", 36,  83.800,  1.89,  2.82  },
{ "La", 57,  138.906,  1.87,  2.79  },
{ "Li", 3,  6.941,  0.68,  1.01  },
{ "Lr", 103,  260.000,  0.00,  0.00  },
{ "Lu", 71,  174.967,  1.72,  2.57  },
{ "Md", 101,  258.000,  0.00,  0.00  },
{ "Mg", 12,  24.305,  1.10,  1.64  },
{ "Mn", 25,  54.938,  1.35,  2.01  },
{ "Mo", 42,  95.940,  1.47,  2.19  },
{ "N", 7,  14.007,  0.68,  1.55  },
{ "Na", 11,  22.990,  0.97,  1.45  },
{ "Nb", 41,  92.906,  1.48,  2.21  },
{ "Nd", 60,  114.240,  1.81,  2.70  },
{ "Ne", 10,  20.180,  0.00,  0.00  },
{ "Ni", 28,  58.690,  1.50,  1.81  },
{ "No", 102,  259.000,  0.00,  0.00  },
{ "Np", 93,  237.048,  1.55,  2.31  },
{ "O", 8,  15.999,  0.68,  1.40  },
{ "Os", 76,  190.200,  1.37,  2.04  },
{ "P", 15,  30.974,  1.05,  1.90  },
{ "Pa", 91,  231.036,  1.61,  2.40  },
{ "Pb", 82,  207.200,  1.54,  2.30  },
{ "Pd", 46,  106.420,  1.50,  2.24  },
{ "Pm", 61,  145.000,  1.80,  2.69  },
{ "Po", 84,  209.000,  1.68,  2.51  },
{ "Pr", 59,  104.908,  1.82,  2.72  },
{ "Pt", 78,  195.080,  1.50,  1.97  },
{ "Pu", 94,  244.000,  1.53,  2.28  },
{ "Ra", 88,  226.025,  1.90,  2.84  },
{ "Rb", 37,  85.468,  1.47,  2.19  },
{ "Re", 75,  186.207,  1.35,  2.01  },
{ "Rh", 45,  102.906,  1.45,  2.16  },
{ "Rn", 86,  222.000,  0.00,  0.00  },
{ "Ru", 44,  101.070,  1.40,  2.09  },
{ "S", 16,  32.066,  1.02,  1.80  },
{ "Sb", 51,  121.750,  1.46,  2.17  },
{ "Sc", 21,  44.956,  1.44,  2.15  },
{ "Se", 34,  78.960,  1.22,  2.00  },
{ "Si", 14,  28.086,  1.20,  2.00  },
{ "Sm", 62,  150.360,  1.80,  2.69  },
{ "Sn", 50,  118.710,  1.46,  2.18  },
{ "Sr", 38,  87.620,  1.12,  1.67  },
{ "Ta", 73,  180.948,  1.43,  2.13  },
{ "Tb", 65,  158.925,  1.76,  2.63  },
{ "Tc", 43,  98.000,  1.35,  2.01  },
{ "Te", 52,  127.600,  1.47,  2.20  },
{ "Th", 90,  232.038,  1.79,  2.67  },
{ "Ti", 22,  47.880,  1.47,  2.19  },
{ "Tl", 81,  204.383,  1.55,  2.31  },
{ "Tm", 69,  168.934,  1.72,  2.57  },
{ "U", 92,  238.029,  1.58,  2.36  },
{ "V", 23,  50.942,  1.33,  1.99  },
{ "W", 74,  183.850,  1.37,  2.04  },
{ "Xe", 54,  131.290,  0.00,  0.00  },
{ "Y", 39,  88.906,  1.78,  2.66  },
{ "Yb", 70,  173.040,  1.94,  2.90  },
{ "Zn", 30,  65.390,  1.45,  2.16  },
{ "Zr", 40,  91.224,  1.56,  2.33  },
};

#else
EXTERNAL int num_bond_list;
EXTERNAL bond_dist bond_table[] ;
EXTERNAL int num_fractions;
EXTERNAL fractions fraction_list[];
EXTERNAL potential_list potent[] ;
EXTERNAL periodic_table_data period_table[];
EXTERNAL potential_info pot_info;

#endif


/******** Some useful #defines to clarify code ********/
/******** Dave Willock started 24th May 1995   ********/

/****** Cost_function *********************************/

#define STERIC_COST 1

/****** Main Code: actions ****************************/
#define NUMBER_OF_ACTIONS 7

#define BUILD_ACTION              0
#define ROTATE_LAST_FRAG_ACTION   1
#define SHAKE_ACTION              2
#define ROCK_ACTION               3
#define TWIST_ACTION              4
#define MINIMIZE_ACTION           5
#define MINIMIZE_INPORE_ACTION    6

/****** Minimizer names *********************************/
#define DISCOVER_MINIMIZER "discover"
#define MOPAC_MINIMIZER "mopac"
#define INTERNAL_MINIMIZER "internal"
#define MINIMIZER_DEFAULT "discover"

#define DEFAULT_MOPAC_PATH "/progs/gulp/mopac"
#define DEFAULT_DISCOVER_PATH "/usr/biosym/237/biosym_bin/discover"

#undef EXTERNAL
