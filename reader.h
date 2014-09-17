enum { PRIME_DIRECTIVE= 1, SECONDARY_DIRECTIVE, TERTIARY_DIRECTIVE, REMAINING_DIRECTIVE, SIESTA_DIRECTIVES, SIESTA_2ND_DIRECTIVES };

enum {TITLE = 1, 
      MASTER_FILE, END_FILE, GROUP_CENTRE, VARIABLE_SITES, AMPLITUDE, \
      LINEAR, NUM_STEPS, TEMPERATURE, \
      MIN_WEIGHT, ASSESS, NUMBER_TO_SET, \
      ANALYSE, NUM_PER_FORMULA_UNIT, START_AT_CODE, POTCAR_FILE, \
      OUTCAR_FILE, NEEDED, INTERPOLATE, MODE, MOL_CENTRE, LATE_CENTRE,\
      SWITCH, MORPH, SUPER, ANGLE_INTER, SHIFT, IMAGES, DOSCAR_FILE,\
      DOS_SMEAR, MDTRAJ, RESTART, MILLER, MD_RUN, COMP_MODES, MINIMAGE_END, INCAR_FILE,
      OUTCARS, EMIN, EMAX };

enum { CAR_FILE_NEEDED = 101, POSCAR_FILE_NEEDED, ARC_FILE_NEEDED,\
       FREQUENCY, FORCES, CONVERGENCE, PDB_FILE_NEEDED, GULP_FILE_NEEDED, \
       SHELLS, DOS_FILE_NEEDED, PART_DOS_NEEDED, ENERGY };

enum {MAYBE = 201, YES, NO, ALL};

enum {SI_TITLE = 301, BLOCK, NUMFREE};

enum {CHEMSPEC = 401, LATTICE, ATMCOORDS };

enum {CENTRE_TYPE=1, ANGLE_TYPE};

enum {DEFAULT = 666};    
enum {UNIT = 998};    
enum {PARSE = 1001};

enum {NODIRECT = 999};
enum {BLANK_DIRECT = 999};

#define PRIME_DIRECTIVE_LIST \
	"titl", TITLE, "mast", MASTER_FILE, "end_", END_FILE, \
        "grou",GROUP_CENTRE, \
        "vari", VARIABLE_SITES, "ampl", AMPLITUDE, \
        "line", LINEAR, "step", NUM_STEPS, \
        "temp", TEMPERATURE, "min_", MIN_WEIGHT, \
        "asse", ASSESS, "num_", NUMBER_TO_SET, \
        "anal", ANALYSE, "npfu", NUM_PER_FORMULA_UNIT, "star", START_AT_CODE, \
        "potc", POTCAR_FILE, "outc", OUTCAR_FILE, \
        "need", NEEDED, "inte", INTERPOLATE, "mode",MODE,"mol_",MOL_CENTRE,\
        "late",LATE_CENTRE,"swit",SWITCH,"morp",MORPH,"supe", SUPER, "angl", ANGLE_INTER,\
        "shift",SHIFT,"imag",IMAGES,"dosc", DOSCAR_FILE,\
        "smea", DOS_SMEAR,"traj",MDTRAJ, "rest", RESTART, \
        "mill", MILLER, "md_r", MD_RUN, "comp", COMP_MODES, \
        "endm", MINIMAGE_END, "inca", INCAR_FILE, "outc", OUTCARS,"emin", EMIN,\
        "emax", EMAX, "",NODIRECT

#define SECOND_DIRECTIVE_LIST \
              "car_", CAR_FILE_NEEDED, "arc_", ARC_FILE_NEEDED, \
              "posc", POSCAR_FILE_NEEDED,\
              "freq", FREQUENCY, "forc", FORCES, "conv", CONVERGENCE,\
              "pdb_",PDB_FILE_NEEDED,"gulp", GULP_FILE_NEEDED,\
              "shel", SHELLS, "dos_", DOS_FILE_NEEDED, "part", PART_DOS_NEEDED,\
              "ener", ENERGY, "",NODIRECT

#define THIRD_DIRECTIVE_LIST \
	   "", NODIRECT

#define NULL_DIRECTIVE_LIST \
        "", NODIRECT

#define SIESTA_DIRECTIVE_LIST \
	"systemname", SI_TITLE, "%blo", BLOCK,"numberoffreeatoms", NUMFREE, \
        "",NODIRECT

#define SIESTA_2ND_DIRECTIVE_LIST \
	"chemicalspecieslabel", CHEMSPEC, "latticeparameters", LATTICE, "atomiccoord", ATMCOORDS, \
        "",NODIRECT

#define END_OF_BLOCK "%endblock" /* End of SIESTA parameter block */

typedef struct 
{
  char *directive;
  int token_index;
}list;

static char  target[BUFFER];
static char  buf[BUFFER];
static char  *line, *last_tok;
extern int   read_new_line, line_no;

