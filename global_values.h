#define BOOLEAN int
#define TRUE 1
#define FALSE 0

#define LINESIZ 200

/********************************************************/
/** Values to tell us what type of GULP file we have ****/
/** Expect read_gulp to always return +ve numbers if ****/
/** a good read occured                              ****/
/********************************************************/

#define GULP_CART    2    
#define GULP_FRACT   3

#define END_OF_INPUT "Thatsit" /* As defined in read_line */

/*** BOND_TOL used by generate_neighbours to set upper bound of bond distances ***/
/***          this is a fractional increment used with standard bond values    ***/
/*** CLASH_TOL used to define when atoms are too close to allow a POSCAR to be ***/
/*** written. This is an absolute value in Angstroms                           ***/
#define BOND_TOL 0.1
#define CLASH_TOL 0.3

