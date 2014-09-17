/**************************************************************/
/* Maths function normalisation and constants set at run time */
/* Dave Willock 7th May 1996                                  */
/**************************************************************/

#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

EXTERNAL double one_sixth;
EXTERNAL double one_third;
EXTERNAL double pi;
EXTERNAL double four_pi;
EXTERNAL double pi_tothehalf;

/**************************************************************/
/**** Error function normalisation will be 2/root pi **********/
/**************************************************************/

EXTERNAL double erf_normalise;

#undef EXTERNAL
