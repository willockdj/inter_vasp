/**************************************************************/
/* Ewald sum constants and accuracy parameter                 */
/* Dave Willock 13th June 1996                                */
/**************************************************************/

#ifdef MAIN
#define EXTERNAL
#else
#define EXTERNAL extern
#endif

EXTERNAL double kappa;
EXTERNAL double kappa_sqrd;
EXTERNAL double two_pi;
EXTERNAL double four_pi_sqrd;
EXTERNAL double four_pi_over_vol;
EXTERNAL double four_kappa_sqrd;
EXTERNAL double pi_sqrd_over_kappa_sqrd;
EXTERNAL double erf_accuracy;
EXTERNAL double ewald_accuracy;
EXTERNAL double f_param;
EXTERNAL double coul_prefactor;
EXTERNAL double recip_sum_max;
EXTERNAL double real_sum_max;
EXTERNAL double recip_sum_max2;
EXTERNAL double real_sum_max2;

#undef EXTERNAL
