
#include <math.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

/* WHEN TOUCHING THIS, MAKE THE CORRESPONDING CHANGE TO THE PYTHON
 * INTERFACE!
 */
typedef struct {
  double phi1;   /* Unit conversion parameter, period 1 */
  double phi2;   /* Unit conversion parameter, period 2 */
  double varphi; /* */

  double k0_1;   /* Baseline impairment k, period 1 */
  double k0_2;   /* Baseline impairment k, period 2 */

  double c0_A;   /* Bank A initial capitalization ratio */
  double c0_B;   /* Bank B initial capitalization ratio */

  double L0_A;   /* Bank A initial loan pfolio value */
  double L0_B;   /* Bank B initial loan pfolio value */

  double rho;    /* Ideosync. vs aggregate */
  double xi;     /* Ideosync. time correlation */
  double mu;     /* Aggreg. factor drift. mu = 0 => X random walk, mu = 1 => X random normal every period */
  double zeta;   /* Aggreg. factor noise */

  double psi;    /* Price impact parameter */
  double chi;    /* Price impact persistence, k_imp_2 = (1-chi)k_0 + chi k_imp_1 */

  double gamma0; /* Period 2 capital weight in objective function */
  double gamma1; /* Profit spread weight in objective function */
  double q0;     /* Profit spread lower quantile */
  double q1;     /* Profit spread upper quantile */

  double gamma2; /* Tail risk weight in objective function */
  double c_bar;  /* Tail risk capitalisation ratio thold. */

  double X1;     /* Period one aggregate shock parameter */

  double dummy;  /* Temp variable to sneak an extra parameter into
		  * some functions. Yucky hack, I know. */
} mexp_pars_t;

typedef struct {
  gsl_root_fsolver *solver;
  int max_iter;

  gsl_integration_workspace *intwork;
  int intwork_size;

  double tol_rel;
  double tol_abs;
} mexp_work_t;

mexp_work_t * 
mexp_work_alloc();

void
mexp_work_free(mexp_work_t *work);

#define MEXP_FIN_FIN 0
#define MEXP_FIN_INF 1
#define MEXP_INF_FIN 2
#define MEXP_INF_INF 3

/* Integrate exp(ax)N(bx+c)f(x) over [alpha,omega], where f and N are
 * standard Gaussian pdf and cdf. */
int
mexp_work_stdinteg(int limtype, double alpha, double omega,
		   double a, double b, double c,
		   double *result,
		   mexp_work_t *work);

int
mexp_map_coefs(double eps_dag_1, double eps_dag_dag_1,
	       double eps_dag_2,
	       double *Iota1, double *Kappa1, double *Lambda1,
	       double *Lambda2,
	       mexp_pars_t *p, mexp_work_t *w);

int
mexp_single_objective(double k_dag_1,
		      double k_dag_dag_1,
		      double k_dag_2,
		      double c_0,
		      double L_0,
		      double *Omega,
		      mexp_pars_t *p, mexp_work_t *w);

int
mexp_k_dags_impacted(double k_dag_dag_1_A,
		     double k_dag_dag_1_B,
		     double *k_dag_1, double *k_dag_2,
		     mexp_pars_t *p);

int
mexp_two_objectives(double k_dag_dag_1_A,
		    double k_dag_dag_1_B,
		    double *k_dag_1, double *k_dag_2,
		    double *Omega_A, double *Omega_B,
		    mexp_pars_t *p, mexp_work_t *w);

/* Find the upper bound for the impairment threshold k (obviously then
 * also the upper bound for foreclosures). This is found by setting
 * the foreclosure thresholds of both banks to the impairment
 * threshold in the expression for the fire-sale impact. */
int
mexp_find_k_fcl_ubound(double *k_fcl_ubound,
		       mexp_pars_t *p, mexp_work_t *w);

int
mexp_find_k_fcl_A_ubound(double *k_fcl_A_ubound,
			 double k_fcl_B,
			 double k_fcl_ubound,
			 mexp_pars_t *p,
			 mexp_work_t *w);

int
mexp_find_k_fcl_B_ubound(double *k_fcl_B_ubound,
			 double k_fcl_A,
			 double k_fcl_ubound,
			 mexp_pars_t *p,
			 mexp_work_t *w);

double
mexp_eps_star(double k, const mexp_pars_t *p);

double
mexp_N_star(double k, const mexp_pars_t *p);

/* Convert tilde fractions (ratio of foreclosures to max. number of
   impaired loans across all strategies of all banks) to values of k. q*/
double
mexp_k_of_p_tilde(double p_tilde, double k_bound,
		  mexp_pars_t *p);


