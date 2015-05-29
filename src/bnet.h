
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

/**
 * Structure for holding model parameters.
 *
 */
/* WHEN MODYIFING THIS STRUCT, REMEMBER TO UPDATE NOT ONLY THE C CODE,
 * SETTER FUNCTION, ETC, BUT THE **PYTHON INTERFACE** AS WELL -- THIS
 * IS IMPORTANT AS THINGS CAN GO SILENTLY WRONG IF THIS IS NOT DONE!
 */
typedef struct {
  double LGD;    /*< Loss given default. */
  double PD;     /*< Unconditional default rate. */

  double rho;    /*< Aggregate/ideosync. shock correlation. */
  double psi;    /*< Price impact coefficient. */

  double gamma;  /*< Risk weight in bank objective. */
  double q0;     /*< Lower bound quantile for objective risk measure. */
  double q1;     /*< Upper bound quantile for objective risk measure. */

  double phi;    /*< Single default revaluation exponent. */
  double phi2;   /*< Double default revaluation exponent. */

  double chi;    /*< Weighing for setting period 2 loan impairment k
		     threshold. If k0 is the baseline, kimp is period
		     1 impairment, then period 2 impairment kimp2 =
		     (1-chi)*k0 + chi*k1. */

} bnet_pars_t;


/**
 * Work structure for the solver routines. This holds a few gsl
 * objects (root solvers and integration workspace) needed by our
 * routines, and also some parameters relevant to the numerical
 * solution of the model equations. These objects are allocated and
 * destroyed using bnet_work_alloc and bnet_work_free.
 */
typedef struct {
  gsl_root_fsolver *solver;           /*< First root solver. */
  gsl_root_fsolver *solver_2;         /*< Second solver. */
  /* We really need two solvers: First used in solving profit-loss
     distribution, second needed for finding quantiles of that
     distribution. If the quantile solver used the same object, its
     state would get corrupted everytime the distribution function was
     evaluated. */

  gsl_integration_workspace *intwork; /*< Integration workspace. For
				          all our numerical
				          integration needs. */
  double tol_rel;   /*< Relative tolerance used in root finding. */   
  double tol_abs;   /*< Absolute tolerance used in root finding. */
  /* For the specific meaning of these tolerances, see GSL manual,
   * root fsolvers, gsl_root_test_interval function. */
  int max_iter;     /*< Max. number of root finding iterations. */

  int intwork_size; /*< Size of the integration work array. */

} bnet_work_t;

bnet_work_t *
bnet_work_alloc();

void
bnet_work_free(bnet_work_t *);

bnet_pars_t *
bnet_pars_alloc();

void
bnet_pars_free(bnet_pars_t *self);

double
bnet_comp_kimp(const double epsfc[2], 
	       const double L[2], double X, 
	       double k0, double psi);

int
bnet_calibrate_phi(bnet_work_t *work,
		   bnet_pars_t *pars,
		   double *phi_calib);

int
bnet_comp_Omega(double c0,
		double X,
		double kfc, double kimp,
		double *Omega, 
		double *c1_out, 
		double *Ec2_out,
		double *ipr_out,
		bnet_pars_t *p,
		bnet_work_t *w);

int
bnet_comp_Omegas(double C[2],
		 double L[2],
		 double X, 
		 double kfc[2],
		 double Omega[2],
		 bnet_pars_t *p,
		 bnet_work_t *w);

double
bnet_comp_kimp(const double epsfc[2], const double L[2], double X, 
	       double k0, double psi);

double
bnet_eps_star(double k, double X, double rho);

double
bnet_kzero(double PD);

