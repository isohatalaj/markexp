
#include "mexp.h"

/**
 * Shorthands for the normal distribution pdf, cdf, and inverse cdf.
 */
inline double  gpdf(double x) { return gsl_ran_ugaussian_pdf(x);  }
inline double  gcdf(double x) { return gsl_cdf_ugaussian_P(x);    }
inline double gicdf(double x) { return gsl_cdf_ugaussian_Pinv(x); }

void 
errhandler(const char * reason, 
	   const char * file, 
	   int line, 
	   int gerrno)
{
  fprintf(stderr, "gsl error @ %s(%d): %s [%d]\n", 
	  file, line, reason, gerrno);
}

mexp_work_t * 
mexp_work_alloc()
{
  mexp_work_t *self;

  gsl_set_error_handler(&errhandler);

  self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  self->intwork = NULL;

  self->tol_rel = 1e-9;
  self->tol_abs = 1e-9;

  self->intwork_size = 512;
  self->intwork = gsl_integration_workspace_alloc(self->intwork_size);
  if (self->intwork == NULL) goto fail;

  return self;
  
 fail:
  mexp_work_free(self);

  return NULL;
}

void
mexp_work_free(mexp_work_t *self)
{
  if (self == NULL) return;
  if (self->intwork != NULL) gsl_integration_workspace_free(self->intwork);

  free(self);
}

/* TODO: This standard form of the integral can be simplified
 * somewhat. For now, keep it simple, stupid. */
static double
stdinteg_fun(double x, void *params)
{
  double *v = params;
  return exp(v[0]*x)*gcdf(v[1]*x+v[2])*gpdf(x);
}

int
mexp_work_stdinteg(int limtype, double alpha, double omega,
		   double a, double b, double c,
		   double *result,
		   mexp_work_t *work)
{
  gsl_function f;
  double v[3];
  v[0] = a;
  v[1] = b;
  v[2] = c;
  f.function = stdinteg_fun;
  f.params = v;

  double abserr;
  int status;

  switch (limtype)
    {
    case MEXP_FIN_FIN:
      status = gsl_integration_qag(&f, alpha, omega,
				   work->tol_abs, work->tol_rel,
				   work->intwork_size,
				   GSL_INTEG_GAUSS31, /* <- may not be best. */
				   work->intwork, result, &abserr);
      break;

    case MEXP_FIN_INF:
      status = gsl_integration_qagiu(&f, alpha,
				     work->tol_abs, work->tol_rel,
				     work->intwork_size,
				     work->intwork,
				     result, &abserr);
      break;

    case MEXP_INF_FIN:
      status = gsl_integration_qagil(&f, omega,
				     work->tol_abs, work->tol_rel,
				     work->intwork_size,
				     work->intwork,
				     result, &abserr);
      break;

    case MEXP_INF_INF:
      status = gsl_integration_qagi(&f,
				    work->tol_abs, work->tol_rel,
				    work->intwork_size,
				    work->intwork,
				    result, &abserr);
      break;
      
    default:
      GSL_ERROR("Invalid integration range spec", GSL_EINVAL);
    }

  if (status) GSL_ERROR("Numerical integration failed", GSL_EFAILED);

  /* TODO: Check that abserr is not too large compared to the value of
   * the integral. */

  return GSL_SUCCESS;
}


int
mexp_map_coefs(double eps_dag_1, double eps_dag_dag_1,
	       double eps_dag_2,
	       double *Iota1, double *Kappa1, double *Lambda1,
	       double *Lambda2,
	       mexp_pars_t *p, mexp_work_t *w)
{
  *Iota1   = exp(p->phi1*(0.5*p->phi1-eps_dag_1))*gcdf(eps_dag_dag_1-p->phi1); //ok
  *Kappa1  = -gcdf(eps_dag_1) + exp(p->phi1*(0.5*p->phi1-eps_dag_1))*gcdf(eps_dag_1-p->phi1); //ok
  *Lambda1 = 1.0 - gcdf(eps_dag_1) - exp(p->phi1*(0.5*p->phi1-eps_dag_1))*(gcdf(eps_dag_dag_1 - p->phi1) - gcdf(eps_dag_1 - p->phi1)); //ok

  double xi_prime = sqrt(1.0 - p->xi*p->xi);
  double LI1, LI2, LI3, LI4;
  int status; 

  status = mexp_work_stdinteg(MEXP_FIN_FIN, eps_dag_dag_1, eps_dag_1,
			      p->phi2*(p->varphi+p->xi),
			      -(p->varphi + p->xi)/xi_prime,
			      (eps_dag_2 + p->varphi*eps_dag_1)/xi_prime - p->phi2*xi_prime,
			      &LI1, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration I returned bad status",
			GSL_EFAILED);

  status = mexp_work_stdinteg(MEXP_FIN_FIN, eps_dag_dag_1, eps_dag_1,
			      0.0,
			      -(p->varphi+p->xi)/xi_prime,
			      (eps_dag_2 + p->varphi*eps_dag_1)/xi_prime,
			      &LI2, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration II returned bad status",
			GSL_EFAILED);
			      
  status = mexp_work_stdinteg(MEXP_FIN_INF, eps_dag_1, 0.0,
			      p->phi2*p->xi,
			      -p->xi/xi_prime,
			      eps_dag_2/xi_prime - p->phi2*xi_prime,
			      &LI3, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration III returned bad status",
			GSL_EFAILED);
			      
  status = mexp_work_stdinteg(MEXP_FIN_INF, eps_dag_1, 0.0,
			      0.0,
			      -p->xi/xi_prime,
			      eps_dag_2/xi_prime,
			      &LI4, w);
  if (status) GSL_ERROR("Failed computing map coefs as integration IV returned bad status",
			GSL_EFAILED);
			      

  *Lambda2 = 
    1.0 - gcdf(eps_dag_dag_1)
    + exp(0.5*p->phi2*p->phi2*(1-p->xi*p->xi)-p->phi2*(eps_dag_2+p->varphi*eps_dag_1))*LI1
    - LI2 
    + exp(0.5*p->phi2*p->phi2*(1-p->xi*p->xi)-p->phi2*eps_dag_2)*LI3
    - LI4
    ; //ok

  //ok'd 18.8. 10:54

  return GSL_SUCCESS;
}

static double
eps_star(double k, double X, double rho)
{
  return (k - rho*X)/sqrt(1 - rho*rho);
}

static double
X2_star(double Z2, const mexp_pars_t *p)
{
  return (1.0 - p->mu)*p->X1 + p->zeta*Z2;
}

int
mexp_single_objective(double k_dag_1,
		      double k_dag_dag_1,
		      double k_dag_2,
		      double *Omega,
		      mexp_pars_t *p, mexp_work_t *w)
{
  int status;

  double Z2_median = 0.0;
  double Z2_q0 = gicdf(p->q0);
  double Z2_q1 = gicdf(p->q1);

  double X2_median = X2_star(Z2_median, p);
  double X2_q0 = X2_star(Z2_q0, p);
  double X2_q1 = X2_star(Z2_q1, p);

  double eps_dag_1 = eps_star(k_dag_1, p->X1, p->rho);
  double eps_dag_dag_1 = eps_star(k_dag_dag_1, p->X1, p->rho);
  double eps_dag_2_median = eps_star(k_dag_2, X2_median, p->rho);
  double eps_dag_2_q0 = eps_star(k_dag_2, X2_q0, p->rho);
  double eps_dag_2_q1 = eps_star(k_dag_2, X2_q1, p->rho);
    
  double Iota_1, Kappa_1, Lambda_1, Lambda_2;
  double Delta_median, Delta_q0, Delta_q1, Delta_spread;

  /* Median */
  status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_median,
			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
			  p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);
  Delta_median = Lambda_2 - Lambda_1;
  

  /* Quantile q0 */
  status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_q0,
			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
			  p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);
  Delta_q0 = Lambda_2 - Lambda_1;

  /* Quantile q1 */
  status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_q1,
			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
			  p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);
  Delta_q1 = Lambda_2 - Lambda_1;


  /* Set output and we're done */
  Delta_spread = fabs(Delta_q1 - Delta_q0);

  *Omega = (1.0 - p->gamma)*Delta_median - p->gamma*Delta_spread;

  return GSL_SUCCESS;
}

int
mexp_k_dags_impacted(double k_dag_dag_1_A,
		     double k_dag_dag_1_B,
		     double *k_dag_1, double *k_dag_2,
		     mexp_pars_t *p)
{
  double eps_dag_dag_1_A = eps_star(k_dag_dag_1_A, p->X1, p->rho);
  double eps_dag_dag_1_B = eps_star(k_dag_dag_1_B, p->X1, p->rho);

  double weight_A = p->L0_A/(p->L0_A + p->L0_B);
  double weight_B = p->L0_B/(p->L0_A + p->L0_B);

  *k_dag_1 = p->k0_1 + p->psi*(weight_A*gcdf(eps_dag_dag_1_A) +
			       weight_B*gcdf(eps_dag_dag_1_B));

  *k_dag_2 = (1.0-p->chi)*p->k0_2 + p->chi*(*k_dag_1);

  return GSL_SUCCESS;
}

int
mexp_two_objectives(double k_dag_dag_1_A,
		    double k_dag_dag_1_B,
		    double *k_dag_1, double *k_dag_2,
		    double *Omega_A, double *Omega_B,
		    mexp_pars_t *p, mexp_work_t *w)
{
  int status;

  status = mexp_k_dags_impacted(k_dag_dag_1_A, k_dag_dag_1_B,
				k_dag_1, k_dag_2,
				p);
  if (status) GSL_ERROR("Failed computing price impact",
			GSL_EFAILED);

  status = mexp_single_objective(*k_dag_1, k_dag_dag_1_A, *k_dag_2,
				 Omega_A, p, w);
  if (status) GSL_ERROR("Failed computing objective function A",
			GSL_EFAILED);

  status = mexp_single_objective(*k_dag_1, k_dag_dag_1_B, *k_dag_2,
				 Omega_B, p, w);
  if (status) GSL_ERROR("Failed computing objective function B",
			GSL_EFAILED);

  return GSL_SUCCESS;				
}

