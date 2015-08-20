
#include "mexp.h"

/**
 * Shorthands for the normal distribution pdf, cdf, and inverse cdf.
 */
inline double  gpdf(double x) { return gsl_ran_ugaussian_pdf(x);  }
inline double  gcdf(double x) { return gsl_cdf_ugaussian_P(x);    }
inline double gicdf(double x) { return gsl_cdf_ugaussian_Pinv(x); }

void 
errhandler(const char *reason, 
	   const char *file, 
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
  self->solver = NULL;

  self->tol_rel = 1e-6;
  self->tol_abs = 1e-6;
  self->max_iter = 64;

  self->intwork_size = 512;
  self->intwork = gsl_integration_workspace_alloc(self->intwork_size);
  if (self->intwork == NULL) goto fail;

  self->solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  if (self->solver == NULL) goto fail;

  return self;
  
 fail:
  mexp_work_free(self);

  return NULL;
}

void
mexp_work_free(mexp_work_t *self)
{
  if (self == NULL) return;
  if (self->solver != NULL) gsl_root_fsolver_free(self->solver);
  if (self->intwork != NULL) gsl_integration_workspace_free(self->intwork);

  free(self);
}

static int
findroot(mexp_work_t *self, 
	 double (*fun)(double, void*), 
	 void *params, 
	 double x0, double x1,
	 double *root)
{
  int status;

  gsl_function f;
  f.function = fun;
  f.params   = params;

  double y0, y1;

  y0 = fun(x0, params);
  y1 = fun(x1, params);

  if (y0 == 0.0)
    {
      *root = x0;
      return GSL_SUCCESS;
    }

  if (y1 == 0.0)
    {
      *root = x1;
      return GSL_SUCCESS;
    }

  status = gsl_root_fsolver_set(self->solver, &f, x0, x1);
  if (status)
    {
      fprintf(stderr, 
	      "# error debug: "
	      "x0 = %lg, y0 = %lg; x1 = %lg, y1 = %lg\n",
	      x0, y0, x1, y1);
      GSL_ERROR("root solver could not be initialized", status);
    }

  double x;
  int iter = 0;
  do
    {
      iter++;
      status = gsl_root_fsolver_iterate(self->solver);
      if (status) GSL_ERROR("root solver iteration failed", status);
      x = gsl_root_fsolver_root(self->solver);
      status = gsl_root_test_interval(gsl_root_fsolver_x_lower(self->solver), 
				      gsl_root_fsolver_x_upper(self->solver),
				      self->tol_abs, 
				      self->tol_rel);
    }
  while (status == GSL_CONTINUE && iter < self->max_iter);

  if (status) GSL_ERROR("root solver iterations failed to converge", 
			status);

  *root = x;

  return GSL_SUCCESS;
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
  double s = 1.0;

  const int stdize = 1;
  if (stdize)
    {
      v[0] = 0.0;
      v[1] = b;
      v[2] = b*a + c;
      alpha -= a;
      omega -= a;
      s = exp(0.5*a*a);
    }
  else
    {
      v[0] = a;
      v[1] = b;
      v[2] = c;
    }

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

  *result *= s;

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

typedef struct {
  mexp_work_t *w;
  mexp_pars_t *p;
  double c_0;
  double L_0;
  double eps_dag_1;
  double eps_dag_dag_1;
} tailrisk_pars_t;

static double
tailrisk_fun(double eps_dag_2, void *pars)
{
  int status;
  tailrisk_pars_t *p = pars;
  double Iota_1, Kappa_1, Lambda_1, Lambda_2;

  status = mexp_map_coefs(p->eps_dag_1, p->eps_dag_dag_1, eps_dag_2,
  			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
  			  p->p, p->w);
  if (status) GSL_ERROR_VAL("Failed computing 0->2 period map coefficients",
			    GSL_EFAILED, GSL_NAN);

  return (p->c_0 + Iota_1 - 1.0)/Lambda_2 + 1 - p->p->c_bar;
  // return p->c_0 + Iota_1 - 1.0 + Lambda_2 - p->p->c_bar;
}

static int
tailrisk(double eps_dag_1, double eps_dag_dag_1, 
	 double k_dag_2,
	 double c_0, double L_0,
	 double *risk,
	 mexp_pars_t *p, mexp_work_t *w)
{
  int status;

  tailrisk_pars_t tp;
  tp.w = w;
  tp.p = p;
  tp.c_0 = c_0;
  tp.L_0 = L_0;
  tp.eps_dag_1 = eps_dag_1;
  tp.eps_dag_dag_1 = eps_dag_dag_1;

  /* TODO: Can we refine this bracket a bit? */
  double eps_dag_2_min = -5.0;
  double eps_dag_2_max = +5.0;
  double eps_dag_2_root;

  double y0 = tailrisk_fun(eps_dag_2_min, &tp);
  double y1 = tailrisk_fun(eps_dag_2_max, &tp);

  if (y0*y1 > 0.0)
    {
      *risk = y0 < 0 ? 1.0 : 0.0;
      return GSL_SUCCESS;
    }

  status = findroot(w, tailrisk_fun, &tp,
		    eps_dag_2_min, eps_dag_2_max,
		    &eps_dag_2_root);
  if (status) GSL_ERROR("Failed finding root when computing tailrisk",
			GSL_EFAILED);

  double Z2_root = ((k_dag_2 - sqrt(1-p->rho*p->rho)*eps_dag_2_root)/p->rho - (1-p->mu)*p->X1)/p->zeta;

  *risk = gcdf(Z2_root);

  return GSL_SUCCESS;
}


static inline double
eps_star(double k, double X, double rho)
{
  return (k - rho*X)/sqrt(1 - rho*rho);
}

double
mexp_eps_star(double k, const mexp_pars_t *p)
{
  return eps_star(k, p->X1, p->rho);
}

double
mexp_N_star(double k, const mexp_pars_t *p)
{
  return gcdf(mexp_eps_star(k, p));
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
		      double c_0,
		      double L_0, 
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
  double C2_median, Delta_q0, Delta_q1, Delta_spread, tail;

  /* Median */
  status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_median,
			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
			  p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);
  C2_median = (c_0 + Iota_1 - 1.0 + Lambda_2)*L_0;
  

  /* Profit quantile q0 */
  status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_q0,
			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
			  p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);
  Delta_q0 = (Lambda_2 - Lambda_1)*L_0;

  /* Profit quantile q1 */
  status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_q1,
			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
			  p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);
  Delta_q1 = (Lambda_2 - Lambda_1)*L_0;

  /* C2 Quantile q2 */
  /* status = mexp_map_coefs(eps_dag_1, eps_dag_dag_1, eps_dag_2_q2, */
  /* 			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2, */
  /* 			  p, w); */
  /* if (status) GSL_ERROR("Failed computing 0->2 period map coefficients", */
  /* 			GSL_EFAILED); */
  /* C2_q2 = Lambda_2*L_0; */
  status = tailrisk(eps_dag_1, eps_dag_dag_1, k_dag_2, c_0, L_0,
		    &tail,
		    p, w);
  if (status) GSL_ERROR("Failed computing tailrisk", GSL_EFAILED);

  /* Set output and we're done */

  *Omega = 
    + p->gamma0*C2_median 
    - p->gamma1*fabs(Delta_q1 - Delta_q0)
    - p->gamma2*tail
    ;

  return GSL_SUCCESS;
}

static double
impact(double k_dag_dag_1_A, double k_dag_dag_1_B,
       mexp_pars_t *p)
{
  double eps_dag_dag_1_A = eps_star(k_dag_dag_1_A, p->X1, p->rho);
  double eps_dag_dag_1_B = eps_star(k_dag_dag_1_B, p->X1, p->rho);

  double weight_A = p->L0_A/(p->L0_A + p->L0_B);
  double weight_B = p->L0_B/(p->L0_A + p->L0_B);

  return p->k0_1 + p->psi*(weight_A*gcdf(eps_dag_dag_1_A) +
			   weight_B*gcdf(eps_dag_dag_1_B));
}

int
mexp_k_dags_impacted(double k_dag_dag_1_A,
		     double k_dag_dag_1_B,
		     double *k_dag_1, double *k_dag_2,
		     mexp_pars_t *p)
{
  *k_dag_1 = impact(k_dag_dag_1_A, k_dag_dag_1_B, p);
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
				 p->c0_A, p->L0_A,
				 Omega_A, p, w);
  if (status) GSL_ERROR("Failed computing objective function A",
			GSL_EFAILED);

  status = mexp_single_objective(*k_dag_1, k_dag_dag_1_B, *k_dag_2,
				 p->c0_B, p->L0_B,
				 Omega_B, p, w);
  if (status) GSL_ERROR("Failed computing objective function B",
			GSL_EFAILED);

  return GSL_SUCCESS;				
}


double
k_fcl_ubound_fun(double k_max, void *pars)
{
  mexp_pars_t *p = pars;
  
  return k_max - p->k0_1 - p->psi*gcdf(eps_star(k_max, p->X1, p->rho));
}

double
k_fcl_A_ubound_fun(double k_max, void *pars)
{
  mexp_pars_t *p = pars;
  double k_fcl_B = p->dummy;
  return k_max - impact(k_max, k_fcl_B, p);
}

double
k_fcl_B_ubound_fun(double k_max, void *pars)
{
  mexp_pars_t *p = pars;
  double k_fcl_A = p->dummy;
  return k_max - impact(k_fcl_A, k_max, p);
}

int
mexp_find_k_fcl_ubound(double *k_fcl_ubound,
		       mexp_pars_t *p, mexp_work_t *w)
{
  int status;

  double k_max_max = p->k0_1 + p->psi;
  double k_max_min = p->k0_1;

  status = findroot(w, k_fcl_ubound_fun, p, k_max_min, k_max_max,
		    k_fcl_ubound);
  if (status) GSL_ERROR("Unable to find k impairment maximum -- this is most likely because we have multiple solutions (psi/sqrt(1-rho^2) <= 2pi appears to be a sufficient condition for unique solution)!",
			GSL_EFAILED);

  return GSL_SUCCESS;
}

int
mexp_find_k_fcl_A_ubound(double *k_fcl_A_ubound,
			 double k_fcl_B,
			 double k_fcl_ubound,
			 mexp_pars_t *p,
			 mexp_work_t *w)
{
  int status;

  p->dummy = k_fcl_B;

  double k_fcl_A_ubound_max = k_fcl_ubound;
  double k_fcl_A_ubound_min = p->k0_1;

  status = findroot(w, k_fcl_A_ubound_fun, p,
		    k_fcl_A_ubound_max,
		    k_fcl_A_ubound_min,
		    k_fcl_A_ubound);
  if (status) GSL_ERROR("Unable to find bank-A foreclosure k maximum",
			GSL_EFAILED);

  return GSL_SUCCESS;
}

int
mexp_find_k_fcl_B_ubound(double *k_fcl_B_ubound,
			 double k_fcl_A,
			 double k_fcl_ubound,
			 mexp_pars_t *p,
			 mexp_work_t *w)
{
  int status;

  p->dummy = k_fcl_A;

  double k_fcl_B_ubound_max = k_fcl_ubound;
  double k_fcl_B_ubound_min = p->k0_1;

  status = findroot(w, k_fcl_B_ubound_fun, p,
		    k_fcl_B_ubound_max,
		    k_fcl_B_ubound_min,
		    k_fcl_B_ubound);
  if (status) GSL_ERROR("Unable to find bank-B foreclosure k maximum",
			GSL_EFAILED);

  return GSL_SUCCESS;
}

double
mexp_k_of_p_tilde(double p_tilde, double k_bound,
		  mexp_pars_t *p)
{
  return p->rho*p->X1 + sqrt(1.0-p->rho*p->rho)*gicdf(p_tilde*gcdf(eps_star(k_bound, p->X1, p->rho)));
}
