
#include <stdio.h>
#include <math.h>

#include "bnet.h"

/**
 * Shorthands for the normal distribution pdf, cdf, and inverse cdf.
 */
inline double  gpdf(double x) { return gsl_ran_ugaussian_pdf(x);  }
inline double  gcdf(double x) { return gsl_cdf_ugaussian_P(x);    }
inline double gicdf(double x) { return gsl_cdf_ugaussian_Pinv(x); }

/**
 * Custom error handler. We use GSLs error handling system, with this
 * as our error handler. Unlike the default, this does not call abort,
 * so we (1) have a chance to recover and (2) we get a trace of the
 * call path that lead to the error condition (useful for debugging).
 */
void 
errhandler(const char * reason, 
	   const char * file, 
	   int line, 
	   int gerrno)
{
  fprintf(stderr, "gsl error @ %s(%d): %s [%d]\n", 
	  file, line, reason, gerrno);
}


bnet_work_t *
bnet_work_alloc()
{
  bnet_work_t *self;

  self = malloc(sizeof(*self));
  if (self == NULL) return NULL;

  /* All uninitialized pointers to NULL -- error case cleanup code,
   * bnet_work_free, frees all non-NULL pointers. */
  self->solver   = NULL;
  self->solver_2 = NULL;
  self->intwork  = NULL;

  /* Default tolerances, currently no "official" way to chance these.
   * These choices should be a rather reasonable first guess, though.
   * TODO: Add functions to adjust these.
   */
  self->tol_rel  = 1e-6;
  self->tol_abs  = 1e-6;
  self->max_iter = 64;
  self->intwork_size = 512;

  /* gsl object initializations */
  self->solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  if (self->solver == NULL) goto fail;

  self->solver_2 = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  if (self->solver_2 == NULL) goto fail;

  self->intwork = gsl_integration_workspace_alloc(self->intwork_size);
  if (self->intwork == NULL) goto fail;


  /* Successful exit. */
  return self;

 fail:
  /* Failure exit. All allocated objects are deleted. */
  bnet_work_free(self);
  return NULL;
}

/* Helper function. Initializes either the first or second root solver
 * in bnet_work_t, and findroots a given function. */
static int
findroot_x(bnet_work_t *self, 
	   int s, /*< Zero for first solver, one for second. */
	   double (*fun)(double, void*), 
	   void *params, 
	   double x0, double x1,
	   double *root)
{
  int status;

  gsl_function f;
  f.function = fun;
  f.params   = params;

  gsl_root_fsolver *solver = s == 0 ? self->solver : self->solver_2;

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

  status = gsl_root_fsolver_set(solver, &f, x0, x1);
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
      status = gsl_root_fsolver_iterate(solver);
      if (status) GSL_ERROR("root solver iteration failed", status);
      x = gsl_root_fsolver_root(solver);
      status = gsl_root_test_interval(gsl_root_fsolver_x_lower(solver), 
				      gsl_root_fsolver_x_upper(solver),
				      self->tol_abs, 
				      self->tol_rel);
    }
  while (status == GSL_CONTINUE && iter < self->max_iter);

  if (status) GSL_ERROR("root solver iterations failed to converge", 
			status);

  *root = x;

  return GSL_SUCCESS;
}

/* Conveniance wrappers for findroot_x. */
static inline int
findroot(bnet_work_t *self, 
	 double (*fun)(double, void*), void *params, 
	 double x0, double x1, double *root)
{
  return findroot_x(self, 0, fun, params, x0, x1, root);
}

static inline int
findroot_2(bnet_work_t *self, 
	   double (*fun)(double, void*), void *params, 
	   double x0, double x1, double *root)
{
  return findroot_x(self, 1, fun, params, x0, x1, root);
}

void
bnet_work_free(bnet_work_t *self)
{
  if (self == NULL) return;
  if (self->solver)   gsl_root_fsolver_free(self->solver);
  if (self->solver_2) gsl_root_fsolver_free(self->solver_2);
  if (self->intwork)  gsl_integration_workspace_free(self->intwork);

  free(self);
}

/* Test if x is y to within eps_rel and eps_abs relative or absolute
 * tolerance. */
static int
fcmp(double x, double y, double eps_rel, double eps_abs)
{
  return fabs(x - y) < eps_rel * (fabs(x) + fabs(y)) + eps_abs;
}

/* Baseline impairment k threshold, given unconditional default rate
 * PD. */
double
bnet_kzero(double PD)
{
  return gicdf(PD);
}

/* Conversion function from some k threshold to ideosync. threshold
 * eps, for given aggregate shock X and correlation rho. */
double
bnet_eps_star(double k, double X, double rho)
{
  return (k - rho * X) / sqrt(1 - rho*rho);
}

/* Inverse function of eps_star -- converts a given ideosync.
 * threshold eps into a k value. */
static double
k_star(double eps, double X, double rho)
{
  return sqrt(1 - rho*rho) * eps + rho * X;
}

static double
X_star(double k, double eps, double rho)
{
  return (k - sqrt(1 - rho*rho) * eps) / rho;
}


/* Since eps's and k's are interchangable, we should for simplicity
 * stick to using only one or the other. This will avoid unnecessary
 * complexity, and possible errors. Convention: FAVOUR eps
 * REPRESENTATION OVER k. */

/* Compute period 0 to 1 mapping coefficients kappa, lambda, and iota
 * given impairment threshold eps, foreclosure threshold epsfc, and
 * the conversion factor phi. */
static int
kli(double phi, double epsfc, double eps,
    double *kappa, double *lambda, double *iota)
{
  const double g0 = gcdf(eps);
  const double g1 = gcdf(eps - phi);
  const double g2 = gcdf(epsfc - phi);

  const double h0 = exp(0.5 * phi * (phi - 2.0 * eps));

  *kappa  = g0 - h0 * g1;
  *lambda = g0 + h0 * (g2 - g1);
  *iota   = h0 * g2;

  return 0;
}

/* Expect period 1 loan revaluation factor. Multiply your period zero
 * loan portfolio accounting value by this to get the period 1
 * value. */
static double
comp_alphaone(double epsimp, double epsfc, double phi)
{
  return
    exp(phi*(0.5*phi - epsimp)) 
      * (gcdf(epsimp - phi) - gcdf(epsfc - phi)) 
    + 1.0 - gcdf(epsimp);
}

/* Next three functions are the expected period two revaluation
 * factors. The three different expressions correspond to three
 * qualitatively different outcomes of the period 2 aggregate shock.
 */
static double
comp_alphaminus(double epsimp, double epsimp2, double epsfc,
		double phi, double phi2)
{
  return
    exp(-phi*(epsimp2 - epsimp) + phi2*(0.5*phi2 - epsimp)) 
      * (gcdf(epsimp - phi2) - gcdf(epsfc - phi2))
    + exp(phi*(0.5*phi - epsimp2)) 
      * (gcdf(epsimp2 - phi) - gcdf(epsimp - phi))
    + 1.0 - gcdf(epsimp2);
}

static double
comp_alphaplus(double epsimp, double epsimp2, double epsfc,
	       double phi, double phi2)
{
  return
    exp(phi2*(0.5*phi2 - epsimp2)) 
      * (gcdf(epsimp2 - phi2) - gcdf(epsfc - phi2))
    + 1.0 - gcdf(epsimp2);
}

static inline double
comp_alphaplusplus(double epsfc)
{
  return
    1.0 - gcdf(epsfc);
}

/* Period 2 expected revaluation, wrapped in single neat package.*/
static double
comp_alphatwo(double epsimp, double epsimp2, double epsfc,
	      double phi, double phi2)
{
  if (epsimp < epsimp2)
    {
      return comp_alphaminus(epsimp, epsimp2, epsfc, phi, phi2);
    }
  
  if (epsfc < epsimp2 && epsimp2 <= epsimp)
    {
      return comp_alphaplus(epsimp, epsimp2, epsfc, phi, phi2);
    }

  return comp_alphaplusplus(epsfc);
}

/* We need to rootfind the the function alpha2, and in order to do
 * that using the gsl rootfinders, we need to wrap that function into
 * a gsl function, with the static parameters packed into a single
 * struct. Which is here: */
typedef struct {
  double epsimp;
  double epsfc;
  double phi;
  double phi2;
  double delta;

  double power;

  double X1;
  double theta;

  /* these two are only needed by alphatwo_expect */
  double kimp2; /*< Period 2 impairment k threshold. */
  double rho;
} alpha_pars_t;

/* Alpha2 times normal pdf with parameters deferred to the
 * alpha_pars_t struct, passed in as a voidptr. */
static double
comp_alphatwo_expect_gfun(double Z2, void *params)
{
  alpha_pars_t *p = params;
  double X2 = p->theta*p->X1 + Z2;
  double epsimp2 = bnet_eps_star(p->kimp2, X2, p->rho);
  double alphatwo = comp_alphatwo(p->epsimp, epsimp2, p->epsfc, 
				  p->phi, p->phi2);

  if (p->power != 1.0) return gpdf(Z2) * pow(alphatwo, p->power);
  return gpdf(Z2) * alphatwo;
}

/* Compute period 1 to 2 profit, can be negative. */
static double
comp_Delta(double epsimp, double epsimp2, double epsfc,
	   double phi, double phi2)
{
  double aone = comp_alphaone(epsimp, epsfc, phi);
  double atwo = comp_alphatwo(epsimp, epsimp2, epsfc, 
			      phi, phi2);

  return atwo - aone;
}

/* Compute maximum profit. This occurs when the second period
 * impairment threshold epsimp2 turns out to be at or below the
 * foreclosure level. */
double
comp_Delta_max(double epsimp, double epsfc, double phi, double phi2)
{
  /* Max. profit when epsimp2 >= epsfc: */
  return comp_Delta(epsimp, epsfc, epsfc, phi, phi2);
}

/* Compute minimum profit, will be negative. This occurs in the rather
 * unlikely scenario when literally all, every single loan defaults =>
 * net revaluation factor alpha2 = 0. */
double
comp_Delta_min(double epsimp, double epsfc, double phi)
{
  return -comp_alphaone(epsimp, epsfc, phi);
}

/* Compute the probability mass at the maximum delta point. Since a
 * good number of X2 outcomes correspond to the same (maximum) profit,
 * this is a nonzero value. */
double
comp_Delta_max_mass(double rho, 
		    double X1,
		    double kimp2, double epsimp, double epsfc, 
		    double phi, double phi2, double theta)
{
  double X2 = X_star(kimp2, epsfc, rho);
  double Z2 = X2 - theta*X1;
  double bulk_mass = gcdf(Z2);

  return 1.0 - bulk_mass;
}

/* Wrapper function to be supplied for a rootfinder -- used for
 * finding epsimp2 such that profit is at given delta. */
double
comp_Delta_diff_gfun(double epsimp2, alpha_pars_t *p)
{
  return 
    comp_Delta(p->epsimp, epsimp2, p->epsfc, p->phi, p->phi2) 
    - p->delta;
}

/* Compute the profit cdf. See the main text for the maths background
 * on this is done. The only stochastic factor in the profit is X2,
 * the outturn of period 2 aggregate shock. Therefore, P[profit(X2) <
 * delta] conditioned on X2, P[profit(X2) < delta | X2], is a step
 * function (profit is decreasing in X2). We of course integrate the
 * conditional probability against X2 pdf, and since it is a step
 * function and X2 is normal, that is easy. The only thing that
 * requires some work is finding the X2 point where profit(X2) ==
 * delta. This is what we are rootfinding here.
 */
int
comp_Delta_cdf(double delta,
	       double rho, 
	       double X1,
	       double kimp2,
	       double epsimp, double epsfc,
	       double phi, double phi2,
	       double theta,
	       double *FDelta,
	       bnet_work_t *work)
{
  double Delta_max = comp_Delta_max(epsimp, epsfc, phi, phi2);
  double Delta_min = comp_Delta_min(epsimp, epsfc, phi);

  if (delta >= Delta_max)
    {
      *FDelta = 1.0;
      return GSL_SUCCESS;
    }

  if (delta <= Delta_min)
    {
      *FDelta = 0.0;
      return GSL_SUCCESS;
    }

  /* First do a coarse bracketing of the root */
  
  double x0 = 0;
  double dx = 0.5;
  double x1 = x0 + dx;
  double y0 = Delta_max - delta;
  double y1 = comp_Delta(epsimp, epsfc + x1, epsfc, phi, phi2) - delta;

  int iter = 0;

  while (y1 > 0.0)
    {
      y0 = y1;
      x0 = x1;

      x1 += dx;
      dx *= 1.5;
      y1 = comp_Delta(epsimp, epsfc + x1, epsfc, phi, phi2) - delta;
      
      iter ++;

      if (iter == work->max_iter) 
	GSL_ERROR("cannot locate root Delta(epsimp2) == delta", 
		  GSL_EFAILED);
    }


  /* Now polish the root */

  double eps2_root;

  /* check if we happen to sit on the root already (needs to be done,
   * it's a typical occurance and can cause the root finder to
   * fail). */
  if (fcmp(y0, 0.0, work->tol_rel, work->tol_abs))
    {
      eps2_root = epsfc;
    }
  else
    {
      alpha_pars_t p;
      p.epsimp = epsimp;
      p.epsfc  = epsfc;
      p.phi    = phi;
      p.phi2   = phi2;
      p.delta  = delta;
      
      int status;
      status = findroot(work, 
			(double (*)(double, void *)) 
			comp_Delta_diff_gfun, &p,
			epsfc + x0, epsfc + x1, &eps2_root);
      
      if (status)
	{
	  GSL_ERROR("failed polishing Delta(epsimp2) == delta root", 
		    GSL_EFAILED);
	}
    }

  double X2 = (-eps2_root * sqrt(1.0 - rho*rho) + kimp2) / rho;
  double Z2 = X2 - theta*X1;
  *FDelta = gcdf(Z2);
  
  return GSL_SUCCESS;
}


/* To calculate interpercentile ranges, or rather more generally
 * quantiles, we need to be able to root find the profit distribution
 * function. So, tediuously, we need to weap comp_Delta_cdf into a gsl
 * function to be passed to the root finders. Gah.
 */
typedef struct {
  double quant;
  double rho;
  double X1;
  double kimp2;
  double epsimp;
  double epsfc;
  double phi;
  double phi2;
  double theta;
  bnet_work_t *work;
} delta_cdf_pars_t;

/* This is just comp_Delta_cdf - quant but with static parameters
 * passed in via the null ptr. */
double
comp_Delta_cdf_gfun(double delta, void *params)
{
  delta_cdf_pars_t *p = params;
  int status;
  double F;

  status = comp_Delta_cdf(delta, p->rho, p->X1, p->kimp2, p->epsimp, p->epsfc,
			  p->phi, p->phi2, p->theta, &F, p->work);
  if (status) return GSL_NAN;

  return F - p->quant;
}

/* Compute a profit quantile. Since we have our Delta cdf function
 * ready and wrapped, all we need to do is some simple rootfinding. */
int
comp_Delta_quantile(double quant,
		    double rho, 
		    double X1,
		    double kimp2,
		    double epsimp, double epsfc,
		    double phi, double phi2,
		    double theta,
		    double *deltap,
		    bnet_work_t *work)
{
  double Delta_max  = comp_Delta_max(epsimp, epsfc, phi, phi2);
  double Delta_min  = comp_Delta_min(epsimp, epsfc, phi);
  double Delta_qmax = 1.0 - comp_Delta_max_mass(rho, X1,
						kimp2, epsimp, epsfc, 
						phi, phi2, theta);

  if (quant <= 0.0)        { *deltap = Delta_min; return GSL_SUCCESS; }
  if (quant >= Delta_qmax) { *deltap = Delta_max; return GSL_SUCCESS; }

  delta_cdf_pars_t p;

  p.quant  = quant;
  p.rho    = rho;
  p.X1     = X1;
  p.kimp2  = kimp2;
  p.epsimp = epsimp;
  p.epsfc  = epsfc;
  p.phi    = phi;
  p.phi2   = phi2;
  p.theta  = theta;
  p.work   = work;

  int status;
  status = findroot_2(work, comp_Delta_cdf_gfun, &p,
		      Delta_min, Delta_max, deltap);
  if (status)
    {
      double y0, y1;

      comp_Delta_cdf(Delta_min, rho, X1, kimp2, epsimp, epsfc, 
		     phi, phi2, theta, &y0, work);
      comp_Delta_cdf(Delta_max, rho, X1, kimp2, epsimp, epsfc, 
		     phi, phi2, theta, &y1, work);

      fprintf(stderr, "# error debug: y0 = %le, y1 = %le\n", y0, y1);
      GSL_ERROR("failed finding root when determining Delta quantile", 
		GSL_EFAILED);
    }

  return GSL_SUCCESS; 
}

/* Compute the expectation value of alpha2, where the expectation is
 * naturally to be taken over the outturns of Z2. */
int
comp_alphatwo_expectation(double rho, 
			  double X1,
			  double kimp2,
			  double epsimp, double epsfc,
			  double phi, double phi2,
			  double theta,
			  double power,
			  double *a2expect,
			  bnet_work_t *w)
{
  int status;
  alpha_pars_t p;

  p.epsimp = epsimp;
  p.epsfc  = epsfc;
  p.phi    = phi;
  p.phi2   = phi2;
  p.kimp2  = kimp2;
  p.rho    = rho;
  p.X1     = X1;
  p.theta  = theta;
  p.power  = power;

  gsl_function f;

  f.function = comp_alphatwo_expect_gfun;
  f.params   = &p;

  double a2_mean_err;
  double a2_mean;

  status = gsl_integration_qagi(&f,
  				w->tol_abs, w->tol_rel,
  				w->intwork_size, w->intwork,
  				&a2_mean, &a2_mean_err);
  if (status) GSL_ERROR("integration failed", GSL_EFAILED);
  
  /* TODO: The above integral can be simplified a bit. */
    
  *a2expect = a2_mean;

  return GSL_SUCCESS;  
}

/* phi calibration code ********************************************************* */
/* A boring piece of code. */

typedef struct {
  double k;
  double rho;
  double phi;
} calib_phi_fun_integ_pars_t;

inline double 
comp_lambda(double phi, double eps)
{
  return gcdf(eps) - exp(0.5 * phi * (phi - 2.0 * eps)) * gcdf(eps - phi);
}

double
calib_phi_fun_integ_fun(double x, void *params)
{
  calib_phi_fun_integ_pars_t *p = params;
  double eps = bnet_eps_star(p->k, x, p->rho);
  return gpdf(x) * comp_lambda(p->phi, eps);
}

typedef struct {
  double PD;
  double LGD;
  double rho;
  double k;
  bnet_work_t *work;
} calib_phi_pars_t;

double 
calib_phi_fun(double phi, void *params)
{
  int status;
  calib_phi_pars_t *p = params;
  gsl_function f;
  calib_phi_fun_integ_pars_t ip;
  
  ip.k   = p->k;
  ip.rho = p->rho;
  ip.phi = phi;

  f.function = calib_phi_fun_integ_fun;
  f.params   = &ip;

  double pd, pderr;

  status = gsl_integration_qagi(&f, 
				p->work->tol_abs, p->work->tol_rel, 
				p->work->intwork_size, p->work->intwork, 
				&pd, &pderr);
  if (status) GSL_ERROR_VAL("integration failed", status, GSL_NAN);

  return pd - p->PD*p->LGD;
}

int
calib_phi(double LGD, double PD, double rho, double *phires, 
	  bnet_work_t *work)
{
  calib_phi_pars_t p;
  p.PD   = PD;
  p.LGD  = LGD;
  p.rho  = rho;
  p.k    = bnet_kzero(PD);
  p.work = work;

  double phi;
  double phi0 = 0.0, phi1 = 0.5;
  double fphi0 = 0.0; // calib_phi_fun(phi0, &p);
  double fphi1 = calib_phi_fun(phi1, &p);
  
  while (fphi0 * fphi1 >= 0.0)
    {
      phi0 = phi1;
      fphi0 = fphi1;

      phi1 *= 2.0;
      fphi1 = calib_phi_fun(phi1, &p);
    }

  int status;
  status = findroot(work, calib_phi_fun, &p,
		    phi0, phi1, &phi);
  if (status) GSL_ERROR("root finder failed", status);

  *phires = phi;
  return GSL_SUCCESS;
}

int
bnet_calibrate_phi(bnet_work_t *work,
		   bnet_pars_t *pars,
		   double *phi_calib)
{
  int status;

  status = calib_phi(pars->LGD, pars->PD, pars->rho, phi_calib, work);
  if (status) GSL_ERROR("failed calibrating phi", GSL_EFAILED);

  return GSL_SUCCESS;
}


/* ****************************************************************************** */

/* Calculate the impairment threshold with price interactions. */
double
bnet_comp_kimp(const double epsfc[2], const double L[2], double X, 
	       double k0, double psi)
{
  return k0 + psi * (L[0] * gcdf(epsfc[0]) + L[1] * gcdf(epsfc[1])) 
                  / (L[0] + L[1]);
}

double
bnet_comp_kimp2(double kimp, const bnet_pars_t *p)
{
  return (1 - p->chi)*bnet_kzero(p->PD) + p->chi*kimp;
}

int
bnet_comp_Omega(double c0,
		double X,
		double kfc, double kimp,
		double *Omega, 
		double *c1_out, 
		double *Ec2_out,
		double *ipr_out,
		bnet_pars_t *p,
		bnet_work_t *w)
{
  int status;
  double kappa1, lambda1, iota1;
  double c1;
  double kimp2 = bnet_comp_kimp2(kimp, p);
  double epsfc = bnet_eps_star(kfc, X, p->rho);
  double epsimp = bnet_eps_star(kimp, X, p->rho);
  double epsimp2 = bnet_eps_star(kimp2, X, p->rho);
  
  status = kli(p->phi, epsfc, epsimp, &kappa1, &lambda1, &iota1);
  if (status) 
    GSL_ERROR("failed computing kli coefficients", GSL_EFAILED);

  c1 = (c0 - kappa1) / (1 - lambda1);

  /* */
  double a1;
  a1 = comp_alphaone(epsimp, epsfc, p->phi);

  double a2expect;
  status = comp_alphatwo_expectation(p->rho, X, kimp2, epsimp, epsfc, 
				     p->phi, p->phi2, p->theta,
				     1.0,
				     &a2expect,
				     w);
  if (status) 
    GSL_ERROR("failed computing alpha2 expectation value", 
	      GSL_EFAILED);

  /* double a2invexpect; */
  /* status = comp_alphatwo_expectation(p->rho, X, kimp2, epsimp, epsfc,  */
  /* 				     p->phi, p->phi2, p->theta, */
  /* 				     -1.0, */
  /* 				     &a2expect, */
  /* 				     w); */
  /* if (status)  */
  /*   GSL_ERROR("failed computing alpha2 inverse expectation value",  */
  /* 	      GSL_EFAILED); */

  /* double delta_median; */
  /* status = comp_Delta_quantile(0.5, p->rho, X, kimp2, epsimp, epsfc,  */
  /* 			       p->phi, p->phi2, p->theta, &delta_median, w); */
  /* if (status) GSL_ERROR("failed finding delta median", GSL_EFAILED); */

  
  double deltaq0, deltaq1;

  status = comp_Delta_quantile(p->q0, p->rho, X, kimp2, epsimp, epsfc, 
			       p->phi, p->phi2, p->theta, &deltaq0, w);
  if (status) GSL_ERROR("failed finding delta quantile q0", GSL_EFAILED);

  status = comp_Delta_quantile(p->q1, p->rho, X, kimp2, epsimp, epsfc, 
			       p->phi, p->phi2, p->theta, &deltaq1, w);
  if (status) GSL_ERROR("failed finding delta quantile q1", GSL_EFAILED);

  double ipr = deltaq1 - deltaq0;
  double EC2overL0 = c0 + iota1 + a2expect - 1.0;
  double EDelta = a2expect - a1;
  /* double Ec2 = (c0 + iota1 - 1.0) * a2invexpect + 1.0; */

  double risk = ipr;
  /* double gain = a2expect;  */
  double gain = EC2overL0;

  *Omega = (1 - p->gamma) * gain - p->gamma * risk;

  *c1_out = c1;
  *ipr_out = risk;
  *Ec2_out = gain;

  return GSL_SUCCESS;
}

int
bnet_comp_Omegas(double C[2],
		 double L[2],
		 double X, 
		 double kfc[2],
		 double Omega[2],
		 bnet_pars_t *p,
		 bnet_work_t *w)
{
  int status;
  double epsfc[2];
  double k0, kimp;

  k0 = bnet_kzero(p->PD);

  epsfc[0] = bnet_eps_star(kfc[0], X, p->rho);
  epsfc[1] = bnet_eps_star(kfc[1], X, p->rho);

  kimp = bnet_comp_kimp(epsfc, L, X, k0, p->psi);

  double c1_out, Ec2_out, ipr_out;

  status = bnet_comp_Omega(C[0] / L[0],
			   X, kfc[0], kimp,
			   &Omega[0], &c1_out, &Ec2_out, &ipr_out,
			   p, w);
  if (status) GSL_ERROR("Failed evaluating OmegaA", GSL_EFAILED);

  status = bnet_comp_Omega(C[1] / L[1],
			   X, kfc[1], kimp,
			   &Omega[1], &c1_out, &Ec2_out, &ipr_out,
			   p, w);
  if (status) GSL_ERROR("Failed evaluating OmegaB", GSL_EFAILED);

  return GSL_SUCCESS;  
}

