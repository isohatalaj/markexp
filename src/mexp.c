/*  libmexp -- C library for solving bank market exposure model of
    Isohatala and Milne.

    Copyright (C) 2015  Jukka Isohätälä, email jukka.isohatala@gmail.com

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.           */

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

static double
exp_ff_integral(int limtype, double alpha, double omega,
		double a, double b, double c)
{
  double t2 = b*b + 1.0;
  double t = sqrt(t2);
  
  double Ialpha, Iomega;
  double Z = exp(a*(0.5*a-b*c)/t2)*gpdf(c/t);
  double Y = (a-b*c)/t;

  switch (limtype)
    {
    case MEXP_FIN_FIN:
      Ialpha = gcdf(t*alpha - Y);
      Iomega = gcdf(t*omega - Y);
      break;

    case MEXP_FIN_INF:
      Ialpha = gcdf(t*alpha - Y);
      Iomega = 1.0;
      break;

    case MEXP_INF_FIN:
      Ialpha = 0.0;
      Iomega = gcdf(t*omega - Y);
      break;

    case MEXP_INF_INF:
      Ialpha = 0.0;
      Iomega = 1.0;
    }

  return Z*(Iomega - Ialpha);
}
		
static int
LI_integrals(double eps_dag_1, double eps_dag_dag_1,
	     double eps_dag_2,
	     double *LI,
	     double *LJ,
	     mexp_pars_t *p, mexp_work_t *w)
{
  int status;
  double xi_prime = sqrt(1.0 - p->xi*p->xi);

  int lt;
  double alpha, omega;
  double a, b, c;

  lt = MEXP_FIN_FIN;
  alpha = eps_dag_dag_1;
  omega = eps_dag_1;
  a = p->phi2*(p->varphi+p->xi);
  b = -(p->varphi + p->xi)/xi_prime;
  c = (eps_dag_2 + p->varphi*eps_dag_1)/xi_prime - p->phi2*xi_prime;
  status = mexp_work_stdinteg(lt, alpha, omega,
			      a, b, c, LI, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration I returned bad status",
			GSL_EFAILED);
  if (LJ) LJ[0] = exp_ff_integral(lt, alpha, omega, a, b, c)/xi_prime;


  lt = MEXP_FIN_FIN;
  alpha = eps_dag_dag_1;
  omega = eps_dag_1;
  a = 0.0;
  b = -(p->varphi+p->xi)/xi_prime;
  c = (eps_dag_2 + p->varphi*eps_dag_1)/xi_prime;
  status = mexp_work_stdinteg(lt, alpha, omega,
			      a, b, c, LI + 1, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration II returned bad status",
			GSL_EFAILED);
  if (LJ) LJ[1] = exp_ff_integral(lt, alpha, omega, a, b, c)/xi_prime;
			      
  lt = MEXP_FIN_INF;
  alpha = eps_dag_1;
  omega = 0.0; // infty
  a = p->phi2*p->xi;
  b = -p->xi/xi_prime;
  c = eps_dag_2/xi_prime - p->phi2*xi_prime;
  status = mexp_work_stdinteg(lt, alpha, omega,
			      a, b, c, LI + 2, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration III returned bad status",
			GSL_EFAILED);
  if (LJ) LJ[2] = exp_ff_integral(lt, alpha, omega, a, b, c)/xi_prime;
		
  lt = MEXP_FIN_INF;
  alpha = eps_dag_1;
  omega = 0.0; // infty
  a = 0.0;
  b = -p->xi/xi_prime;
  c = eps_dag_2/xi_prime;
  status = mexp_work_stdinteg(lt, alpha, omega,
			      a, b, c, LI + 3, w); //ok
  if (status) GSL_ERROR("Failed computing map coefs as integration IV returned bad status",
			GSL_EFAILED);
  if (LJ) LJ[3] = exp_ff_integral(lt, alpha, omega, a, b, c)/xi_prime;

  return GSL_SUCCESS;
}

int
mexp_Lambda2(double eps_dag_1, double eps_dag_dag_1,
	     double eps_dag_2,
	     double *Lambda2, double *Lambda2_prime,
	     mexp_pars_t *p, mexp_work_t *w)
{
  double LI[4], LJ[4];
  int status; 

  status = LI_integrals(eps_dag_1, eps_dag_dag_1, eps_dag_2, 
			LI, Lambda2_prime ? LJ : NULL,
			p, w);
  if (status) GSL_ERROR("Failed doing Lambda_2 integrals",
			GSL_EFAILED);

  /* TODO: Optimisation opportunities here. */
  double A0 = exp(0.5*p->phi2*p->phi2*(1-p->xi*p->xi)-p->phi2*(eps_dag_2+p->varphi*eps_dag_1));
  double A2 = exp(0.5*p->phi2*p->phi2*(1-p->xi*p->xi)-p->phi2*eps_dag_2);
  double B0 = 1.0 - gcdf(eps_dag_dag_1);

  *Lambda2 = B0 + A0*LI[0] - LI[1] + A2*LI[2] - LI[3];

  if (Lambda2_prime)
    {
      *Lambda2_prime = 
	A0*(-p->phi2*LI[0] + LJ[0]) 
	- LJ[1] 
	+ A2*(-p->phi2*LI[2] + LJ[2])
	- LJ[3];
    }

  return GSL_SUCCESS;
}

int
mexp_map_coefs_1(double eps_dag_1, double eps_dag_dag_1,
		 double eps_dag_2,
		 double *Iota1, double *Kappa1, double *Lambda1,
		 mexp_pars_t *p, mexp_work_t *w)
{
  /* TODO: Optimisation opportunities here. */

  *Iota1   = exp(p->phi1*(0.5*p->phi1-eps_dag_1))*gcdf(eps_dag_dag_1-p->phi1); //ok
  *Kappa1  = -gcdf(eps_dag_1) + exp(p->phi1*(0.5*p->phi1-eps_dag_1))*gcdf(eps_dag_1-p->phi1); //ok
  *Lambda1 = 1.0 - gcdf(eps_dag_1) - exp(p->phi1*(0.5*p->phi1-eps_dag_1))*(gcdf(eps_dag_dag_1 - p->phi1) - gcdf(eps_dag_1 - p->phi1)); //ok

  return GSL_SUCCESS;
}

int
mexp_map_coefs(double eps_dag_1, double eps_dag_dag_1,
	       double eps_dag_2,
	       double *Iota1, double *Kappa1, double *Lambda1,
	       double *Lambda2,
	       mexp_pars_t *p, mexp_work_t *w)
{
  int status;

  status = mexp_map_coefs_1(eps_dag_1, eps_dag_dag_1, eps_dag_2,
			    Iota1, Kappa1, Lambda1,
			    p, w);
  if (status) GSL_ERROR("Failed computing period 1 IKL coefficients",
			GSL_EFAILED);
  
  status = mexp_Lambda2(eps_dag_1, eps_dag_dag_1, eps_dag_2,
			Lambda2, NULL, p, w);
  if (status) GSL_ERROR("Failed computing Lambda2",
			GSL_EFAILED);

  return GSL_SUCCESS;
}

typedef struct {
  int C_or_c;
  mexp_work_t *w;
  mexp_pars_t *p;
  double c_0;
  double L_0;
  double eps_dag_1;
  double eps_dag_dag_1;
  double x_target;
} dag2finder_pars_t;

static double
dag2finder_fun(double eps_dag_2, void *pars)
{
  int status;
  dag2finder_pars_t *p = pars;
  double Iota_1, Kappa_1, Lambda_1, Lambda_2;

  status = mexp_map_coefs(p->eps_dag_1, p->eps_dag_dag_1, eps_dag_2,
  			  &Iota_1, &Kappa_1, &Lambda_1, &Lambda_2,
  			  p->p, p->w);
  if (status) GSL_ERROR_VAL("Failed computing 0->2 period map coefficients",
			    GSL_EFAILED, GSL_NAN);

  if (p->C_or_c == 0)
    {
      return (p->c_0 + Iota_1 - 1.0 + Lambda_2)*p->L_0 - p->x_target;
    }
  else
    {
      return (p->c_0 + Iota_1 - 1.0)/Lambda_2 + 1 - p->x_target;
    }
}

int
mexp_find_eps_dag_2_given_c(double eps_dag_1, double eps_dag_dag_1,
			    int C_or_c, double c_0, double L_0,
			    double x_target,
			    double *eps_dag_2_root,
			    mexp_pars_t *p, mexp_work_t *w)
{
  int status;

  dag2finder_pars_t tp;
  tp.C_or_c = C_or_c;
  tp.w = w;
  tp.p = p;
  tp.c_0 = c_0;
  tp.L_0 = L_0;
  tp.eps_dag_1 = eps_dag_1;
  tp.eps_dag_dag_1 = eps_dag_dag_1;
  tp.x_target = x_target;

  /* TODO: Can we refine this bracket a bit? */
  double eps_dag_2_min = -8.0;
  double eps_dag_2_max = +8.0;

  double y0 = dag2finder_fun(eps_dag_2_min, &tp);
  double y1 = dag2finder_fun(eps_dag_2_max, &tp);

  if (y0*y1 > 0.0)
    {
      *eps_dag_2_root = y0 > 0 ? eps_dag_2_max : eps_dag_2_min;
    }
  else
    {
      status = findroot(w, dag2finder_fun, &tp,
			eps_dag_2_min, eps_dag_2_max,
			eps_dag_2_root);
      if (status) GSL_ERROR("Failed finding root when computing tailrisk",
			    GSL_EFAILED);
    }

  return GSL_SUCCESS;
}

static double
eps_dag_2_to_Z2(double k_dag_2, double eps_dag_2, mexp_pars_t *p)
{
  return ((k_dag_2 - sqrt(1-p->rho*p->rho)*eps_dag_2)/p->rho - (1-p->mu)*p->X1)/p->zeta;
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
  double C2_median, Delta_q0, Delta_q1, tail;

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

  /* P of c2 being under cbar */
  double dummy;
  status = mexp_cr2_pdf_cdf(eps_dag_1, eps_dag_dag_1, k_dag_2,
			    1, c_0, L_0, p->c_bar, &dummy, &tail, p, w);
  if (status) GSL_ERROR("Failed computing tailrisk", GSL_EFAILED);

  /* Set output and we're done */

  *Omega = 
    + p->gamma0*C2_median 
    - p->gamma1*fabs(Delta_q1 - Delta_q0)
    - p->gamma2*tail
    ;

  return GSL_SUCCESS;
}

int
mexp_cr2_pdf_cdf(double eps_dag_1, double eps_dag_dag_1, double k_dag_2,
		 int C_or_c, double c_0, double L_0,
		 double x, double *fx, double *Fx,
		 mexp_pars_t *p, mexp_work_t *w)
{
  int status;


  double eps_dag_2_root;

  status = mexp_find_eps_dag_2_given_c(eps_dag_1, eps_dag_dag_1,
  				       C_or_c, c_0, L_0, x, &eps_dag_2_root, 
				       p, w);
  if (status) GSL_ERROR("Could not find eps_dag_2 corresponding to given "
  			"period 2 capitalisation",
  			GSL_EFAILED);

  double Iota_1, Kappa_1, Lambda_1;

  status = mexp_map_coefs_1(eps_dag_1, eps_dag_dag_1, eps_dag_2_root,
			    &Iota_1, &Kappa_1, &Lambda_1,
			    p, w);
  if (status) GSL_ERROR("Failed computing 0->2 period map coefficients",
			GSL_EFAILED);

  double Lambda_2, Lambda_2_prime;

  status = mexp_Lambda2(eps_dag_1, eps_dag_dag_1, eps_dag_2_root,
			&Lambda_2, &Lambda_2_prime, p, w);
  if (status) GSL_ERROR("Failed computing Lambda 2 and its derivative",
			GSL_EFAILED);

  double x_prime;
  
  // c2 = (p->c_0 + Iota_1 - 1.0)/Lambda_2 + 1;
  // C2 = c_0*L_0 + (Io1 - 1 + Lam2)*L0
  if (C_or_c == 0)
    {
      x_prime = L_0*Lambda_2_prime;
    }
  else
    {
      x_prime = -(c_0 + Iota_1 - 1.0)*Lambda_2_prime/(Lambda_2*Lambda_2);
    }

  double dx_dZ2;
  double deps_2_dag_dZ2 = -p->rho*p->zeta/sqrt(1.0-p->rho*p->rho);
  
  dx_dZ2 = x_prime * deps_2_dag_dZ2;

  double Z2_root;

  Z2_root = eps_dag_2_to_Z2(k_dag_2, eps_dag_2_root, p);

  double fZ = gpdf(Z2_root);
  
  *fx = fZ/fabs(dx_dZ2);
  *Fx = gcdf(Z2_root);

  return GSL_SUCCESS;
}

int
mexp_cap2_pdf_cdf(double k_dag_dag_1_A,
		  double k_dag_dag_1_B,
		  int bank,
		  int C_or_c,
		  double x, double *fx, double *Fx,
		  mexp_pars_t *p, mexp_work_t *w)
{
  double eps_dag_1, eps_dag_dag_1;
  double k_dag_dag_1_X;
  double k_dag_1;
  double k_dag_2;
  double c0, L0;
  
  mexp_k_dags_impacted(k_dag_dag_1_A,
		       k_dag_dag_1_B,
		       &k_dag_1, &k_dag_2,
		       p);

  eps_dag_1 = eps_star(k_dag_1, p->X1, p->rho);

  if (bank == 0)
    {
      k_dag_dag_1_X = k_dag_dag_1_A;
      c0 = p->c0_A;
      L0 = p->L0_A;
    }
  else
    {
      k_dag_dag_1_X = k_dag_dag_1_B;
      c0 = p->c0_B;
      L0 = p->L0_B;
    }

  eps_dag_dag_1 = eps_star(k_dag_dag_1_X, p->X1, p->rho);

  return mexp_cr2_pdf_cdf(eps_dag_1, eps_dag_dag_1, k_dag_2,
			  C_or_c, c0, L0,
			  x, fx, Fx,
			  p, w);
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
  if (status) GSL_ERROR("Unable to find k impairment maximum -- "
			"this is most likely because we have multiple "
			"solutions.",
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
  return p->rho*p->X1 
    + sqrt(1.0-p->rho*p->rho)
      *gicdf(p_tilde*gcdf(eps_star(k_bound, p->X1, p->rho)));
}
