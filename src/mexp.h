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

#include <math.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_integration.h>

/* NOTE: When touching this, make the corresponding change to the
 * python interface!
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

/* Main solver routines *********************************************** */

/** Compute impairment thresholds for periods 1 and 2, k_dag_1 and
 *  k_dag_2, given bank A and B foreclosure thresholds k_dag_dag_1_A
 *  and k_dag_dag_1_B.
 */
int
mexp_k_dags_impacted(double k_dag_dag_1_A,
		     double k_dag_dag_1_B,
		     double *k_dag_1, double *k_dag_2,
		     mexp_pars_t *p);

/** Compute period 0 to 1 transition coefficients Iota_1, Kappa_1,
 *  and Lambda_1, defined so that
 *
 *    L_1 = Lambda_1 L_0
 *    C_1 = C_0 + Kappa_1 L_0
 *
 *    Iota_1 = Kappa_1 - Lambda_1 + 1
 *
 *  Minimally, these can be computed with given ideosync period 1
 *  impairment, foreclosure, and period 2 impairment thresholds.
 */
int
mexp_map_coefs_period_1(double eps_dag_1, double eps_dag_dag_1,
			double eps_dag_2,
			double *Iota1, double *Kappa1, double *Lambda1,
			mexp_pars_t *p, mexp_work_t *w);

/** Compute period 1 to 2 transition coefficient Lambda_2, defined so
 *  that L_2 = Lambda_2 L_0. If needed, will also calculate Lambda_2
 *  derivative w.r.t. eps_dag_2 -- the derivative is used in
 *  evaluating the period 2 capital probability distribution
 *  functions. Set Lambda2_prime to NULL if the derivative is not
 *  needed.
 */
int
mexp_map_coefs_period_2(double eps_dag_1, double eps_dag_dag_1,
			double eps_dag_2,
			double *Lambda2, double *Lambda2_prime,
			mexp_pars_t *p, mexp_work_t *w);

/** Compute period 0 to 2 transition coefficients. Simply wraps the
 *  period 0->1 and 1->2 functions together. If Lambda2 is not needed,
 *  use the period 0->1 function, if Lambda2 derivative is not needed,
 *  set Lambda2_prime to NULL. */
int
mexp_map_coefs(double eps_dag_1, double eps_dag_dag_1,
	       double eps_dag_2,
	       double *Iota1, double *Kappa1, double *Lambda1,
	       double *Lambda2, double *Lambda2_prime,
	       mexp_pars_t *p, mexp_work_t *w);

/* TODO: Add here functions for directly computing period two capital
 * or capitalisation ratio. */

/* Compute the eps_dag_2 value corresponding to period 2 capital C2
 * (C_or_c == 0) or capitalisation ratio c2 = C2/L2 (C_or_c == 1).
 * This is needed in computing the C2 or c2 probability distribution
 * functions. */
int
mexp_find_eps_dag_2_given_c(double eps_dag_1, double eps_dag_dag_1,
			    int C_or_c, double c_0, double L_0,
			    double x_target,
			    double *eps_dag_2_root,
			    mexp_pars_t *p, mexp_work_t *w);

/* Compute capitalisation (absolute C2 if C_or_c == 0, ratio C2/L2 if
 * C_or_c == 1) of a bank (A if bank == 0, B if bank == 1) probability
 * distribution functions. This function is little more than a
 * conveniance wrapper for mexp_raw_cap2_pdf_cdf.
 */
int
mexp_cap2_pdf_cdf(double k_dag_dag_1_A,
		  double k_dag_dag_1_B,
		  int bank,
		  int C_or_c,
		  double x, double *fx, double *Fx,
		  mexp_pars_t *p, mexp_work_t *w);

/* Compute capitalisation (absolute C2 if C_or_c == 0, ratio C2/L2 if
 * C_or_c == 1) of a bank (A if bank == 0, B if bank == 1) probability
 * distribution functions. At simplest, we can do this if we are given
 * the period 1 ideosync impairment and foreclosure thresholds, and
 * the period 2 impairment threshold value k_dag_2.
 */
int
mexp_raw_cap2_pdf_cdf(double eps_dag_1, double eps_dag_dag_1, double k_dag_2,
		      int C_or_c, double c_0, double L_0,
		      double x, double *fx, double *Fx,
		      mexp_pars_t *p, mexp_work_t *w);

/* Compute the objective function. This is a minimal definition for
 * the function, mexp_two_objectives is a wrapper to evaluate Omegas
 * given just basic parameters and foreclosure choice.
 */
int
mexp_single_objective(double k_dag_1,
		      double k_dag_dag_1,
		      double k_dag_2,
		      double c_0,
		      double L_0,
		      double *Omega,
		      mexp_pars_t *p, mexp_work_t *w);

/* Compute the objective functions of banks A and B given the period 1
 * foreclosure threshold values k_dag_dag_1_A and k_dag_dag_1_B.  On
 * output also gives the period 1 and 2 impairment threshold values
 * k_dag_1 and k_dag_2.
 */
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
mexp_find_k_fcl_global_ubound(double *k_fcl_ubound,
			      mexp_pars_t *p, mexp_work_t *w);

/* Find the upper bound for a bank's foreclosure threshold value k_fcl
 * given the other bank's foreclosure percentage p_other. Need to have
 * the global k_fcl bound precomputed.
 */
int
mexp_find_k_fcl_bank_ubound(int bank, double p_other,
			    double *k_fcl_bank_ubound,
			    double k_fcl_global_ubound,
			    mexp_pars_t *p, mexp_work_t *w);

/* Given foreclosure percentages of banks A and B, compute the
 * corresponding foreclosure value thresholds k_fcl_A and k_fcl_B. To
 * compute this, one needs the upper bounds for k_fcl_A and k_fcl_B --
 * computed under the assumption that the other bank holds its p
 * constant. This is done by the function
 * mexp_find_k_fcl_bank_ubound. */
int
mexp_k_fcls_from_ps_given_bounds(double p_A, double p_B,
				 double k_fcl_A_ubound,
				 double k_fcl_B_ubound,
				 double *k_fcl_A, double *k_fcl_B,
				 mexp_pars_t *p, mexp_work_t *w);

/* Given foreclosure percentages of banks A and B, compute the
 * corresponding foreclosure value thresholds k_fcl_A and
 * k_fcl_B. Foreclosure percentage of bank A is the ratio of loans
 * foreclosed to max. number of loans impaired, given that the other
 * bank holds its foreclosure percentage fixed.
 *
 * This is just a wrapper for mexp_k_fcls_from_ps_given_bounds that
 * calls mexp_find_k_fcl_bank_ubound.
 * If large numbers of k_fcl pairs are needed, for instance, if
 * evaluating this function over a regular grid, it may be more
 * effective to use a different function. TODO: What is that function?
 */
int
mexp_k_fcls_from_ps(double p_A, double p_B,
		    double *k_fcl_A, double *k_fcl_B,
		    mexp_pars_t *p, mexp_work_t *w);


int
mexp_compute_Omega_data(int n,
			double *p_A, double *p_B,
			double *k_fcl_A_ubound,
			double *k_fcl_B_ubound,
			double *k_fcl_A, double *k_fcl_B,
			double *Omega_A, double *Omega_B,
			mexp_pars_t *p, mexp_work_t *w);


/* Map given period 1 value k threshold to corresponding ideosync
 * threshold using the period 1 shock given in parameters. */
double
mexp_eps_star(double k, const mexp_pars_t *p);

/* Convenience function: Evaluate N(eps_star(k)), where N is the
 * normal c.d.f. Uses X1 from parameters for shock.*/
double
mexp_N_star(double k, const mexp_pars_t *p);

/* Convenience function: Inverse function of mexp_N_star. */
double
mexp_N_star_inv(double N, const mexp_pars_t *p);
