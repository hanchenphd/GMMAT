/*  GMMAT : An R Package for Generalized linear Mixed Model Association Tests
 *  Copyright (C) 2014--2020  Han Chen, Matthew P. Conomos, Duy T. Pham
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*  Function use_brent was modified from Brent_fmin function
 *  by Han Chen on 03/11/2014, changes marked.
 *  Brent_fmin is the C function used by R function optimize.
 *  The copyright information of the original Brent_fmin function
 *  is shown below.
 */

/*
 *  R : A Computer Language for Statistical Data Analysis
 *  Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 *  Copyright (C) 2003-2004  The R Foundation
 *  Copyright (C) 1998--2013  The R Core Team
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  http://www.r-project.org/Licenses/
 */

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <zlib.h>
#include <bzlib.h>
#define STRICT_R_HEADERS
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include "read_bgen.h"
#include "zstd-1.4.5/lib/zstd.h"
#include "libdeflate-1.6/libdeflate.h"
using namespace std;
using namespace arma;
using namespace Rcpp;


typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

#define DBL_EPSILON 2.2204460492503131e-16;





/* Function use_brent was modified from Brent_fmin function
   by Han Chen on 03/11/2014, changes marked below.
   Brent_fmin is the C function used by R function optimize. */

/* Formerly in src/appl/fmim.c */

/* fmin.f -- translated by f2c (version 19990503).
*/

/* R's  optimize() :   function	fmin(ax,bx,f,tol)
   =    ==========		~~~~~~~~~~~~~~~~~

        an approximation  x  to the point where  f  attains a minimum  on
    the interval  (ax,bx)  is determined.

    INPUT..

    ax    left endpoint of initial interval
    bx    right endpoint of initial interval
    f     function which evaluates  f(x, info)  for any  x
          in the interval  (ax,bx)
    tol   desired length of the interval of uncertainty of the final
          result ( >= 0.)

    OUTPUT..

    fmin  abcissa approximating the point where  f  attains a minimum

        The method used is a combination of  golden  section  search  and
    successive parabolic interpolation.  convergence is never much slower
    than  that  for  a  Fibonacci search.  If  f  has a continuous second
    derivative which is positive at the minimum (which is not  at  ax  or
    bx),  then  convergence  is  superlinear, and usually of the order of
    about  1.324....
        The function  f  is never evaluated at two points closer together
    than  eps*abs(fmin)+(tol/3), where eps is  approximately  the  square
    root  of  the  relative  machine  precision.   if   f   is a unimodal
    function and the computed values of   f   are  always  unimodal  when
    separated  by  at least  eps*abs(x)+(tol/3), then  fmin  approximates
    the abcissa of the global minimum of  f  on the interval  ax,bx  with
    an error less than  3*eps*abs(fmin)+tol.  if   f   is  not  unimodal,
    then fmin may approximate a local, but perhaps non-global, minimum to
    the same accuracy.
        This function subprogram is a slightly modified  version  of  the
    Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
    Minimization without Derivatives, Prentice-Hall, Inc. (1973).
*/

double use_brent(double ax, double bx, double (*f)(double, void *),
		  void *info, double tol, double &fx)
{
    /* added argument: double &fx, a reference for fx - HC 03/11/2014 */

    /*  c is the squared inverse of the golden ratio */
    const double c = (3. - sqrt(5.)) * .5;

    /* Local variables */
    double a, b, d, e, p, q, r, u, v, w, x;
    double t2, fu, fv, fw, xm, eps, tol1, tol3;
    /* double t2, fu, fv, fw, fx, xm, eps, tol1, tol3; - HC 03/11/2014 */

/*  eps is approximately the square root of the relative machine precision. */
    eps = DBL_EPSILON;
    tol1 = eps + 1.;/* the smallest 1.000... > 1 */
    eps = sqrt(eps);

    a = ax;
    b = bx;
    v = a + c * (b - a);
    w = v;
    x = v;

    d = 0.;/* -Wall */
    e = 0.;
    fx = (*f)(x, info);
    fv = fx;
    fw = fx;
    tol3 = tol / 3.;

/*  main loop starts here ----------------------------------- */

    for(;;) {
	xm = (a + b) * .5;
	tol1 = eps * fabs(x) + tol3;
	t2 = tol1 * 2.;

	/* check stopping criterion */

	if (fabs(x - xm) <= t2 - (b - a) * .5) break;
	p = 0.;
	q = 0.;
	r = 0.;
	if (fabs(e) > tol1) { /* fit parabola */

	    r = (x - w) * (fx - fv);
	    q = (x - v) * (fx - fw);
	    p = (x - v) * q - (x - w) * r;
	    q = (q - r) * 2.;
	    if (q > 0.) p = -p; else q = -q;
	    r = e;
	    e = d;
	}

	if (fabs(p) >= fabs(q * .5 * r) ||
	    p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

	    if (x < xm) e = b - x; else e = a - x;
	    d = c * e;
	}
	else { /* a parabolic-interpolation step */

	    d = p / q;
	    u = x + d;

	    /* f must not be evaluated too close to ax or bx */

	    if (u - a < t2 || b - u < t2) {
		d = tol1;
		if (x >= xm) d = -d;
	    }
	}

	/* f must not be evaluated too close to x */

	if (fabs(d) >= tol1)
	    u = x + d;
	else if (d > 0.)
	    u = x + tol1;
	else
	    u = x - tol1;

	fu = (*f)(u, info);

	/*  update  a, b, v, w, and x */

	if (fu <= fx) {
	    if (u < x) b = x; else a = x;
	    v = w;    w = x;   x = u;
	    fv = fw; fw = fx; fx = fu;
	} else {
	    if (u < x) a = u; else b = u;
	    if (fu <= fw || w == x) {
		v = w; fv = fw;
		w = u; fw = fu;
	    } else if (fu <= fv || v == x || v == w) {
		v = u; fv = fu;
	    }
	}
    }
    /* end of main loop */

    return x;
}
/* end of use_brent */

// For calculation of UtXY: mat(n, (p+1)*(p+2)/2); 0 <= i <= j <= p
size_t index2 (const size_t p, const size_t i, const size_t j) {
	return (2*p+3-i)*i/2+j-i;
}

// For calculation of XPY: vec((p+1)*(p+2)*(p+3)/6); 0 <= i <= j <= k <= p
// last term: index2(p-i, j-i, k-i)
size_t index3 (const size_t p, const size_t i, const size_t j, const size_t k) {
	return ((p+1)*(p+2)*(p+3)-(p+1-i)*(p+2-i)*(p+3-i))/6+(2*p+3-i-j)*(j-i)/2+k-j;
}

class Params
{
public:
	mat UtXY;
	vec eigval;
	const size_t p;
	const char method;
	const char dispersion;
};

double loglikelihood (mat &UtXY, vec &eval, const size_t p, const char method, const char dispersion)
{
	double logdetSigma = -sum(log(eval));
	size_t n = eval.n_elem;
	size_t ind2, ind3;
	double a, b, c, d, l;
	vec temp((p+1)*(p+2)*(p+3)/6);
	for(size_t j=0; j<p+1; ++j) {
		for(size_t k=j; k<p+1; ++k) {
			ind3 = index3(p, 0, j, k);
			ind2 = index2(p, j, k);
			temp[ind3] = dot(eval, UtXY.col(ind2));
		}
	}
	for(size_t i=1; i<p+1; ++i) {
		d = temp[index3(p, i-1, i-1, i-1)];
		for(size_t j=i; j<p+1; ++j) {
			c = temp[index3(p, i-1, i-1, j)];
			for(size_t k=j; k<p+1; ++k) {
				ind3 = index3(p, i, j, k);
				a = temp[index3(p, i-1, j, k)];
				b = temp[index3(p, i-1, i-1, k)];
				temp[ind3] = a-b*c/d;
			}
		}
	}
	double yPy = temp[index3(p, p, p, p)];
	if(method!='R') { //ML
		if(dispersion=='Y') { //overdispersion parameter profiled
			l = logdetSigma + (double)n*log(yPy);
		} else { //overdispersion parameter = 1
			l = logdetSigma + yPy;
		}
	} else { //REML
		double logdetXSigmaiX = 0.0;
		for(size_t i=0; i<p; ++i) {
			logdetXSigmaiX += log(temp[index3(p, i, i, i)]);
		}
		if(dispersion=='Y') { //overdispersion parameter profiled
			l = logdetSigma + logdetXSigmaiX + (double)(n-p)*log(yPy);
		} else { //overdispersion parameter = 1
			l = logdetSigma + logdetXSigmaiX + yPy;
		}
	}
	return l;
}

double Loglikelihood (double x, void *params)
{
	Params *par = (Params *) params;	
	vec eigval = par->eigval;
	eigval = 1.0 / (1.0 + eigval * exp(x));
	double l = loglikelihood(par->UtXY, eigval, par->p, par->method, par->dispersion);
	return l;
}

class Params2
{
public:
	vec Y;
        mat X;
        vec W;
        Rcpp::List Phi;
	const char method;
	const char dispersion;
        mat U;
        vec eval;
        mat UtX;
        vec UtY;
        mat cov;
        const vec tau;
        const uvec fixtau;
};

double loglikelihood2 (vec &eval, mat &UtX, vec &UtY, mat &cov, const char method, const char dispersion)
{
        double l, logdetSigma = -sum(log(eval));
	size_t n = eval.n_elem;
	size_t p = UtX.n_cols;
	mat XtSigmaiX = UtX.t() * diagmat(eval) * UtX;
	vec XtSigmaiY = UtX.t() * diagmat(eval) * UtY;
	mat U2;
	vec eval2;
	eig_sym(eval2, U2, XtSigmaiX, "dc");
	cov = U2 * diagmat(1.0 / eval2) * U2.t();
	double yPy = sum(UtY % eval % UtY) - as_scalar(XtSigmaiY.t() * cov * XtSigmaiY);
	if(method!='R') { //ML
		if(dispersion=='Y') { //overdispersion parameter profiled
			l = logdetSigma + (double)n*log(yPy);
		} else { //overdispersion parameter = 1
			l = logdetSigma + yPy;
		}
	} else { //REML
	        double logdetXSigmaiX = sum(log(eval2));
		if(dispersion=='Y') { //overdispersion parameter profiled
			l = logdetSigma + logdetXSigmaiX + (double)(n-p)*log(yPy);
		} else { //overdispersion parameter = 1
			l = logdetSigma + logdetXSigmaiX + yPy;
		}
	}
	return l;
}

double Loglikelihood2 (void *params)
{
	Params2 *par = (Params2 *) params;
	vec tau = par->tau;
	size_t i, q = tau.n_elem;
	mat Sigma = diagmat(1.0 / par->W);
	for(i=0; i<q; ++i) {
		Sigma = Sigma + tau[i] * as<mat>(par->Phi[i]);
	}
	eig_sym(par->eval, par->U, Sigma, "dc");
	par->eval = 1.0 / par->eval;
	par->UtX = par->U.t() * par->X;
	par->UtY = par->U.t() * par->Y;
	double l = loglikelihood2(par->eval, par->UtX, par->UtY, par->cov, par->method, par->dispersion);
	return l;
}

double Loglikelihood2 (int q2, double *x, void *params)
{
	Params2 *par = (Params2 *) params;
	vec tau = par->tau;
	int i, q = tau.n_elem;
	const uvec idxtau = find(par->fixtau == 0);
	for(i=0; i<q2; ++i) {
		if(x[i]<0.0) {x[i] = 0.0;}
	        tau[idxtau[i]] = x[i];
	}
	mat Sigma = diagmat(1.0 / par->W);
	for(i=0; i<q; ++i) {
		Sigma = Sigma + tau[i] * as<mat>(par->Phi[i]);
	}
	eig_sym(par->eval, par->U, Sigma, "dc");
	par->eval = 1.0 / par->eval;
	par->UtX = par->U.t() * par->X;
	par->UtY = par->U.t() * par->Y;
	double l = loglikelihood2(par->eval, par->UtX, par->UtY, par->cov, par->method, par->dispersion);
	return l;
}

extern "C" 
{
  SEXP fitglmm_brent(SEXP Y_in, SEXP X_in, SEXP Phi_in, SEXP sqrtW_in, SEXP method_in, SEXP dispersion_in, SEXP tau_in, SEXP fixtau_in, SEXP tol_in, SEXP taumin_in, SEXP taumax_in, SEXP tauregion_in)
{
	try {
	        mat X = as<mat>(X_in);
		size_t n = X.n_rows, p = X.n_cols;
		Rcpp::NumericVector Y_r(Y_in);
		Rcpp::NumericVector sqrtW_r(sqrtW_in);
		mat Phi = as<mat>(Phi_in);
		arma::vec Y(Y_r.begin(), Y_r.size(), false);
		arma::vec sqrtW(sqrtW_r.begin(), sqrtW_r.size(), false);
		const char method = Rcpp::as<char>(method_in);
		const char dispersion = Rcpp::as<char>(dispersion_in);
		double tau = Rcpp::as<double>(tau_in);
		const size_t fixtau = Rcpp::as<size_t>(fixtau_in);
		const double tol = Rcpp::as<double>(tol_in);
		const double taumin = Rcpp::as<double>(taumin_in);
		const double taumax = Rcpp::as<double>(taumax_in);
		const size_t tauregion = Rcpp::as<size_t>(tauregion_in);
		size_t ind2, i, j;
		double taua, taub, tautmp, l, ltmp;

		Y %= sqrtW;
		X.each_col() %= sqrtW;
		Phi.each_col() %= sqrtW;
		Phi.each_row() %= sqrtW.t();
		mat U;
		vec eval;
		eig_sym(eval, U, Phi, "dc");
		double tol2 = tol*tol;
		for(i=0; i<n; ++i) {
			if(eval[i] < tol2) {
				eval[i] = 0.0;
			}
		}
		mat UtX = U.t() * X;
		vec UtY = U.t() * Y;
		if(fixtau == 0) {
		        mat UtXY(n, (p+1)*(p+2)/2);
			for(i=0; i<p; ++i) {
			        for(j=i; j<p+1; ++j) {
				        ind2 = index2(p, i, j);
					if(j<p) {
					        UtXY.col(ind2) = UtX.col(i) % UtX.col(j);
					} else {
					        UtXY.col(ind2) = UtX.col(i) % UtY;
					}
				}
			}
			ind2 = index2(p, p, p);
			UtXY.col(ind2) = UtY % UtY;
			Params parameters = {UtXY, eval, p, method, dispersion};

			for (i=0; i<tauregion; ++i) {
			        taua = log(taumin) + log(taumax/taumin) * (double)i / (double)tauregion;
				taub = log(taumin) + log(taumax/taumin) * ((double)i + 1.0) / (double)tauregion;
				tautmp = use_brent(taua, taub, &Loglikelihood, &parameters, tol, ltmp);
				if(i==0) {
				        tau = tautmp;
					l = ltmp;
				} else {
				        if(ltmp<l) {
					        tau = tautmp;
						l = ltmp;
					}
				}
			}
			tau = exp(tau);
		}
		eval = 1.0 / (1.0 + eval * tau);
		mat cov = inv_sympd(UtX.t() * diagmat(eval) * UtX);
		vec alpha = cov * UtX.t() * diagmat(eval) * UtY;
		vec eta = U * (UtY - diagmat(eval) * (UtY - UtX * alpha)) / sqrtW;
		U.each_col() %= sqrtW;
		return List::create(Named("tau") = tau, Named("U") = U, Named("eval") = eval, Named("UtX") = UtX, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

typedef double optimfn(int, double *, void *);
void nmmin(int n, double *Bvec, double *X, double *Fmin, optimfn fn,
           int *fail, double abstol, double intol, void *ex,
           double alpha, double bet, double gamm, int trace,
           int *fncount, int maxit);

  SEXP fitglmm_nm(SEXP Y_in, SEXP X_in, SEXP q_in, SEXP Phi_in, SEXP W_in, SEXP method_in, SEXP dispersion_in, SEXP tau_in, SEXP fixtau_in, SEXP maxiter_in, SEXP tol_in)
{
	try {
		Rcpp::NumericMatrix X_r(X_in);
		size_t n = X_r.nrow(), p = X_r.ncol();
		Rcpp::NumericVector Y_r(Y_in);
		Rcpp::NumericVector W_r(W_in);
		arma::mat X(X_r.begin(), n, p, false);
		arma::vec Y(Y_r.begin(), Y_r.size(), false);
		arma::vec W(W_r.begin(), W_r.size(), false);
		const size_t q = Rcpp::as<size_t>(q_in);
		const Rcpp::List Phi(Phi_in);
		const char method = Rcpp::as<char>(method_in);
		const char dispersion = Rcpp::as<char>(dispersion_in);
		vec tau = as<vec>(tau_in);
		const uvec fixtau = as<uvec>(fixtau_in);
		const size_t q2 = sum(fixtau == 0);
		const size_t maxiter = Rcpp::as<size_t>(maxiter_in);
		const double tol = Rcpp::as<double>(tol_in);
		mat U, UtX, cov;
		vec eval, UtY;
		Params2 parameters = {Y, X, W, Phi, method, dispersion, U, eval, UtX, UtY, cov, tau, fixtau};
		double Fmin = 0.0;

		if(q2 > 0) {
		        const uvec idxtau = find(fixtau == 0);
		        int fail = 0, fncount = 0;
			//double xin[q2], x[q2];
			double *xin = new double[q2];
			double *x = new double[q2];
			for(size_t i=0; i<q2; ++i) {
			        xin[i] = tau[idxtau[i]];
				x[i] = tau[idxtau[i]];
			}
			nmmin(q, xin, x, &Fmin, Loglikelihood2, &fail, -INFINITY, tol, &parameters, 1.0, 0.5, 2.0, 0, &fncount, maxiter);
			if(fail != 0) {
			        Rcout << "Warning: convergence may not be reached!\n";
			}
			for(size_t i=0; i<q2; ++i) {
			        tau[idxtau[i]] = x[i];
			}
			delete [] xin;
			delete [] x;
		} else {
		        Fmin = Loglikelihood2(&parameters);
		}
		vec alpha = parameters.cov * parameters.UtX.t() * diagmat(parameters.eval) * parameters.UtY;
		vec eta = Y - (parameters.U * diagmat(parameters.eval) * (parameters.UtY - parameters.UtX * alpha)) / W;
		return List::create(Named("tau") = tau, Named("U") = parameters.U, Named("eval") = parameters.eval, Named("UtX") = parameters.UtX, Named("cov") = parameters.cov, Named("alpha") = alpha, Named("eta") = eta);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

  //  SEXP fitglmm_ai(SEXP Y_in, SEXP X_in, SEXP q_in, SEXP Phi_in, SEXP ng_in, SEXP group_in, SEXP W_in, SEXP tau_in, SEXP fixtau_in, SEXP tol_in)
  SEXP fitglmm_ai(SEXP Y_in, SEXP X_in, SEXP q_in, SEXP Phi_in, SEXP ng_in, SEXP group_in, SEXP W_in, SEXP tau_in, SEXP fixtau_in)
{
	try {
		Rcpp::NumericMatrix X_r(X_in);
		size_t n = X_r.nrow(), p = X_r.ncol();
		Rcpp::NumericVector Y_r(Y_in);
		Rcpp::NumericVector W_r(W_in);
		arma::mat X(X_r.begin(), n, p, false);
		arma::vec Y(Y_r.begin(), Y_r.size(), false);
		arma::vec W(W_r.begin(), W_r.size(), false);
		const size_t q = Rcpp::as<size_t>(q_in);
		const size_t ng = Rcpp::as<size_t>(ng_in);
		const Rcpp::List Phi(Phi_in);
		const Rcpp::List group(group_in);
		vec tau = as<vec>(tau_in);
		const uvec fixtau = as<uvec>(fixtau_in);
		const size_t q2 = sum(fixtau == 0);
		//const double tol = Rcpp::as<double>(tol_in);
		//uvec ZERO = (tau < tol);
		mat cov(p, p);
		vec alpha(p), eta(n);
		vec diagP = zeros<vec>(n);
		for(size_t i=1; i<=ng; ++i) {
		        uvec group_idx = as<uvec>(group[i-1]) - 1;
			diagP.elem( group_idx ) = tau[i-1] / W.elem( group_idx );
		}
		mat P = diagmat(diagP);
		for(size_t i=1; i<=q; ++i) {
			P = P + tau[i+ng-1] * as<mat>(Phi[i-1]);
		}
		mat Sigma_i = inv_sympd(P);
		mat Sigma_iX = Sigma_i * X;
		cov = inv_sympd(X.t() * Sigma_iX);
		P = Sigma_i - Sigma_iX * cov * Sigma_iX.t();
		alpha = cov * Sigma_iX.t() * Y;
		eta = Y - diagP % (Sigma_i * (Y - X * alpha));
       		if(q2 > 0) {
		        const uvec idxtau = find(fixtau == 0);
		        mat AI(q2, q2);
		        vec PY = P * Y;
			vec score(q2), PAPY;
			vec APY = PY / W;
			diagP = diagvec(P) / W;
			for(size_t i=0; i<q2; ++i) {
			        if(idxtau[i] < ng) {
				        uvec group_idx = as<uvec>(group[idxtau[i]]) - 1;
					score[i] = dot(APY.elem( group_idx ), PY.elem( group_idx )) - sum(diagP.elem( group_idx ));
					for(size_t j=0; j<=i; ++j) {
					        uvec group_idx2 = as<uvec>(group[idxtau[j]]) - 1;
						AI(i,j) = dot(APY.elem( group_idx ), P.submat( group_idx, group_idx2) * APY.elem( group_idx2 ));
						if(j!=i) {AI(j,i) = AI(i,j);}
					}
				} else {
				        PAPY = P * as<mat>(Phi[idxtau[i]-ng]) * PY;
					score[i] = dot(Y, PAPY) - accu(P % as<mat>(Phi[idxtau[i]-ng]));
					for(size_t j=0; j<=i; ++j) {
					        if(idxtau[j] < ng) {
						        uvec group_idx = as<uvec>(group[idxtau[j]]) - 1;
						        AI(i,j) = dot(APY.elem( group_idx ), PAPY.elem( group_idx ));
							AI(j,i) = AI(i,j);
						} else {
						        AI(i,j) = dot(PY, as<mat>(Phi[idxtau[j]-ng]) * PAPY);
							if(j!=i) {AI(j,i) = AI(i,j);}
						}
					}
				}
			}
			vec Dtau = solve(AI, score);
			//vec tau0 = tau;
			//tau.elem( idxtau ) = tau0.elem( idxtau ) + Dtau;
			//tau.elem( find(ZERO % (tau < tol)) ).zeros();
			//double step = 1.0;
			//while(any(tau < 0.0)) {
			//        step *= 0.5;
			//	tau.elem( idxtau ) = tau0.elem( idxtau ) + step * Dtau;
			//	tau.elem( find(ZERO % (tau < tol)) ).zeros();
			//}
			//tau.elem( find(tau < tol) ).zeros();
			return List::create(Named("Dtau") = Dtau, Named("P") = P, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
		} else {
			return List::create(Named("Dtau") = R_NilValue, Named("P") = P, Named("cov") = cov, Named("alpha") = alpha, Named("eta") = eta);
		}
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

  SEXP glmm_score_text(SEXP res_in, SEXP P_in, SEXP infile_in, SEXP outfile_in, SEXP tol_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP p_in, SEXP infile_nrow_skip_in, SEXP infile_sep_in, SEXP infile_na_in, SEXP infile_ncol_skip_in, SEXP infile_ncol_print_in, SEXP infile_header_print_in, SEXP nperbatch_in, SEXP select_in) {
	try {
		Rcpp::NumericMatrix P_r(P_in);
		Rcpp::NumericVector res_r(res_in);
		size_t n_P = P_r.nrow(), k_P = P_r.ncol();
		arma::mat P(P_r.begin(), n_P, k_P, false);
		arma::vec res(res_r.begin(), res_r.size(), false);
		const double tol = Rcpp::as<double>(tol_in);
		const char center = Rcpp::as<char>(center_in);
		const double minmaf = Rcpp::as<double>(minmaf_in);
		const double maxmaf = Rcpp::as<double>(maxmaf_in);
		const double missrate = Rcpp::as<double>(missrate_in);
		const char miss_method = Rcpp::as<char>(miss_method_in);
		string infile = Rcpp::as<string>(infile_in);
		string outfile = Rcpp::as<string>(outfile_in);
		const size_t p = Rcpp::as<size_t>(p_in);
		const size_t infile_nrow_skip = Rcpp::as<size_t>(infile_nrow_skip_in);
		string infile_sep = Rcpp::as<string>(infile_sep_in);
		string infile_na = Rcpp::as<string>(infile_na_in);
		const int infile_ncol_skip = Rcpp::as<int>(infile_ncol_skip_in);
		Rcpp::IntegerVector infile_ncol_print(infile_ncol_print_in);
		Rcpp::CharacterVector infile_header_print(infile_header_print_in);
		const size_t npb = Rcpp::as<size_t>(nperbatch_in);
		Rcpp::IntegerVector select(select_in);
		string snp;
		char *cp;
		size_t n = res.n_elem, ns = select.size();
		vec g(n);
		uvec gmiss(n), snp_skip = zeros<uvec>(npb);
		mat G(n, npb);
		string *tmpout = new string[npb];
		double gmean, geno, gmax, gmin, num, denom, pval;
		//size_t nmiss, infile_ncol_print_idx, npbidx = 0;
		size_t nmiss, npbidx = 0;
		int infile_ncol_print_idx;
		clock_t time0;
		double compute_time = 0.0;
		ofstream writefile (outfile.c_str(), ofstream::out);
		if (!writefile) {Rcout << "Error writing file: " << outfile << "\n"; return R_NilValue;}
		if (infile_ncol_print[0] != 0 && strcmp(infile_header_print[0], infile_na.c_str()) != 0) {
		      for(int i=0; i<infile_header_print.size(); ++i) {
			    writefile << infile_header_print[i] << "\t";
		      }
		}
		writefile << "N" << "\t" << "AF" << "\t" << "SCORE" << "\t" << "VAR" << "\t" << "PVAL" << "\n";
		size_t tmppos = infile.rfind('.');
		string ext = infile.substr(tmppos);
		if (strcmp(ext.c_str(), ".bz2")==0) { // bzip2 infile
		      FILE* fp = fopen(infile.c_str(), "rb");
		      int bzerror;
		      BZFILE* bfp = BZ2_bzReadOpen(&bzerror, fp, 0, 0, NULL, 0);;
		      int buffer_length = 100 * (int(ns) + infile_ncol_skip);
		      if(bzerror != BZ_OK) {
			    BZ2_bzReadClose(&bzerror, bfp);
			    fclose(fp);
			    Rcout << "File " << infile << " appears not to be compressed by bzip2\n";
			    return R_NilValue;
		      }
		      for(size_t i=0; i<p; ++i) {
			    int nread = 0;
			    char *buffer = new char[buffer_length];
			    while (bzerror == BZ_OK) {
			          int nn = BZ2_bzRead(&bzerror, bfp, buffer + nread, 1);
				  if (bzerror == BZ_STREAM_END) {
				        char *unused, *next_unused = NULL;
					int nUnused;
					BZ2_bzReadGetUnused(&bzerror, bfp, (void**) &unused, &nUnused);
					if (bzerror == BZ_OK) {
					      if (nUnused > 0) {
						    next_unused = (char*) malloc(nUnused);
						    if (!next_unused) {
						          BZ2_bzReadClose(&bzerror, bfp);
							  fclose(fp);
						          Rcout << "Allocation of overflow buffer for bzfile failed\n";
							  return R_NilValue;
						    }
						    memcpy(next_unused, unused, nUnused);
					      }
					      if (nUnused > 0 || !feof(fp)) {
						    BZ2_bzReadClose(&bzerror, bfp);
						    bfp = BZ2_bzReadOpen(&bzerror, fp, 0, 0, next_unused, nUnused);
						    if(bzerror != BZ_OK) {
						          BZ2_bzReadClose(&bzerror, bfp);
							  fclose(fp);
						          Rcout << "File " << infile << " has trailing content that appears not to be compressed by bzip2\n";
							  return R_NilValue;
						    }
					      }
					      if (next_unused) free(next_unused);
					}
				  } else if (bzerror != BZ_OK) {
				        nread += nn;
					break;
				  }
				  nread += nn;
				  if(nread > 0 && buffer[nread-1] == '\n') {
				        buffer[nread-1] = '\0';
					break;
				  }
			    }
			    if(i<infile_nrow_skip) {
			          delete [] buffer;
				  continue;
			    }
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    stringstream writeout;
			    cp=strtok (buffer, infile_sep.c_str());
			    if(infile_ncol_skip == 0) {
			          if(select[0] > 0) {
				        if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[0]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[0]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }
				  for (size_t j=1; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
				 	      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    } else {
			          infile_ncol_print_idx = 0;
				  if (infile_ncol_print[0] == 1) {
				        snp=cp;
					writeout << snp << "\t";
					infile_ncol_print_idx++;
				  }
				  if (infile_ncol_skip > 1) {
				        for(int k=1; k<infile_ncol_skip; ++k) {
				              cp=strtok (NULL, infile_sep.c_str());
					      if (infile_ncol_print_idx<infile_ncol_print.size()) {
					            if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						          snp=cp;
							  writeout << snp << "\t";
							  infile_ncol_print_idx++;
						    }
					      }
					}
				  }
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    }
			    gmean/=(double)(n-nmiss);
			    for (size_t j=0; j<n; ++j) {
			          if (gmiss[j]==1) {
					g[j] = gmean;
					if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
				  }
				  if (center=='c') {
					g[j] -= gmean;
				  }
			    }
			    gmean/=2.0; // convert mean to allele freq
			    writeout << n-nmiss << "\t" << gmean << "\t";
			    tmpout[npbidx] = writeout.str();
			    writeout.clear();
			    if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			          snp_skip[npbidx] = 1;
			    } else {
			          G.col(npbidx) = g;
			    }
			    npbidx++;
			    if((i+1 == p) || (npbidx == npb)) {
			          uvec snp_idx = find(snp_skip == 0);
				  time0 = clock();
				  vec allnum = G.cols(snp_idx).t() * res;
				  mat alldenom = G.cols(snp_idx).t() * P * G.cols(snp_idx);
				  compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
				  for(size_t j=0; j<npbidx; ++j) {
				        if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
					      writefile << tmpout[j] << "NA\tNA\tNA\n";
					} else {
					      uvec snp_idx2 = find(snp_idx == j);
					      num = allnum[snp_idx2[0]];
					      denom = alldenom(snp_idx2[0], snp_idx2[0]);
					      if(denom < tol) { // G collinearity with X
						    writefile << tmpout[j] << "0\t0\tNA\n";
					      } else {
						    pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
						    writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					      }
					}
				  }
				  if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
			    }
			    if((i+1) % 100000 == 0) {writefile << flush;}
			    delete [] buffer;
		      }
		      if(p % 100000 != 0) {writefile << flush;}
		      writefile.close();
		      writefile.clear();
		      BZ2_bzReadClose(&bzerror, bfp);
		      fclose(fp);
		} else if (strcmp(ext.c_str(), ".gz")==0) { // gzip infile
		      gzFile fp = gzopen(infile.c_str(), "rb");
		      int buffer_length = 100 * (int(ns) + infile_ncol_skip);
		      if(!fp) {
			    Rcout << "Error: cannot open gzipped file " << infile << "\n";
			    return R_NilValue;
		      }
		      for(size_t i=0; i<p; ++i) {
			    char *buffer = new char[buffer_length];
			    gzgets(fp, buffer, buffer_length);
			    if(i<infile_nrow_skip) {
			          delete [] buffer;
				  continue;
			    }
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    stringstream writeout;
			    cp=strtok (buffer, infile_sep.c_str());
			    if(infile_ncol_skip == 0) {
			          if(select[0] > 0) {
				        if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[0]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[0]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }
				  for (size_t j=1; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
				 	      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    } else {
			          infile_ncol_print_idx = 0;
				  if (infile_ncol_print[0] == 1) {
				        snp=cp;
					writeout << snp << "\t";
					infile_ncol_print_idx++;
				  }
				  if (infile_ncol_skip > 1) {
				        for(int k=1; k<infile_ncol_skip; ++k) {
				              cp=strtok (NULL, infile_sep.c_str());
					      if (infile_ncol_print_idx<infile_ncol_print.size()) {
					            if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						          snp=cp;
							  writeout << snp << "\t";
							  infile_ncol_print_idx++;
						    }
					      }
					}
				  }
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    }
			    gmean/=(double)(n-nmiss);
			    for (size_t j=0; j<n; ++j) {
			          if (gmiss[j]==1) {
					g[j] = gmean;
					if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
				  }
				  if (center=='c') {
					g[j] -= gmean;
				  }
			    }
			    gmean/=2.0; // convert mean to allele freq
			    writeout << n-nmiss << "\t" << gmean << "\t";
			    tmpout[npbidx] = writeout.str();
			    writeout.clear();
			    if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			          snp_skip[npbidx] = 1;
			    } else {
			          G.col(npbidx) = g;
			    }
			    npbidx++;
			    if((i+1 == p) || (npbidx == npb)) {
			          uvec snp_idx = find(snp_skip == 0);
				  time0 = clock();
				  vec allnum = G.cols(snp_idx).t() * res;
				  mat alldenom = G.cols(snp_idx).t() * P * G.cols(snp_idx);
				  compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
				  for(size_t j=0; j<npbidx; ++j) {
				        if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
					      writefile << tmpout[j] << "NA\tNA\tNA\n";
					} else {
					      uvec snp_idx2 = find(snp_idx == j);
					      num = allnum[snp_idx2[0]];
					      denom = alldenom(snp_idx2[0], snp_idx2[0]);
					      if(denom < tol) { // G collinearity with X
						    writefile << tmpout[j] << "0\t0\tNA\n";
					      } else {
						    pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
						    writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					      }
					}
				  }
				  if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
			    }
			    if((i+1) % 100000 == 0) {writefile << flush;}
			    delete [] buffer;
		      }
		      if(p % 100000 != 0) {writefile << flush;}
		      writefile.close();
		      writefile.clear();
		      gzclose(fp);
		} else { // plain text infile
		      ifstream readfile (infile.c_str(), ifstream::in);
		      if (!readfile) {Rcout << "Error reading genotype file: " << infile << "\n"; return R_NilValue;}
		      string line;
		      for(size_t i=0; i<p; ++i) {
			    getline(readfile, line);
			    if(i<infile_nrow_skip) {continue;}
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    stringstream writeout;
			    cp=strtok ((char *)line.c_str(), infile_sep.c_str());
			    if(infile_ncol_skip == 0) {
			          if(select[0] > 0) {
				        if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[0]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[0]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }
				  for (size_t j=1; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
				 	      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    } else {
			          infile_ncol_print_idx = 0;
				  if (infile_ncol_print[0] == 1) {
				        snp=cp;
					writeout << snp << "\t";
					infile_ncol_print_idx++;
				  }
				  if (infile_ncol_skip > 1) {
				        for(int k=1; k<infile_ncol_skip; ++k) {
				              cp=strtok (NULL, infile_sep.c_str());
					      if (infile_ncol_print_idx<infile_ncol_print.size()) {
					            if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						          snp=cp;
							  writeout << snp << "\t";
							  infile_ncol_print_idx++;
						    }
					      }
					}
				  }
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    }
			    gmean/=(double)(n-nmiss);
			    for (size_t j=0; j<n; ++j) {
			          if (gmiss[j]==1) {
					g[j] = gmean;
					if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
				  }
				  if (center=='c') {
					g[j] -= gmean;
				  }
			    }
			    gmean/=2.0; // convert mean to allele freq
			    writeout << n-nmiss << "\t" << gmean << "\t";
			    tmpout[npbidx] = writeout.str();
			    writeout.clear();
			    if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			          snp_skip[npbidx] = 1;
			    } else {
			          G.col(npbidx) = g;
			    }
			    npbidx++;
			    if((i+1 == p) || (npbidx == npb)) {
			          uvec snp_idx = find(snp_skip == 0);
				  time0 = clock();
				  vec allnum = G.cols(snp_idx).t() * res;
				  mat alldenom = G.cols(snp_idx).t() * P * G.cols(snp_idx);
				  compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
				  for(size_t j=0; j<npbidx; ++j) {
				        if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
					      writefile << tmpout[j] << "NA\tNA\tNA\n";
					} else {
					      uvec snp_idx2 = find(snp_idx == j);
					      num = allnum[snp_idx2[0]];
					      denom = alldenom(snp_idx2[0], snp_idx2[0]);
					      if(denom < tol) { // G collinearity with X
						    writefile << tmpout[j] << "0\t0\tNA\n";
					      } else {
						    pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
						    writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					      }
					}
				  }
				  if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
			    }
			    if((i+1) % 100000 == 0) {writefile << flush;}
		      }
		      if(p % 100000 != 0) {writefile << flush;}
		      writefile.close();
		      writefile.clear();
		      readfile.close();
		      readfile.clear();
		}
		delete [] tmpout;
		return wrap(compute_time);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

  SEXP glmm_score_text_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP infile_in, SEXP outfile_in, SEXP tol_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP p_in, SEXP infile_nrow_skip_in, SEXP infile_sep_in, SEXP infile_na_in, SEXP infile_ncol_skip_in, SEXP infile_ncol_print_in, SEXP infile_header_print_in, SEXP nperbatch_in, SEXP select_in) {
	try {
		Rcpp::NumericVector res_r(res_in);
		arma::vec res(res_r.begin(), res_r.size(), false);
		arma::sp_mat Sigma_i = Rcpp::as<arma::sp_mat>(Sigma_i_in);
		arma::sp_mat Sigma_iX = Rcpp::as<arma::sp_mat>(Sigma_iX_in);
		arma::sp_mat cov = Rcpp::as<arma::sp_mat>(cov_in);
		const double tol = Rcpp::as<double>(tol_in);
		const char center = Rcpp::as<char>(center_in);
		const double minmaf = Rcpp::as<double>(minmaf_in);
		const double maxmaf = Rcpp::as<double>(maxmaf_in);
		const double missrate = Rcpp::as<double>(missrate_in);
		const char miss_method = Rcpp::as<char>(miss_method_in);
		string infile = Rcpp::as<string>(infile_in);
		string outfile = Rcpp::as<string>(outfile_in);
		const size_t p = Rcpp::as<size_t>(p_in);
		const size_t infile_nrow_skip = Rcpp::as<size_t>(infile_nrow_skip_in);
		string infile_sep = Rcpp::as<string>(infile_sep_in);
		string infile_na = Rcpp::as<string>(infile_na_in);
		const int infile_ncol_skip = Rcpp::as<int>(infile_ncol_skip_in);
		Rcpp::IntegerVector infile_ncol_print(infile_ncol_print_in);
		Rcpp::CharacterVector infile_header_print(infile_header_print_in);
		const size_t npb = Rcpp::as<size_t>(nperbatch_in);
		Rcpp::IntegerVector select(select_in);
		string snp;
		char *cp;
		size_t n = res.n_elem, ns = select.size();
		vec g(n);
		uvec gmiss(n), snp_skip = zeros<uvec>(npb);
		mat G(n, npb);
		string *tmpout = new string[npb];
		double gmean, geno, gmax, gmin, num, denom, pval;
		//size_t nmiss, infile_ncol_print_idx, npbidx = 0;
		size_t nmiss, npbidx = 0;
		int infile_ncol_print_idx;
		clock_t time0;
		double compute_time = 0.0;
		ofstream writefile (outfile.c_str(), ofstream::out);
		if (!writefile) {Rcout << "Error writing file: " << outfile << "\n"; return R_NilValue;}
		if (infile_ncol_print[0] != 0 && strcmp(infile_header_print[0], infile_na.c_str()) != 0) {
		      for(int i=0; i<infile_header_print.size(); ++i) {
			    writefile << infile_header_print[i] << "\t";
		      }
		}
		writefile << "N" << "\t" << "AF" << "\t" << "SCORE" << "\t" << "VAR" << "\t" << "PVAL" << "\n";
		size_t tmppos = infile.rfind('.');
		string ext = infile.substr(tmppos);
		if (strcmp(ext.c_str(), ".bz2")==0) { // bzip2 infile
		      FILE* fp = fopen(infile.c_str(), "rb");
		      int bzerror;
		      BZFILE* bfp = BZ2_bzReadOpen(&bzerror, fp, 0, 0, NULL, 0);;
		      int buffer_length = 100 * (int(ns) + infile_ncol_skip);
		      if(bzerror != BZ_OK) {
			    BZ2_bzReadClose(&bzerror, bfp);
			    fclose(fp);
			    Rcout << "File " << infile << " appears not to be compressed by bzip2\n";
			    return R_NilValue;
		      }
		      for(size_t i=0; i<p; ++i) {
			    int nread = 0;
			    char *buffer = new char[buffer_length];
			    while (bzerror == BZ_OK) {
			          int nn = BZ2_bzRead(&bzerror, bfp, buffer + nread, 1);
				  if (bzerror == BZ_STREAM_END) {
				        char *unused, *next_unused = NULL;
					int nUnused;
					BZ2_bzReadGetUnused(&bzerror, bfp, (void**) &unused, &nUnused);
					if (bzerror == BZ_OK) {
					      if (nUnused > 0) {
						    next_unused = (char*) malloc(nUnused);
						    if (!next_unused) {
						          BZ2_bzReadClose(&bzerror, bfp);
							  fclose(fp);
						          Rcout << "Allocation of overflow buffer for bzfile failed\n";
							  return R_NilValue;
						    }
						    memcpy(next_unused, unused, nUnused);
					      }
					      if (nUnused > 0 || !feof(fp)) {
						    BZ2_bzReadClose(&bzerror, bfp);
						    bfp = BZ2_bzReadOpen(&bzerror, fp, 0, 0, next_unused, nUnused);
						    if(bzerror != BZ_OK) {
						          BZ2_bzReadClose(&bzerror, bfp);
							  fclose(fp);
						          Rcout << "File " << infile << " has trailing content that appears not to be compressed by bzip2\n";
							  return R_NilValue;
						    }
					      }
					      if (next_unused) free(next_unused);
					}
				  } else if (bzerror != BZ_OK) {
				        nread += nn;
					break;
				  }
				  nread += nn;
				  if(nread > 0 && buffer[nread-1] == '\n') {
				        buffer[nread-1] = '\0';
					break;
				  }
			    }
			    if(i<infile_nrow_skip) {
			          delete [] buffer;
				  continue;
			    }
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    stringstream writeout;
			    cp=strtok (buffer, infile_sep.c_str());
			    if(infile_ncol_skip == 0) {
			          if(select[0] > 0) {
				        if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[0]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[0]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }
				  for (size_t j=1; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
				 	      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    } else {
			          infile_ncol_print_idx = 0;
				  if (infile_ncol_print[0] == 1) {
				        snp=cp;
					writeout << snp << "\t";
					infile_ncol_print_idx++;
				  }
				  if (infile_ncol_skip > 1) {
				        for(int k=1; k<infile_ncol_skip; ++k) {
				              cp=strtok (NULL, infile_sep.c_str());
					      if (infile_ncol_print_idx<infile_ncol_print.size()) {
					            if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						          snp=cp;
							  writeout << snp << "\t";
							  infile_ncol_print_idx++;
						    }
					      }
					}
				  }
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    }
			    gmean/=(double)(n-nmiss);
			    for (size_t j=0; j<n; ++j) {
			          if (gmiss[j]==1) {
					g[j] = gmean;
					if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
				  }
				  if (center=='c') {
					g[j] -= gmean;
				  }
			    }
			    gmean/=2.0; // convert mean to allele freq
			    writeout << n-nmiss << "\t" << gmean << "\t";
			    tmpout[npbidx] = writeout.str();
			    writeout.clear();
			    if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			          snp_skip[npbidx] = 1;
			    } else {
			          G.col(npbidx) = g;
			    }
			    npbidx++;
			    if((i+1 == p) || (npbidx == npb)) {
			          uvec snp_idx = find(snp_skip == 0);
				  time0 = clock();
				  vec allnum = G.cols(snp_idx).t() * res;
				  sp_mat Gsp(G.cols(snp_idx));
				  sp_mat XSigma_iG = Sigma_iX.t() * Gsp;
				  sp_mat alldenom = Gsp.t() * Sigma_i * Gsp - XSigma_iG.t() * cov * XSigma_iG;
				  compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
				  for(size_t j=0; j<npbidx; ++j) {
				        if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
					      writefile << tmpout[j] << "NA\tNA\tNA\n";
					} else {
					      uvec snp_idx2 = find(snp_idx == j);
					      num = allnum[snp_idx2[0]];
					      denom = alldenom(snp_idx2[0], snp_idx2[0]);
					      if(denom < tol) { // G collinearity with X
						    writefile << tmpout[j] << "0\t0\tNA\n";
					      } else {
						    pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
						    writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					      }
					}
				  }
				  if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
			    }
			    if((i+1) % 100000 == 0) {writefile << flush;}
			    delete [] buffer;
		      }
		      if(p % 100000 != 0) {writefile << flush;}
		      writefile.close();
		      writefile.clear();
		      BZ2_bzReadClose(&bzerror, bfp);
		      fclose(fp);
		} else if (strcmp(ext.c_str(), ".gz")==0) { // gzip infile
		      gzFile fp = gzopen(infile.c_str(), "rb");
		      int buffer_length = 100 * (int(ns) + infile_ncol_skip);
		      if(!fp) {
			    Rcout << "Error: cannot open gzipped file " << infile << "\n";
			    return R_NilValue;
		      }
		      for(size_t i=0; i<p; ++i) {
			    char *buffer = new char[buffer_length];
			    gzgets(fp, buffer, buffer_length);
			    if(i<infile_nrow_skip) {
			          delete [] buffer;
				  continue;
			    }
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    stringstream writeout;
			    cp=strtok (buffer, infile_sep.c_str());
			    if(infile_ncol_skip == 0) {
			          if(select[0] > 0) {
				        if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[0]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[0]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }
				  for (size_t j=1; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
				 	      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    } else {
			          infile_ncol_print_idx = 0;
				  if (infile_ncol_print[0] == 1) {
				        snp=cp;
					writeout << snp << "\t";
					infile_ncol_print_idx++;
				  }
				  if (infile_ncol_skip > 1) {
				        for(int k=1; k<infile_ncol_skip; ++k) {
				              cp=strtok (NULL, infile_sep.c_str());
					      if (infile_ncol_print_idx<infile_ncol_print.size()) {
					            if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						          snp=cp;
							  writeout << snp << "\t";
							  infile_ncol_print_idx++;
						    }
					      }
					}
				  }
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    }
			    gmean/=(double)(n-nmiss);
			    for (size_t j=0; j<n; ++j) {
			          if (gmiss[j]==1) {
					g[j] = gmean;
					if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
				  }
				  if (center=='c') {
					g[j] -= gmean;
				  }
			    }
			    gmean/=2.0; // convert mean to allele freq
			    writeout << n-nmiss << "\t" << gmean << "\t";
			    tmpout[npbidx] = writeout.str();
			    writeout.clear();
			    if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			          snp_skip[npbidx] = 1;
			    } else {
			          G.col(npbidx) = g;
			    }
			    npbidx++;
			    if((i+1 == p) || (npbidx == npb)) {
			          uvec snp_idx = find(snp_skip == 0);
				  time0 = clock();
				  vec allnum = G.cols(snp_idx).t() * res;
				  sp_mat Gsp(G.cols(snp_idx));
				  sp_mat XSigma_iG = Sigma_iX.t() * Gsp;
				  sp_mat alldenom = Gsp.t() * Sigma_i * Gsp - XSigma_iG.t() * cov * XSigma_iG;
				  compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
				  for(size_t j=0; j<npbidx; ++j) {
				        if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
					      writefile << tmpout[j] << "NA\tNA\tNA\n";
					} else {
					      uvec snp_idx2 = find(snp_idx == j);
					      num = allnum[snp_idx2[0]];
					      denom = alldenom(snp_idx2[0], snp_idx2[0]);
					      if(denom < tol) { // G collinearity with X
						    writefile << tmpout[j] << "0\t0\tNA\n";
					      } else {
						    pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
						    writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					      }
					}
				  }
				  if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
			    }
			    if((i+1) % 100000 == 0) {writefile << flush;}
			    delete [] buffer;
		      }
		      if(p % 100000 != 0) {writefile << flush;}
		      writefile.close();
		      writefile.clear();
		      gzclose(fp);
		} else { // plain text infile
		      ifstream readfile (infile.c_str(), ifstream::in);
		      if (!readfile) {Rcout << "Error reading genotype file: " << infile << "\n"; return R_NilValue;}
		      string line;
		      for(size_t i=0; i<p; ++i) {
			    getline(readfile, line);
			    if(i<infile_nrow_skip) {continue;}
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    stringstream writeout;
			    cp=strtok ((char *)line.c_str(), infile_sep.c_str());
			    if(infile_ncol_skip == 0) {
			          if(select[0] > 0) {
				        if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[0]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[0]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }
				  for (size_t j=1; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
				 	      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    } else {
			          infile_ncol_print_idx = 0;
				  if (infile_ncol_print[0] == 1) {
				        snp=cp;
					writeout << snp << "\t";
					infile_ncol_print_idx++;
				  }
				  if (infile_ncol_skip > 1) {
				        for(int k=1; k<infile_ncol_skip; ++k) {
				              cp=strtok (NULL, infile_sep.c_str());
					      if (infile_ncol_print_idx<infile_ncol_print.size()) {
					            if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						          snp=cp;
							  writeout << snp << "\t";
							  infile_ncol_print_idx++;
						    }
					      }
					}
				  }
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
			    }
			    gmean/=(double)(n-nmiss);
			    for (size_t j=0; j<n; ++j) {
			          if (gmiss[j]==1) {
					g[j] = gmean;
					if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
				  }
				  if (center=='c') {
					g[j] -= gmean;
				  }
			    }
			    gmean/=2.0; // convert mean to allele freq
			    writeout << n-nmiss << "\t" << gmean << "\t";
			    tmpout[npbidx] = writeout.str();
			    writeout.clear();
			    if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			          snp_skip[npbidx] = 1;
			    } else {
			          G.col(npbidx) = g;
			    }
			    npbidx++;
			    if((i+1 == p) || (npbidx == npb)) {
			          uvec snp_idx = find(snp_skip == 0);
				  time0 = clock();
				  vec allnum = G.cols(snp_idx).t() * res;
				  sp_mat Gsp(G.cols(snp_idx));
				  sp_mat XSigma_iG = Sigma_iX.t() * Gsp;
				  sp_mat alldenom = Gsp.t() * Sigma_i * Gsp - XSigma_iG.t() * cov * XSigma_iG;
				  compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
				  for(size_t j=0; j<npbidx; ++j) {
				        if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
					      writefile << tmpout[j] << "NA\tNA\tNA\n";
					} else {
					      uvec snp_idx2 = find(snp_idx == j);
					      num = allnum[snp_idx2[0]];
					      denom = alldenom(snp_idx2[0], snp_idx2[0]);
					      if(denom < tol) { // G collinearity with X
						    writefile << tmpout[j] << "0\t0\tNA\n";
					      } else {
						    pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
						    writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					      }
					}
				  }
				  if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
			    }
			    if((i+1) % 100000 == 0) {writefile << flush;}
		      }
		      if(p % 100000 != 0) {writefile << flush;}
		      writefile.close();
		      writefile.clear();
		      readfile.close();
		      readfile.clear();
		}
		delete [] tmpout;
		return wrap(compute_time);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

  SEXP glmm_score_bed(SEXP res_in, SEXP P_in, SEXP bimfile_in, SEXP bedfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in) {
	try {
		Rcpp::NumericMatrix P_r(P_in);
		Rcpp::NumericVector res_r(res_in);
		size_t n_P = P_r.nrow(), k_P = P_r.ncol();
		arma::mat P(P_r.begin(), n_P, k_P, false);
		arma::vec res(res_r.begin(), res_r.size(), false);
	  const char center = Rcpp::as<char>(center_in);
		const double minmaf = Rcpp::as<double>(minmaf_in);
		const double maxmaf = Rcpp::as<double>(maxmaf_in);
		const double missrate = Rcpp::as<double>(missrate_in);
		const char miss_method = Rcpp::as<char>(miss_method_in);
		string bimfile = Rcpp::as<string>(bimfile_in);
		string bedfile = Rcpp::as<string>(bedfile_in);
		string outfile = Rcpp::as<string>(outfile_in);
		const size_t npb = Rcpp::as<size_t>(nperbatch_in);
		Rcpp::IntegerVector select(select_in);
		string line, snp;
		char *cp;
		size_t n = res.n_elem, ns = select.size();
		vec g(n);
		uvec gmiss(n), snp_skip = zeros<uvec>(npb);
		mat G(n, npb);
		string *tmpout = new string[npb];
		vector <string> biminfo;
		double gmean, geno, gmax, gmin, num, denom, pval;
		const double tol = 1e-5;
		size_t ncount, nmiss, npbidx = 0, p = 0;
		clock_t time0;
		double compute_time = 0.0;
		ifstream readfile (bimfile.c_str(), ifstream::in);
		if (!readfile) {Rcout << "Error reading bimfile: " << bimfile << "\n"; return R_NilValue;}
		ofstream writefile (outfile.c_str(), ofstream::out);
		if (!writefile) {Rcout << "Error writing outfile: " << outfile << "\n"; return R_NilValue;}
		while(getline(readfile, line)) {
		      cp=strtok ((char *)line.c_str(), " \t");
		      stringstream writeout;
		      snp=cp;
		      writeout << snp << "\t";
		      for(size_t i=0; i<5; ++i) {
			    cp=strtok(NULL, " \t");
			    snp=cp;
			    writeout << snp << "\t";
		      }
		      biminfo.push_back(writeout.str());
		      writeout.clear();
		      p++;
		}
		readfile.close();
		readfile.clear();
		ifstream readbedfile (bedfile.c_str(), ios::binary);
		if (!readbedfile) {Rcout << "Error reading bedfile: " << bedfile << "\n"; return R_NilValue;}
		int nblocks = (ns + 3) / 4, pos;
		unsigned char magic[3], temp[2];
		unsigned char *buffer = new unsigned char[nblocks];
		readbedfile.read((char *)magic, 3);
		if(magic[0] != 0x6C || magic[1] != 0x1B) {Rcout << "Error: " << bedfile << " is not a plink binary file!\n"; return R_NilValue;}
		if(magic[2] != 0x1) {Rcout << "Error: SNP major mode expected in " << bedfile << "\n"; return R_NilValue;}
		writefile << "CHR\tSNP\tcM\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
		for(size_t i=0; i<p; ++i) {
		      stringstream writeout;
		      gmean=0.0;
		      gmax=-100.0;
		      gmin=100.0;
		      nmiss=0;
		      gmiss.zeros();
		      readbedfile.seekg(i*nblocks+3);
		      ncount = 0;
		      readbedfile.read((char *)buffer, nblocks);
		      for(int j=0; j<nblocks; ++j) {
			    pos = 0;
			    for(size_t k=0; k<4; ++k) {
			          if(ncount == ns && j == nblocks - 1) {break;}
				  for(size_t l=0; l<2; ++l) {
				        temp[l] = (buffer[j] >> pos) & 1;
					pos++;
				  }
				  if(select[ncount] <= 0) {ncount++; continue;}
				  if(temp[0] ==0 && temp[1] ==0){
				        geno = 0.0;
				  } else if(temp[0] ==1 && temp[1] ==1){
				        geno = 2.0;
				  } else if(temp[0] ==0 && temp[1] ==1){
				        geno = 1.0;
				  } else {
				        gmiss[select[ncount]-1] = 1;
					nmiss++;
					ncount++;
					continue;
				  }
				  g[select[ncount]-1] = geno;
				  gmean += geno;
				  if(geno>gmax) {gmax=geno;}
				  if(geno<gmin) {gmin=geno;}
				  ncount++;
			    }
		      }
		      gmean/=(double)(n-nmiss);
		      for (size_t j=0; j<n; ++j) {
			    if (gmiss[j]==1) {
			          g[j] = gmean;
				  if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
			    }
			    if (center=='c') {
			          g[j] -= gmean;
			    }
		      }
		      gmean/=2.0; // convert mean to allele freq
		      writeout << biminfo[i] << n-nmiss << "\t" << gmean << "\t";
		      tmpout[npbidx] = writeout.str();
		      writeout.clear();
		      if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			    snp_skip[npbidx] = 1;
		      } else {
			    G.col(npbidx) = g;
		      }
		      npbidx++;
		      if((i+1 == p) || (npbidx == npb)) {
			    uvec snp_idx = find(snp_skip == 0);
			    time0 = clock();
			    vec allnum = G.cols(snp_idx).t() * res;
			    mat alldenom = G.cols(snp_idx).t() * P * G.cols(snp_idx);
			    compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
			    for(size_t j=0; j<npbidx; ++j) {
			          if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
				        writefile << tmpout[j] << "NA\tNA\tNA\n";
				  } else {
				        uvec snp_idx2 = find(snp_idx == j);
					num = allnum[snp_idx2[0]];
					denom = alldenom(snp_idx2[0], snp_idx2[0]);
					if(denom < tol) { // G collinearity with X
					      writefile << tmpout[j] << "0\t0\tNA\n";
					} else {
					      pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
					      writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					}
				  }
			    }
			    if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
		      }
		      if((i+1) % 100000 == 0) {writefile << flush;}
		}
		if(p % 100000 != 0) {writefile << flush;}
		delete [] buffer;
		delete [] tmpout;
		writefile.close();
		writefile.clear();
		readbedfile.close();
		readbedfile.clear();
		return wrap(compute_time);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

  SEXP glmm_score_bed_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP bimfile_in, SEXP bedfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in) {
	try {
		Rcpp::NumericVector res_r(res_in);
		arma::vec res(res_r.begin(), res_r.size(), false);
		arma::sp_mat Sigma_i = Rcpp::as<arma::sp_mat>(Sigma_i_in);
		arma::sp_mat Sigma_iX = Rcpp::as<arma::sp_mat>(Sigma_iX_in);
		arma::sp_mat cov = Rcpp::as<arma::sp_mat>(cov_in);
	        const char center = Rcpp::as<char>(center_in);
		const double minmaf = Rcpp::as<double>(minmaf_in);
		const double maxmaf = Rcpp::as<double>(maxmaf_in);
		const double missrate = Rcpp::as<double>(missrate_in);
		const char miss_method = Rcpp::as<char>(miss_method_in);
		string bimfile = Rcpp::as<string>(bimfile_in);
		string bedfile = Rcpp::as<string>(bedfile_in);
		string outfile = Rcpp::as<string>(outfile_in);
		const size_t npb = Rcpp::as<size_t>(nperbatch_in);
		Rcpp::IntegerVector select(select_in);
		string line, snp;
		char *cp;
		size_t n = res.n_elem, ns = select.size();
		vec g(n);
		uvec gmiss(n), snp_skip = zeros<uvec>(npb);
		mat G(n, npb);
		string *tmpout = new string[npb];
		vector <string> biminfo;
		double gmean, geno, gmax, gmin, num, denom, pval;
		const double tol = 1e-5;
		size_t ncount, nmiss, npbidx = 0, p = 0;
		clock_t time0;
		double compute_time = 0.0;
		ifstream readfile (bimfile.c_str(), ifstream::in);
		if (!readfile) {Rcout << "Error reading bimfile: " << bimfile << "\n"; return R_NilValue;}
		ofstream writefile (outfile.c_str(), ofstream::out);
		if (!writefile) {Rcout << "Error writing outfile: " << outfile << "\n"; return R_NilValue;}
		while(getline(readfile, line)) {
		      cp=strtok ((char *)line.c_str(), " \t");
		      stringstream writeout;
		      snp=cp;
		      writeout << snp << "\t";
		      for(size_t i=0; i<5; ++i) {
			    cp=strtok(NULL, " \t");
			    snp=cp;
			    writeout << snp << "\t";
		      }
		      biminfo.push_back(writeout.str());
		      writeout.clear();
		      p++;
		}
		readfile.close();
		readfile.clear();
		ifstream readbedfile (bedfile.c_str(), ios::binary);
		if (!readbedfile) {Rcout << "Error reading bedfile: " << bedfile << "\n"; return R_NilValue;}
		int nblocks = (ns + 3) / 4, pos;
		unsigned char magic[3], temp[2];
		unsigned char *buffer = new unsigned char[nblocks];
		readbedfile.read((char *)magic, 3);
		if(magic[0] != 0x6C || magic[1] != 0x1B) {Rcout << "Error: " << bedfile << " is not a plink binary file!\n"; return R_NilValue;}
		if(magic[2] != 0x1) {Rcout << "Error: SNP major mode expected in " << bedfile << "\n"; return R_NilValue;}
		writefile << "CHR\tSNP\tcM\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
		for(size_t i=0; i<p; ++i) {
		      stringstream writeout;
		      gmean=0.0;
		      gmax=-100.0;
		      gmin=100.0;
		      nmiss=0;
		      gmiss.zeros();
		      readbedfile.seekg(i*nblocks+3);
		      ncount = 0;
		      readbedfile.read((char *)buffer, nblocks);
		      for(int j=0; j<nblocks; ++j) {
			    pos = 0;
			    for(size_t k=0; k<4; ++k) {
			          if(ncount == ns && j == nblocks - 1) {break;}
				  for(size_t l=0; l<2; ++l) {
				        temp[l] = (buffer[j] >> pos) & 1;
					pos++;
				  }
				  if(select[ncount] <= 0) {ncount++; continue;}
				  if(temp[0] ==0 && temp[1] ==0){
				        geno = 0.0;
				  } else if(temp[0] ==1 && temp[1] ==1){
				        geno = 2.0;
				  } else if(temp[0] ==0 && temp[1] ==1){
				        geno = 1.0;
				  } else {
				        gmiss[select[ncount]-1] = 1;
					nmiss++;
					ncount++;
					continue;
				  }
				  g[select[ncount]-1] = geno;
				  gmean += geno;
				  if(geno>gmax) {gmax=geno;}
				  if(geno<gmin) {gmin=geno;}
				  ncount++;
			    }
		      }
		      gmean/=(double)(n-nmiss);
		      for (size_t j=0; j<n; ++j) {
			    if (gmiss[j]==1) {
			          g[j] = gmean;
				  if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
			    }
			    if (center=='c') {
			          g[j] -= gmean;
			    }
		      }
		      gmean/=2.0; // convert mean to allele freq
		      writeout << biminfo[i] << n-nmiss << "\t" << gmean << "\t";
		      tmpout[npbidx] = writeout.str();
		      writeout.clear();
		      if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
			    snp_skip[npbidx] = 1;
		      } else {
			    G.col(npbidx) = g;
		      }
		      npbidx++;
		      if((i+1 == p) || (npbidx == npb)) {
			    uvec snp_idx = find(snp_skip == 0);
			    time0 = clock();
			    vec allnum = G.cols(snp_idx).t() * res;
			    sp_mat Gsp(G.cols(snp_idx));
			    sp_mat XSigma_iG = Sigma_iX.t() * Gsp;
			    sp_mat alldenom = Gsp.t() * Sigma_i * Gsp - XSigma_iG.t() * cov * XSigma_iG;
			    compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
			    for(size_t j=0; j<npbidx; ++j) {
			          if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
				        writefile << tmpout[j] << "NA\tNA\tNA\n";
				  } else {
				        uvec snp_idx2 = find(snp_idx == j);
					num = allnum[snp_idx2[0]];
					denom = alldenom(snp_idx2[0], snp_idx2[0]);
					if(denom < tol) { // G collinearity with X
					      writefile << tmpout[j] << "0\t0\tNA\n";
					} else {
					      pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
					      writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
					}
				  }
			    }
			    if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
		      }
		      if((i+1) % 100000 == 0) {writefile << flush;}
		}
		if(p % 100000 != 0) {writefile << flush;}
		delete [] buffer;
		delete [] tmpout;
		writefile.close();
		writefile.clear();
		readbedfile.close();
		readbedfile.clear();
		return wrap(compute_time);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

SEXP glmm_wald_text(SEXP n_in, SEXP snp_in, SEXP infile_in, SEXP tol_in, SEXP center_in, SEXP miss_method_in, SEXP p_in, SEXP infile_nrow_skip_in, SEXP infile_sep_in, SEXP infile_na_in, SEXP infile_ncol_skip_in, SEXP infile_ncol_print_in, SEXP snp_col_in, SEXP select_in) {
	try {
	        const size_t n = Rcpp::as<size_t>(n_in);
		string snp = Rcpp::as<string>(snp_in);
		const double tol = Rcpp::as<double>(tol_in);
		const char center = Rcpp::as<char>(center_in);
		const char miss_method = Rcpp::as<char>(miss_method_in);
		string infile = Rcpp::as<string>(infile_in);
		const size_t p = Rcpp::as<size_t>(p_in);
		const size_t infile_nrow_skip = Rcpp::as<size_t>(infile_nrow_skip_in);
		string infile_sep = Rcpp::as<string>(infile_sep_in);
		string infile_na = Rcpp::as<string>(infile_na_in);
		const int infile_ncol_skip = Rcpp::as<int>(infile_ncol_skip_in);
		Rcpp::IntegerVector infile_ncol_print(infile_ncol_print_in);
		const size_t ncol_print = infile_ncol_print.size();
		Rcpp::CharacterVector infile_header_print(ncol_print);
		const int snp_col = Rcpp::as<int>(snp_col_in);
		Rcpp::IntegerVector select(select_in);
		char *cp;
		vec g(n);
		uvec gmiss(n);
		double gmean, geno, gmax, gmin;
		size_t nmiss, infile_ncol_print_idx, skip = 0, ns = select.size();
		bool snpfound = false;
		size_t tmppos = infile.rfind('.');
		string ext = infile.substr(tmppos);
		if (strcmp(ext.c_str(), ".bz2")==0) { // bzip2 infile
		      FILE* fp = fopen(infile.c_str(), "rb");
		      int bzerror;
		      BZFILE* bfp = BZ2_bzReadOpen(&bzerror, fp, 0, 0, NULL, 0);;
		      int buffer_length = 100 * (int(ns) + infile_ncol_skip);
		      if(bzerror != BZ_OK) {
			    BZ2_bzReadClose(&bzerror, bfp);
			    fclose(fp);
			    Rcout << "File " << infile << " appears not to be compressed by bzip2\n";
			    return R_NilValue;
		      }
		      for(size_t i=0; i<p; ++i) {
			    int nread = 0;
			    char *buffer = new char[buffer_length];
			    while (bzerror == BZ_OK) {
			          int nn = BZ2_bzRead(&bzerror, bfp, buffer + nread, 1);
				  if (bzerror == BZ_STREAM_END) {
				        char *unused, *next_unused = NULL;
					int nUnused;
					BZ2_bzReadGetUnused(&bzerror, bfp, (void**) &unused, &nUnused);
					if (bzerror == BZ_OK) {
					      if (nUnused > 0) {
						    next_unused = (char*) malloc(nUnused);
						    if (!next_unused) {
						          BZ2_bzReadClose(&bzerror, bfp);
							  fclose(fp);
						          Rcout << "Allocation of overflow buffer for bzfile failed\n";
							  return R_NilValue;
						    }
						    memcpy(next_unused, unused, nUnused);
					      }
					      if (nUnused > 0 || !feof(fp)) {
						    BZ2_bzReadClose(&bzerror, bfp);
						    bfp = BZ2_bzReadOpen(&bzerror, fp, 0, 0, next_unused, nUnused);
						    if(bzerror != BZ_OK) {
						          BZ2_bzReadClose(&bzerror, bfp);
							  fclose(fp);
						          Rcout << "File " << infile << " has trailing content that appears not to be compressed by bzip2\n";
							  return R_NilValue;
						    }
					      }
					      if (next_unused) free(next_unused);
					}
				  } else if (bzerror != BZ_OK) {
				        nread += nn;
					break;
				  }
				  nread += nn;
				  if(nread > 0 && buffer[nread-1] == '\n') {
				        buffer[nread-1] = '\0';
					break;
				  }
			    }
			    if(i<infile_nrow_skip) {
			          delete [] buffer;
				  continue;
			    }
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    cp=strtok (buffer, infile_sep.c_str());
			    infile_ncol_print_idx = 0;
			    if (infile_ncol_print[0] == 1) {
				  infile_header_print[0] = cp;
				  if (snp_col == 1 && strcmp(cp, snp.c_str()) == 0) {
				        snpfound = true;
				  }
				  infile_ncol_print_idx++;
			    }
			    if (infile_ncol_skip > 1) {
			          for(int k=1; k<infile_ncol_skip; ++k) {
				        cp=strtok (NULL, infile_sep.c_str());
					if (snp_col == k+1 && strcmp(cp, snp.c_str()) == 0) {
					      snpfound = true;
					}
					if (infile_ncol_print_idx<ncol_print) {
					      if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						    infile_header_print[infile_ncol_print_idx] = cp;
						    infile_ncol_print_idx++;
					      }
					}
				  }
			    }
			    if (snpfound) {
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
				  gmean/=(double)(n-nmiss);
				  for (size_t j=0; j<n; ++j) {
				        if (gmiss[j]==1) {
					      if(miss_method=='o') {
						    g[j] = -999.9;
					      } else {
						    g[j] = gmean;
					      }
					}
					if (center=='c' && (miss_method!='o' || gmiss[j]!=1)) {
					      g[j] -= gmean;
					}
				  }
				  gmean/=2.0; // convert mean to allele freq
				  if(gmax-gmin<tol) { // monomorphic
				        skip = 1;
				  }
				  delete [] buffer;
				  break;
			    }
			    delete [] buffer;
		      }
		      BZ2_bzReadClose(&bzerror, bfp);
		      fclose(fp);
		} else if (strcmp(ext.c_str(), ".gz")==0) { // gzip infile
		      gzFile fp = gzopen(infile.c_str(), "rb");
		      int buffer_length = 100 * (int(ns) + infile_ncol_skip);
		      if(!fp) {
			    Rcout << "Error: cannot open gzipped file " << infile << "\n";
			    return R_NilValue;
		      }
		      for(size_t i=0; i<p; ++i) {
			    char *buffer = new char[buffer_length];
			    gzgets(fp, buffer, buffer_length);
			    if(i<infile_nrow_skip) {
			          delete [] buffer;
				  continue;
			    }
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    cp=strtok (buffer, infile_sep.c_str());
			    infile_ncol_print_idx = 0;
			    if (infile_ncol_print[0] == 1) {
				  infile_header_print[0] = cp;
				  if (snp_col == 1 && strcmp(cp, snp.c_str()) == 0) {
				        snpfound = true;
				  }
				  infile_ncol_print_idx++;
			    }
			    if (infile_ncol_skip > 1) {
			          for(int k=1; k<infile_ncol_skip; ++k) {
				        cp=strtok (NULL, infile_sep.c_str());
					if (snp_col == k+1 && strcmp(cp, snp.c_str()) == 0) {
					      snpfound = true;
					}
					if (infile_ncol_print_idx<ncol_print) {
					      if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						    infile_header_print[infile_ncol_print_idx] = cp;
						    infile_ncol_print_idx++;
					      }
					}
				  }
			    }
			    if (snpfound) {
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
				  gmean/=(double)(n-nmiss);
				  for (size_t j=0; j<n; ++j) {
				        if (gmiss[j]==1) {
					      if(miss_method=='o') {
						    g[j] = -999.9;
					      } else {
						    g[j] = gmean;
					      }
					}
					if (center=='c' && (miss_method!='o' || gmiss[j]!=1)) {
					      g[j] -= gmean;
					}
				  }
				  gmean/=2.0; // convert mean to allele freq
				  if(gmax-gmin<tol) { // monomorphic
				        skip = 1;
				  }
				  delete [] buffer;
				  break;
			    }
			    delete [] buffer;
		      }
		      gzclose(fp);
		} else { // plain text infile
		      ifstream readfile (infile.c_str(), ifstream::in);
		      if (!readfile) {Rcout << "Error reading genotype file: " << infile << "\n"; return R_NilValue;}
		      string line;
		      for(size_t i=0; i<p; ++i) {
			    getline(readfile, line);
			    if(i<infile_nrow_skip) {continue;}
			    gmean=0.0;
			    gmax=-100.0;
			    gmin=100.0;
			    nmiss=0;
			    gmiss.zeros();
			    cp=strtok ((char *)line.c_str(), infile_sep.c_str());
			    infile_ncol_print_idx = 0;
			    if (infile_ncol_print[0] == 1) {
				  infile_header_print[0] = cp;
				  if (snp_col == 1 && strcmp(cp, snp.c_str()) == 0) {
				        snpfound = true;
				  }
				  infile_ncol_print_idx++;
			    }
			    if (infile_ncol_skip > 1) {
			          for(int k=1; k<infile_ncol_skip; ++k) {
				        cp=strtok (NULL, infile_sep.c_str());
					if (snp_col == k+1 && strcmp(cp, snp.c_str()) == 0) {
					      snpfound = true;
					}
					if (infile_ncol_print_idx<ncol_print) {
					      if(infile_ncol_print[infile_ncol_print_idx] == k+1) {
						    infile_header_print[infile_ncol_print_idx] = cp;
						    infile_ncol_print_idx++;
					      }
					}
				  }
			    }
			    if (snpfound) {
				  for (size_t j=0; j<ns; ++j) {
				        cp=strtok (NULL, infile_sep.c_str());
					if(select[j] <= 0) {continue;}
					if (strcmp(cp, infile_na.c_str())==0) {
					      gmiss[select[j]-1] = 1;
					      nmiss++;
					} else {
					      geno = atof(cp); 				
					      g[select[j]-1] = geno;
					      gmean += geno;
					      if(geno>gmax) {gmax=geno;}
					      if(geno<gmin) {gmin=geno;}
					}
				  }	
				  gmean/=(double)(n-nmiss);
				  for (size_t j=0; j<n; ++j) {
				        if (gmiss[j]==1) {
					      if(miss_method=='o') {
						    g[j] = -999.9;
					      } else {
						    g[j] = gmean;
					      }
					}
					if (center=='c' && (miss_method!='o' || gmiss[j]!=1)) {
					      g[j] -= gmean;
					}
				  }
				  gmean/=2.0; // convert mean to allele freq
				  if(gmax-gmin<tol) { // monomorphic
				        skip = 1;
				  }
				  break;
			    }
		      }
		      readfile.close();
		      readfile.clear();
		}
		if (!snpfound) {skip = 2;}
		return List::create(Named("snpinfo") = infile_header_print, Named("N") = n-nmiss, Named("AF") = gmean, Named("G") = g, Named("skip") = skip);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}

SEXP glmm_wald_bed(SEXP n_in, SEXP snp_in, SEXP bimfile_in, SEXP bedfile_in, SEXP center_in, SEXP miss_method_in, SEXP select_in) {
	try {
	        const size_t n = Rcpp::as<size_t>(n_in);
		string snp = Rcpp::as<string>(snp_in);
		const char center = Rcpp::as<char>(center_in);
		const char miss_method = Rcpp::as<char>(miss_method_in);
		string bimfile = Rcpp::as<string>(bimfile_in);
		string bedfile = Rcpp::as<string>(bedfile_in);
		Rcpp::IntegerVector select(select_in);
		Rcpp::CharacterVector infile_header_print(6);
		string line;
		char *cp;
		vec g(n);
		uvec gmiss(n);
		double gmean, geno, gmax, gmin;
		const double tol = 1e-5;
		size_t ncount, nmiss=0, p=0, skip = 0, ns = select.size();
		bool snpfound = false;
		ifstream readfile (bimfile.c_str(), ifstream::in);
		if (!readfile) {Rcout << "Error reading bimfile: " << bimfile << "\n"; return R_NilValue;}
		while(getline(readfile, line)) {
		      cp=strtok ((char *)line.c_str(), " \t");
		      infile_header_print[0] = cp;
		      for(size_t i=0; i<5; ++i) {
			    cp=strtok(NULL, " \t");
			    infile_header_print[i+1]=cp;
			    if(i == 0 && strcmp(cp, snp.c_str()) == 0) {snpfound = true;}
		      }
		      if(snpfound) {break;}
		      p++;
		}
		readfile.close();
		readfile.clear();
		if(snpfound) {
		      ifstream readfile (bedfile.c_str(), ios::binary);
		      if (!readfile) {Rcout << "Error reading bedfile: " << bedfile << "\n"; return R_NilValue;}
		      int nblocks = (ns + 3) / 4, pos;
		      unsigned char magic[3], temp[2];
		      unsigned char *buffer = new unsigned char[nblocks];
		      readfile.read((char *)magic, 3);
		      if(magic[0] != 0x6C || magic[1] != 0x1B) {Rcout << "Error: " << bedfile << " is not a plink binary file!\n"; return R_NilValue;}
		      if(magic[2] != 0x1) {Rcout << "Error: SNP major mode expected in " << bedfile << "\n"; return R_NilValue;}
		      gmean=0.0;
		      gmax=-100.0;
		      gmin=100.0;
		      nmiss=0;
		      gmiss.zeros();
		      readfile.seekg(p*nblocks+3);
		      ncount = 0;
		      readfile.read((char *)buffer, nblocks);
		      for(int j=0; j<nblocks; ++j) {
			    pos = 0;
			    for(size_t k=0; k<4; ++k) {
			          if(ncount == ns && j == nblocks - 1) {break;}
				  for(size_t l=0; l<2; ++l) {
				        temp[l] = (buffer[j] >> pos) & 1;
					pos++;
				  }
				  if(select[ncount] <= 0) {ncount++; continue;}
				  if(temp[0] ==0 && temp[1] ==0){
				        geno = 0.0;
				  } else if(temp[0] ==1 && temp[1] ==1){
				        geno = 2.0;
				  } else if(temp[0] ==0 && temp[1] ==1){
				        geno = 1.0;
				  } else {
				        gmiss[select[ncount]-1] = 1;
					nmiss++;
					ncount++;
					continue;
				  }
				  g[select[ncount]-1] = geno;
				  gmean += geno;
				  if(geno>gmax) {gmax=geno;}
				  if(geno<gmin) {gmin=geno;}
				  ncount++;
			    }
		      }
		      delete [] buffer;
		      gmean/=(double)(n-nmiss);
		      for (size_t j=0; j<n; ++j) {
			    if (gmiss[j]==1) {
			          if(miss_method=='o') {
				        g[j] = -999.9;
				  } else {
				        g[j] = gmean;
				  }
			    }
			    if (center=='c' && (miss_method!='o' || gmiss[j]!=1)) {
			          g[j] -= gmean;
			    }
		      }
		      gmean/=2.0; // convert mean to allele freq
		      if(gmax-gmin<tol) { // monomorphic
			    skip = 1;
		      }
		      readfile.close();
		      readfile.clear();
		} else {skip = 2;}
		return List::create(Named("snpinfo") = infile_header_print, Named("N") = n-nmiss, Named("AF") = gmean, Named("G") = g, Named("skip") = skip);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}
  
  
    
  SEXP glmm_score_bgen13(SEXP res_in, SEXP P_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
    try{
      
      Rcpp::NumericMatrix P_r(P_in);
      Rcpp::NumericVector res_r(res_in);
    	size_t n_P = P_r.nrow(), k_P = P_r.ncol();
    	arma::mat P(P_r.begin(), n_P, k_P, false);
      arma::vec res(res_r.begin(), res_r.size(), false);
      const char center = Rcpp::as<char>(center_in);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const char miss_method = Rcpp::as<char>(miss_method_in);
      
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      string outfile  = Rcpp::as<string>(outfile_in);
      
      const size_t npb = Rcpp::as<size_t>(nperbatch_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n), snp_skip = zeros<uvec>(npb);
      mat G(n, npb);
      string* tmpout = new string[npb];
      vector <string> biminfo;
      double gmean, geno, gmax, gmin, num, denom, pval;
      const double tol = 1e-5;
      size_t ncount, nmiss, npbidx = 0;
      clock_t time0;
      double compute_time = 0.0;
      ofstream writefile(outfile.c_str(), ofstream::out);
      
      

      uint maxLA = 65536;
      uint maxLB = 65536;
      std::vector<uchar> zBuf12;
      std::vector<uchar> shortBuf12;
      uint compression = Rcpp::as<uint>(compression_in);
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];
      struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
      
      uint begin = Rcpp::as<uint>(begin_in);
      uint end   = Rcpp::as<uint>(end_in);
      long long unsigned int byte = Rcpp::as<long long unsigned int>(pos_in);
      FILE* fp = fopen(bgenfile.c_str(), "rb");
      fseek(fp, byte, SEEK_SET);
      
      bool isMultiThread = Rcpp::as<bool>(isMultiThread_in);
      if (!isMultiThread){
        writefile << "SNP\tRSID\tCHR\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
      }
      
      int ret;
      for (uint m = begin; m < end; m++) {
        stringstream writeout;
        maxLA = 65536;
        maxLB = 65536;
        
        ushort LS; ret = fread(&LS, 2, 1, fp);
        ret = fread(snpID, 1, LS, fp); snpID[LS] = '\0';
        
        ushort LR; ret = fread(&LR, 2, 1, fp);
        ret = fread(rsID, 1, LR, fp);  rsID[LR] = '\0';
        
        ushort LC; ret = fread(&LC, 2, 1, fp);
        ret = fread(chrStr, 1, LC, fp); chrStr[LC] = '\0';
        
        uint physpos; ret = fread(&physpos, 4, 1, fp);
	string physpos_tmp = to_string(physpos);
        ushort LKnum; ret = fread(&LKnum, 2, 1, fp);
        if (LKnum != 2) {Rcout << "Error reading BGEN file: There are non-bi-allelic variants. \n"; return R_NilValue;}
        uint32_t LA; ret = fread(&LA, 4, 1, fp);
        if (LA > maxLA) {
          maxLA = 2 * LA;
          delete[] allele1;
          allele1 = new char[maxLA + 1];
        }
        ret = fread(allele1, 1, LA, fp); allele1[LA] = '\0';

        uint32_t LB; ret = fread(&LB, 4, 1, fp);
        if (LB > maxLB) {
          maxLB = 2 * LB;
          delete[] allele0;
          allele0 = new char[maxLB + 1];
        }
        ret = fread(allele0, 1, LB, fp); allele0[LB] = '\0';
        
        uint cLen; ret = fread(&cLen, 4, 1, fp);
        
        uchar* bufAt;
        if (compression == 1) {
          zBuf12.resize(cLen - 4);
          uint dLen; ret = fread(&dLen, 4, 1, fp);
          ret = fread(&zBuf12[0], 1, cLen - 4, fp);
          shortBuf12.resize(dLen);

          uLongf destLen = dLen;
          if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
			  Rcout << "Error reading bgen file: Decompressing variant block failed with libdeflate. \n"; return R_NilValue;
          }
          bufAt = &shortBuf12[0];
        }
        else if (compression == 2) {
          zBuf12.resize(cLen - 4);
          uint dLen; ret = fread(&dLen, 4, 1, fp);
          ret = fread(&zBuf12[0], 1, cLen - 4, fp);
          shortBuf12.resize(dLen);

          uLongf destLen = dLen;
          size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
          if (ret > destLen) {
            if (ZSTD_isError(ret)) {
              Rcout << "Error reading bgen file: Decompressing genotype block failed. \n"; return R_NilValue;
            }
          }
          bufAt = &shortBuf12[0];
        }
        else {
          zBuf12.resize(cLen);
          ret = fread(&zBuf12[0], 1, cLen, fp);
          bufAt = &zBuf12[0];
        }
      
        uint32_t N; memcpy(&N, bufAt, sizeof(int32_t));
        uint16_t K; memcpy(&K, &(bufAt[4]), sizeof(int16_t));
        if (K != 2) {Rcout << "Error reading bgen file: There are variants with more than 2 alleles. \n"; return R_NilValue;}
        const uint32_t min_ploidy = bufAt[6];
        if (min_ploidy != 2) {Rcout << "Error reading bgen file: Minimum ploidy should be 2. \n"; return R_NilValue;}
        const uint32_t max_ploidy = bufAt[7];
        if (max_ploidy != 2) {Rcout << "Error reading bgen file: Maximum ploidy should be 2. \n"; return R_NilValue;}
        
        const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
        const unsigned char* probs_start = &(bufAt[10 + N]);
        const uint32_t is_phased = probs_start[-2];
        if (is_phased < 0 || is_phased > 1) {Rcout << "Error reading bgen file: Phased value must be 0 or 1. \n"; return R_NilValue;}
        
        const uint32_t B = probs_start[-1];
        
	if (B != 8 && B!= 16 && B !=24 && B != 32) {
            Rcout << "Error reading bgen file: Bits to store probabilities must be 8, 16, 24, or 32. \n"; return R_NilValue;
	}
	
	const uintptr_t numer_mask = (1U << B) - 1;
        const uintptr_t probs_offset = B / 8;
        
        gmean=0.0;
        gmax=-100.0;
        gmin=100.0;
        nmiss=0;
        ncount = 0;
        
        if (!is_phased) {
          for (size_t i = 0; i < N; i++) {
            const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
            uintptr_t numer_aa;
            uintptr_t numer_ab;
            
            if (missing_and_ploidy == 2){
                Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
                probs_start += (probs_offset * 2);
                      
            } else if (missing_and_ploidy == 130){
                probs_start += (probs_offset * 2);
                gmiss[select[ncount]-1] = 1;
                nmiss++;
                ncount++;
                continue;
              
            } else {
		Rcout << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n"; return R_NilValue;    
            }
            
            
            if (select[ncount] > 0){ 
                double p11 = numer_aa / double(1.0 * numer_mask);
                double p10 = numer_ab / double(1.0 * numer_mask);
                geno = 2 * (1 - p11 - p10) + p10;
                gmiss[select[ncount]-1] = 0;
                g[select[ncount]-1] = geno;
                gmean += geno;
                if (geno > gmax) { gmax = geno; }
                if (geno < gmin) { gmin = geno; }
            } 
            ncount++;
            
          }
          
       } else {
         for (size_t i = 0; i < N; i++) {
           const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
           uintptr_t numer_aa;
           uintptr_t numer_ab;
           
           if (missing_and_ploidy == 2){
               Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
               probs_start += (probs_offset * 2);
             
           } else if (missing_and_ploidy == 130){
               probs_start += (probs_offset * 2);
               gmiss[select[ncount]-1] = 1;
               nmiss++;
               ncount++;
               continue;
             
           } else {
	       Rcout << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n"; return R_NilValue;
           }
           
           if (select[ncount] > 0){ 
               double p11 = numer_aa / double(1.0 * numer_mask);
               double p10 = numer_ab / double(1.0 * numer_mask);
               geno = double(1.0 * missing_and_ploidy) - (p11 + p10);
               gmiss[select[ncount]-1] = 0;
               g[select[ncount]-1] = geno;
               gmean += geno;
               if (geno > gmax) { gmax = geno; }
               if (geno < gmin) { gmin = geno; }
           } 
           ncount++;
           
         }
       }
        
       gmean/=(double)(n-nmiss);
       for (size_t j=0; j<n; ++j) {
        if (gmiss[j]==1) {
          g[j] = gmean;
          if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
        }
        if (center=='c') {
          g[j] -= gmean;
        }
       }

       gmean /= 2.0; // convert mean to allele freq
       writeout << snpID << "\t" << rsID << "\t" << chrStr << "\t" << physpos_tmp << "\t" << allele1 << "\t" << allele0 << "\t" << (n-nmiss) << "\t" << gmean << "\t";
       tmpout[npbidx] = writeout.str();
       writeout.clear();
       if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
        snp_skip[npbidx] = 1;
       } else {
        G.col(npbidx) = g;
       }
       npbidx++;


       if((m+1 == end) || (npbidx == npb)) {
        uvec snp_idx = find(snp_skip == 0);
        time0 = clock();
        vec allnum = G.cols(snp_idx).t() * res;
        mat alldenom = G.cols(snp_idx).t() * P * G.cols(snp_idx);
        compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
        for(size_t j=0; j<npbidx; ++j) {
          if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
            writefile << tmpout[j] << "NA\tNA\tNA\n";
          } else {
            uvec snp_idx2 = find(snp_idx == j);
            num = allnum[snp_idx2[0]];
            denom = alldenom(snp_idx2[0], snp_idx2[0]);
            if(denom < tol) { // G collinearity with X
              writefile << tmpout[j] << "0\t0\tNA\n";
            } else {
              pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
              writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
            }
          }
        }
        if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
       }
        if((m+1) % 100000 == 0) {writefile << flush;}
     }

     if(end % 100000 != 0) {writefile << flush;}
     (void)ret;
     libdeflate_free_decompressor(decompressor);
     delete [] tmpout;
     delete [] snpID;
     delete [] rsID;
     delete [] chrStr;
     delete [] allele0;
     delete [] allele1;
 
     writefile.close();
     writefile.clear();
     fclose(fp);
     return wrap(compute_time);
    } 
     catch( std::exception &ex ) {
     forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;
  }  
  
  
  SEXP glmm_score_bgen13_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
    try{
      
      Rcpp::NumericVector res_r(res_in);
      arma::vec res(res_r.begin(), res_r.size(), false);
      arma::sp_mat Sigma_i = Rcpp::as<arma::sp_mat>(Sigma_i_in);
      arma::sp_mat Sigma_iX = Rcpp::as<arma::sp_mat>(Sigma_iX_in);
      arma::sp_mat cov = Rcpp::as<arma::sp_mat>(cov_in);
      const char center = Rcpp::as<char>(center_in);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const char miss_method = Rcpp::as<char>(miss_method_in);
      
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      string outfile  = Rcpp::as<string>(outfile_in);
      
      const size_t npb = Rcpp::as<size_t>(nperbatch_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n), snp_skip = zeros<uvec>(npb);
      mat G(n, npb);
      string* tmpout = new string[npb];
      vector <string> biminfo;
      double gmean, geno, gmax, gmin, num, denom, pval;
      const double tol = 1e-5;
      size_t ncount, nmiss, npbidx = 0;
      clock_t time0;
      double compute_time = 0.0;
      ofstream writefile(outfile.c_str(), ofstream::out);
      
      
      
      uint maxLA = 65536;
      uint maxLB = 65536;
      std::vector<uchar> zBuf12;
      std::vector<uchar> shortBuf12;
      uint compression = Rcpp::as<uint>(compression_in);
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];
      struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
      
      uint begin = Rcpp::as<uint>(begin_in);
      uint end   = Rcpp::as<uint>(end_in);
      long long unsigned int byte = Rcpp::as<long long unsigned int>(pos_in);
      FILE* fp = fopen(bgenfile.c_str(), "rb");
      fseek(fp, byte, SEEK_SET);
      
      bool isMultiThread = Rcpp::as<bool>(isMultiThread_in);
      if (!isMultiThread){
        writefile << "SNP\tRSID\tCHR\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
      }
      
      int ret;
      for (uint m = begin; m < end; m++) {
        stringstream writeout;
        maxLA = 65536;
        maxLB = 65536;
        
        ushort LS; ret = fread(&LS, 2, 1, fp);
        ret = fread(snpID, 1, LS, fp); snpID[LS] = '\0';
        
        ushort LR; ret = fread(&LR, 2, 1, fp);
        ret = fread(rsID, 1, LR, fp);  rsID[LR] = '\0';
        
        ushort LC; ret = fread(&LC, 2, 1, fp);
        ret = fread(chrStr, 1, LC, fp); chrStr[LC] = '\0';
        
        uint physpos; ret = fread(&physpos, 4, 1, fp);
      	string physpos_tmp = to_string(physpos);
        ushort LKnum; ret = fread(&LKnum, 2, 1, fp);
        if (LKnum != 2) {Rcout << "Error reading BGEN file: There are non-bi-allelic variants. \n"; return R_NilValue;}
        uint32_t LA; ret = fread(&LA, 4, 1, fp);
        if (LA > maxLA) {
          maxLA = 2 * LA;
          delete[] allele1;
          allele1 = new char[maxLA + 1];
        }
        ret = fread(allele1, 1, LA, fp); allele1[LA] = '\0';
        
        uint32_t LB; ret = fread(&LB, 4, 1, fp);
        if (LB > maxLB) {
          maxLB = 2 * LB;
          delete[] allele0;
          allele0 = new char[maxLB + 1];
        }
        ret = fread(allele0, 1, LB, fp); allele0[LB] = '\0';
        
        uint cLen; ret = fread(&cLen, 4, 1, fp);
        
        uchar* bufAt;
        if (compression == 1) {
          zBuf12.resize(cLen - 4);
          uint dLen; ret = fread(&dLen, 4, 1, fp);
          ret = fread(&zBuf12[0], 1, cLen - 4, fp);
          shortBuf12.resize(dLen);
          
          uLongf destLen = dLen;
          if (libdeflate_zlib_decompress(decompressor, &zBuf12[0], cLen - 4, &shortBuf12[0], destLen, NULL) != LIBDEFLATE_SUCCESS) {
	      Rcout << "Error reading bgen file: Decompressing variant block failed with libdeflate. \n"; return R_NilValue;
          }
          bufAt = &shortBuf12[0];
        }
        else if (compression == 2) {
          zBuf12.resize(cLen - 4);
          uint dLen; ret = fread(&dLen, 4, 1, fp);
          ret = fread(&zBuf12[0], 1, cLen - 4, fp);
          shortBuf12.resize(dLen);
          
          uLongf destLen = dLen;
          size_t ret = ZSTD_decompress(&shortBuf12[0], destLen, &zBuf12[0], cLen - 4);
          if (ret > destLen) {
            if (ZSTD_isError(ret)) {
              Rcout << "Error reading bgen file: Decompressing genotype block failed. \n"; return R_NilValue;
            }
          }
          bufAt = &shortBuf12[0];
        }
        else {
          zBuf12.resize(cLen);
          ret = fread(&zBuf12[0], 1, cLen, fp);
          bufAt = &zBuf12[0];
        }
        
        uint32_t N; memcpy(&N, bufAt, sizeof(int32_t));
        uint16_t K; memcpy(&K, &(bufAt[4]), sizeof(int16_t));
        if (K != 2) {Rcout << "Error reading bgen file: There are variants with more than 2 alleles. \n"; return R_NilValue;}
        const uint32_t min_ploidy = bufAt[6];
        if (min_ploidy != 2) {Rcout << "Error reading bgen file: Minimum ploidy must be 2. \n"; return R_NilValue;}
        const uint32_t max_ploidy = bufAt[7];
        if (max_ploidy != 2) {Rcout << "Error reading bgen file: Maximum ploidy must be 2. \n"; return R_NilValue;}
        
        const unsigned char* missing_and_ploidy_info = &(bufAt[8]);
        const unsigned char* probs_start = &(bufAt[10 + N]);
        const uint32_t is_phased = probs_start[-2];
        if (is_phased < 0 || is_phased > 1) {Rcout << "Error reading bgen file: Phased value must be 0 or 1. \n"; return R_NilValue;}
        
        const uint32_t B = probs_start[-1];
        if (B != 8 && B!= 16 && B !=24 && B != 32) {
          Rcout << "Error reading bgen file: Bits to store probabilities must be 8, 16, 24, or 32. \n"; return R_NilValue;
        }
        const uintptr_t numer_mask = (1U << B) - 1;
        const uintptr_t probs_offset = B / 8;

        gmean=0.0;
        gmax=-100.0;
        gmin=100.0;
        nmiss=0;
        ncount = 0;
        
        if (!is_phased) {
            for (size_t i = 0; i < N; i++) {
              const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
              uintptr_t numer_aa;
              uintptr_t numer_ab;
              
              if (missing_and_ploidy == 2){
                  Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
                  probs_start += (probs_offset * 2);
                
              } else if (missing_and_ploidy == 130){
                  probs_start += (probs_offset * 2);
                  gmiss[select[ncount]-1] = 1;
                  nmiss++;
                  ncount++;
                  continue;
              } else {
                  Rcout << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n"; return R_NilValue;
	      }
              
              
              if (select[ncount] > 0){ 
                  double p11 = numer_aa / double(1.0 * numer_mask);
                  double p10 = numer_ab / double(1.0 * numer_mask);
                  geno = 2 * (1 - p11 - p10) + p10;
                  gmiss[select[ncount]-1] = 0;
                  g[select[ncount]-1] = geno;
                  gmean += geno;
                  if (geno > gmax) { gmax = geno; }
                  if (geno < gmin) { gmin = geno; }
              } 
              ncount++;
              
            }
            
        } else {
          for (size_t i = 0; i < N; i++) {
            const uint32_t missing_and_ploidy = missing_and_ploidy_info[i];
            uintptr_t numer_aa;
            uintptr_t numer_ab;
            
            if (missing_and_ploidy == 2){
              Bgen13GetTwoVals(probs_start, B, probs_offset, &numer_aa, &numer_ab);
              probs_start += (probs_offset * 2);
              
            } else if (missing_and_ploidy == 130){
              probs_start += (probs_offset * 2);
              gmiss[select[ncount]-1] = 1;
              nmiss++;
              ncount++;
              continue;
              
            } else {
	      Rcout << "Error reading bgen file: Ploidy value " << missing_and_ploidy << " is unsupported. Must be 2 or 130. \n"; return R_NilValue;
            }
            
            if (select[ncount] > 0){ 
                double p11 = numer_aa / double(1.0 * numer_mask);
                double p10 = numer_ab / double(1.0 * numer_mask);
                geno = double(1.0 * missing_and_ploidy) - (p11 + p10);
                gmiss[select[ncount]-1] = 0;
                g[select[ncount]-1] = geno;
                gmean += geno;
                if (geno > gmax) { gmax = geno; }
                if (geno < gmin) { gmin = geno; }
            } 
            ncount++;
            
          }
        }
        
        gmean/=(double)(n-nmiss);
        for (size_t j=0; j<n; ++j) {
          if (gmiss[j]==1) {
            g[j] = gmean;
            if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
          }
          if (center=='c') {
            g[j] -= gmean;
          }
        }
        
        gmean /= 2.0; // convert mean to allele freq
        writeout << snpID << "\t" << rsID << "\t" << chrStr << "\t" << physpos_tmp << "\t" << allele1 << "\t" << allele0 << "\t" << (n-nmiss) << "\t" << gmean << "\t";
        tmpout[npbidx] = writeout.str();
        writeout.clear();
        if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
          snp_skip[npbidx] = 1;
        } else {
          G.col(npbidx) = g;
        }
        npbidx++;
        
        
        if((m+1 == end) || (npbidx == npb)) {
          uvec snp_idx = find(snp_skip == 0);
          time0 = clock();
          vec allnum = G.cols(snp_idx).t() * res;
          sp_mat Gsp(G.cols(snp_idx));
          sp_mat XSigma_iG = Sigma_iX.t() * Gsp;
          sp_mat alldenom = Gsp.t() * Sigma_i * Gsp - XSigma_iG.t() * cov * XSigma_iG;
          compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
          for(size_t j=0; j<npbidx; ++j) {
            if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
              writefile << tmpout[j] << "NA\tNA\tNA\n";
            } else {
              uvec snp_idx2 = find(snp_idx == j);
              num = allnum[snp_idx2[0]];
              denom = alldenom(snp_idx2[0], snp_idx2[0]);
              if(denom < tol) { // G collinearity with X
                writefile << tmpout[j] << "0\t0\tNA\n";
              } else {
                pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
                writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
              }
            }
          }
          if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
        }
        if((m+1) % 100000 == 0) {writefile << flush;}
      }
      
      if(end % 100000 != 0) {writefile << flush;}
      (void)ret;
      libdeflate_free_decompressor(decompressor);
      delete [] tmpout;
      delete [] snpID;
      delete [] rsID;
      delete [] chrStr;
      delete [] allele0;
      delete [] allele1;
      writefile.close();
      writefile.clear();
      fclose(fp);
      return wrap(compute_time);
    } 
    catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;
  }  
  
  
  SEXP glmm_score_bgen11(SEXP res_in, SEXP P_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
    try{
      
      Rcpp::NumericMatrix P_r(P_in);
      Rcpp::NumericVector res_r(res_in);
      size_t n_P = P_r.nrow(), k_P = P_r.ncol();
      arma::mat P(P_r.begin(), n_P, k_P, false);
      arma::vec res(res_r.begin(), res_r.size(), false);
      const char center = Rcpp::as<char>(center_in);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const char miss_method = Rcpp::as<char>(miss_method_in);
      
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      string outfile  = Rcpp::as<string>(outfile_in);
      
      const size_t npb = Rcpp::as<size_t>(nperbatch_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n), snp_skip = zeros<uvec>(npb);
      mat G(n, npb);
      string* tmpout = new string[npb];
      vector <string> biminfo;
      double gmean, geno, gmax, gmin, num, denom, pval;
      const double tol = 1e-5;
      size_t ncount, nmiss, npbidx = 0;
      clock_t time0;
      double compute_time = 0.0;
      ofstream writefile(outfile.c_str(), ofstream::out);
      
      
        
      uint maxLA = 65536;
      uint maxLB = 65536;
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];
      uint Nbgen = Rcpp::as<uint>(nbgen_in);
      uint compression = Rcpp::as<uint>(compression_in);
      uchar* zBuf1  = (unsigned char*)malloc(6 * Nbgen);
      uint16_t* shortBuf1 = (uint16_t*)malloc(6 * Nbgen);
      struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
      
      uint begin = Rcpp::as<uint>(begin_in);
      uint end   = Rcpp::as<uint>(end_in);
      long long unsigned int byte = Rcpp::as<long long unsigned int>(pos_in);
      FILE* fp = fopen(bgenfile.c_str(), "rb");
      fseek(fp, byte, SEEK_SET);
      
      bool isMultiThread = Rcpp::as<bool>(isMultiThread_in);
    
      if (!isMultiThread){
        writefile << "SNP\tRSID\tCHR\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
      }
      
      int ret;

      for (uint m = begin; m < end; m++) {
        stringstream writeout;
        maxLA = 65536;
        maxLB = 65536;
        
        uint Nprob; ret = fread(&Nprob, 4, 1, fp); 
        if (Nprob != Nbgen) {
          Rcout << "Error reading bgen file: Number of samples with genotype probabilities does not match number of samples in BGEN file. \n"; return R_NilValue;
        }
        ushort LS; ret = fread(&LS, 2, 1, fp);
        ret = fread(snpID, 1, LS, fp); snpID[LS] = '\0';
        
        ushort LR; ret = fread(&LR, 2, 1, fp);
        ret = fread(rsID, 1, LR, fp);  rsID[LR] = '\0';
        
        ushort LC; ret = fread(&LC, 2, 1, fp);
        ret = fread(chrStr, 1, LC, fp); chrStr[LC] = '\0';
        
        uint physpos; ret = fread(&physpos, 4, 1, fp);
	      string physpos_tmp = to_string(physpos);

        uint32_t LA; ret = fread(&LA, 4, 1, fp);
        if (LA > maxLA) {
          maxLA = 2 * LA;
          delete[] allele1;
          allele1 = new char[maxLA + 1];
        }
        ret = fread(allele1, 1, LA, fp); allele1[LA] = '\0';
        
        uint32_t LB; ret = fread(&LB, 4, 1, fp);
        if (LB > maxLB) {
          maxLB = 2 * LB;
          delete[] allele0;
          allele0 = new char[maxLB + 1];
        }
        ret = fread(allele0, 1, LB, fp); allele0[LB] = '\0';
        

        uLongf destLen1 = 6 * Nbgen;
        
        if (compression == 1) {
          uint cLen; ret = fread(&cLen, 4, 1, fp);
          ret = fread(zBuf1, 1, cLen, fp);
          
          if (libdeflate_zlib_decompress(decompressor, &zBuf1[0], cLen, &shortBuf1[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
			  Rcout << "Error reading bgen file: Decompressing variant block failed with libdeflate. \n"; return R_NilValue;
          }
        }
        else {
          ret = fread(zBuf1, 1, destLen1, fp);
          shortBuf1 = reinterpret_cast<uint16_t*>(zBuf1);
        }
        
        const double scale = 1.0 / 32768;
        

        gmean=0.0;
        gmax=-100.0;
        gmin=100.0;
        nmiss=0;
        ncount = 0;
        for (size_t i = 0; i < Nbgen; i++) {
          if (select[ncount] > 0) {
            
            double p11 = shortBuf1[3 * i] * scale;
            double p10 = shortBuf1[3 * i + 1] * scale;
            double p00 = shortBuf1[3 * i + 2] * scale;
            if (p11 == 0.0 && p10 == 0.0 && p00 == 0.0){
              gmiss[select[ncount]-1] = 1;
              nmiss++;
        
            } else {
              double pTot = p11 + p10 + p00;
              geno = (2 * p00 + p10) / pTot;
              gmiss[select[ncount]-1] = 0;
              g[select[ncount]-1] = geno;
              gmean += geno;
              if (geno > gmax) { gmax = geno; }
              if (geno < gmin) { gmin = geno; }
            } 
            
          }
          ncount++;
        }
        
        gmean/=(double)(n-nmiss);
        for (size_t j=0; j<n; ++j) {
          if (gmiss[j]==1) {
            g[j] = gmean;
            if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
          }
          if (center=='c') {
            g[j] -= gmean;
          }
        }
        
        gmean /= 2.0; // convert mean to allele freq
        writeout << snpID << "\t" << rsID << "\t" << chrStr << "\t" << physpos_tmp << "\t" << allele1 << "\t" << allele0 << "\t" << (n-nmiss) << "\t" << gmean << "\t";
        tmpout[npbidx] = writeout.str();
        writeout.clear();
        if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
          snp_skip[npbidx] = 1;
        } else {
          G.col(npbidx) = g;
        }
        npbidx++;
        
        
        if((m+1 == end) || (npbidx == npb)) {
          uvec snp_idx = find(snp_skip == 0);
          time0 = clock();
          vec allnum = G.cols(snp_idx).t() * res;
          mat alldenom = G.cols(snp_idx).t() * P * G.cols(snp_idx);
          compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
          for(size_t j=0; j<npbidx; ++j) {
            if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
              writefile << tmpout[j] << "NA\tNA\tNA\n";
            } else {
              uvec snp_idx2 = find(snp_idx == j);
              num = allnum[snp_idx2[0]];
              denom = alldenom(snp_idx2[0], snp_idx2[0]);
              if(denom < tol) { // G collinearity with X
                writefile << tmpout[j] << "0\t0\tNA\n";
              } else {
                pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
                writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
              }
            }
          }
          if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
        }
        if((m+1) % 100000 == 0) {writefile << flush;}
      }
      
      if(end % 100000 != 0) {writefile << flush;}
      (void)ret;
      delete [] tmpout;
      delete [] snpID;
      delete [] rsID;
      delete [] chrStr;
      delete [] allele0;
      delete [] allele1;
      writefile.close();
      writefile.clear();
      fclose(fp);
      free(zBuf1);
      free(shortBuf1);
      return wrap(compute_time);
    } 
    catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;
  }  
  
  SEXP glmm_score_bgen11_sp(SEXP res_in, SEXP Sigma_i_in, SEXP Sigma_iX_in, SEXP cov_in, SEXP bgenfile_in, SEXP outfile_in, SEXP center_in, SEXP minmaf_in, SEXP maxmaf_in, SEXP missrate_in, SEXP miss_method_in, SEXP nperbatch_in, SEXP select_in, SEXP begin_in, SEXP end_in, SEXP pos_in, SEXP nbgen_in, SEXP compression_in, SEXP isMultiThread_in) {
    try{
      
      Rcpp::NumericVector res_r(res_in);
      arma::vec res(res_r.begin(), res_r.size(), false);
      arma::sp_mat Sigma_i = Rcpp::as<arma::sp_mat>(Sigma_i_in);
      arma::sp_mat Sigma_iX = Rcpp::as<arma::sp_mat>(Sigma_iX_in);
      arma::sp_mat cov = Rcpp::as<arma::sp_mat>(cov_in);
      const char center = Rcpp::as<char>(center_in);
      const double minmaf = Rcpp::as<double>(minmaf_in);
      const double maxmaf = Rcpp::as<double>(maxmaf_in);
      const double missrate = Rcpp::as<double>(missrate_in);
      const char miss_method = Rcpp::as<char>(miss_method_in);
      
      string bgenfile = Rcpp::as<string>(bgenfile_in);
      string outfile  = Rcpp::as<string>(outfile_in);
      
      const size_t npb = Rcpp::as<size_t>(nperbatch_in);
      Rcpp::IntegerVector select(select_in);
      string line, snp;
      size_t n = res.n_elem;
      vec g(n);
      uvec gmiss(n), snp_skip = zeros<uvec>(npb);
      mat G(n, npb);
      string* tmpout = new string[npb];
      vector <string> biminfo;
      double gmean, geno, gmax, gmin, num, denom, pval;
      const double tol = 1e-5;
      size_t ncount, nmiss, npbidx = 0;
      clock_t time0;
      double compute_time = 0.0;
      ofstream writefile(outfile.c_str(), ofstream::out);
      
      
      
      uint maxLA = 65536;
      uint maxLB = 65536;
      char* snpID   = new char[maxLA + 1];
      char* rsID    = new char[maxLA + 1];
      char* chrStr  = new char[maxLA + 1];
      char* allele1 = new char[maxLA + 1];
      char* allele0 = new char[maxLA + 1];
      uint Nbgen = Rcpp::as<uint>(nbgen_in);
      uint compression = Rcpp::as<uint>(compression_in);
      uchar* zBuf1  = (unsigned char*)malloc(6 * Nbgen);
      uint16_t* shortBuf1 = (uint16_t*)malloc(6 * Nbgen);
      struct libdeflate_decompressor* decompressor = libdeflate_alloc_decompressor();
      
      uint begin = Rcpp::as<uint>(begin_in);
      uint end   = Rcpp::as<uint>(end_in);
      long long unsigned int byte = Rcpp::as<long long unsigned int>(pos_in);
      FILE* fp = fopen(bgenfile.c_str(), "rb");
      fseek(fp, byte, SEEK_SET);
      
      bool isMultiThread = Rcpp::as<bool>(isMultiThread_in);
      
      if (!isMultiThread){
        writefile << "SNP\tRSID\tCHR\tPOS\tA1\tA2\tN\tAF\tSCORE\tVAR\tPVAL\n";
      }
     
      int ret; 
      for (uint m = begin; m < end; m++) {
        stringstream writeout;
        maxLA = 65536;
        maxLB = 65536;
        
        uint Nprob; ret = fread(&Nprob, 4, 1, fp); 
        if (Nprob != Nbgen) {
          Rcout << "Error reading bgen file: Number of samples with genotype probabilities does not match number of samples in BGEN file. \n"; return R_NilValue;
        }
        ushort LS; ret = fread(&LS, 2, 1, fp);
        ret = fread(snpID, 1, LS, fp); snpID[LS] = '\0';
        
        ushort LR; ret = fread(&LR, 2, 1, fp);
        ret = fread(rsID, 1, LR, fp);  rsID[LR] = '\0';
        
        ushort LC; ret = fread(&LC, 2, 1, fp);
        ret = fread(chrStr, 1, LC, fp); chrStr[LC] = '\0';
        
        uint physpos; ret = fread(&physpos, 4, 1, fp);
	string physpos_tmp = to_string(physpos);

        uint32_t LA; ret = fread(&LA, 4, 1, fp);
        if (LA > maxLA) {
          maxLA = 2 * LA;
          delete[] allele1;
          allele1 = new char[maxLA + 1];
        }
        ret = fread(allele1, 1, LA, fp); allele1[LA] = '\0';
        
        uint32_t LB; ret = fread(&LB, 4, 1, fp);
        if (LB > maxLB) {
          maxLB = 2 * LB;
          delete[] allele0;
          allele0 = new char[maxLB + 1];
        }
        ret = fread(allele0, 1, LB, fp); allele0[LB] = '\0';
        
        
        uLongf destLen1 = 6 * Nbgen;
        
        if (compression == 1) {
          uint cLen; ret = fread(&cLen, 4, 1, fp);
          ret = fread(zBuf1, 1, cLen, fp);
          
          if (libdeflate_zlib_decompress(decompressor, &zBuf1[0], cLen, &shortBuf1[0], destLen1, NULL) != LIBDEFLATE_SUCCESS) {
	    Rcout << "Error reading bgen file: Decompressing variant block failed with libdeflate. \n"; return R_NilValue;
          }
        }
        else {
          ret = fread(zBuf1, 1, destLen1, fp);
          shortBuf1 = reinterpret_cast<uint16_t*>(zBuf1);
        }
        
        const double scale = 1.0 / 32768;
        
        
        gmean=0.0;
        gmax=-100.0;
        gmin=100.0;
        nmiss=0;
        ncount = 0;
        for (size_t i = 0; i < Nbgen; i++) {
          if (select[ncount] > 0) {
            
            double p11 = shortBuf1[3 * i] * scale;
            double p10 = shortBuf1[3 * i + 1] * scale;
            double p00 = shortBuf1[3 * i + 2] * scale;
            if (p11 == 0.0 && p10 == 0.0 && p00 == 0.0){
              gmiss[select[ncount]-1] = 1;
              nmiss++;
              
            } else {
              double pTot = p11 + p10 + p00;
              geno = (2 * p00 + p10) / pTot;
              gmiss[select[ncount]-1] = 0;
              g[select[ncount]-1] = geno;
              gmean += geno;
              if (geno > gmax) { gmax = geno; }
              if (geno < gmin) { gmin = geno; }
            } 
            
          }
          ncount++;
        }
        
        gmean/=(double)(n-nmiss);
        for (size_t j=0; j<n; ++j) {
          if (gmiss[j]==1) {
            g[j] = gmean;
            if (center=='n' && miss_method=='o') {g[j] = 0.0;} // remove missing genotypes
          }
          if (center=='c') {
            g[j] -= gmean;
          }
        }
        
        gmean /= 2.0; // convert mean to allele freq
        writeout << snpID << "\t" << rsID << "\t" << chrStr << "\t" << physpos_tmp << "\t" << allele1 << "\t" << allele0 << "\t" << (n-nmiss) << "\t" << gmean << "\t";
        tmpout[npbidx] = writeout.str();
        writeout.clear();
        if((gmax-gmin<tol) || ((double)nmiss/n>missrate) || ((gmean<minmaf || gmean>maxmaf) && (gmean<1-maxmaf || gmean>1-minmaf))) { // monomorphic, missrate, MAF
          snp_skip[npbidx] = 1;
        } else {
          G.col(npbidx) = g;
        }
        npbidx++;
        
        
        if((m+1 == end) || (npbidx == npb)) {
          uvec snp_idx = find(snp_skip == 0);
          time0 = clock();
          vec allnum = G.cols(snp_idx).t() * res;
          sp_mat Gsp(G.cols(snp_idx));
          sp_mat XSigma_iG = Sigma_iX.t() * Gsp;
          sp_mat alldenom = Gsp.t() * Sigma_i * Gsp - XSigma_iG.t() * cov * XSigma_iG;
          compute_time += (clock()-time0)/double(CLOCKS_PER_SEC);
          for(size_t j=0; j<npbidx; ++j) {
            if(snp_skip[j] == 1) { // monomorphic, missrate, MAF
              writefile << tmpout[j] << "NA\tNA\tNA\n";
            } else {
              uvec snp_idx2 = find(snp_idx == j);
              num = allnum[snp_idx2[0]];
              denom = alldenom(snp_idx2[0], snp_idx2[0]);
              if(denom < tol) { // G collinearity with X
                writefile << tmpout[j] << "0\t0\tNA\n";
              } else {
                pval = Rf_pchisq(num*num/denom, 1.0, 0, 0);
                writefile << tmpout[j] << num << "\t" << denom << "\t" << pval << "\n";
              }
            }
          }
          if(npbidx == npb) {npbidx = 0; snp_skip.zeros();}
        }
        if((m+1) % 100000 == 0) {writefile << flush;}
      }
      
      if(end % 100000 != 0) {writefile << flush;}
      (void)ret;
      delete [] tmpout;
      delete [] snpID;
      delete [] rsID;
      delete [] chrStr;
      delete [] allele0;
      delete [] allele1;

      writefile.close();
      writefile.clear();
      fclose(fp);
      free(zBuf1);
      free(shortBuf1);
      return wrap(compute_time);
    } 
    catch( std::exception &ex ) {
      forward_exception_to_r( ex );
    } catch(...) {
      ::Rf_error( "C++ exception (unknown reason)..." );
    }
    return R_NilValue;
  }  
  
} // end of extern "C"
