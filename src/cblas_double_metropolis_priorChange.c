/**
 * @file double_metropolis.c
 * Initial attempt at the double metropolis algorithm
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include "Rmath.h"
#include "cblas_negpotential.h"
#include "ms_regular_metropolis.h"
#include "cblas_double_metropolis.h"


double dm_step1(double theta0, pdf target, double var, double *target_params, double *bounds_theta)
{
	//in double-metropolis step 1,;
	//target_pdf is the prior distribution of theta;
	//target_pdf should return log scale
	//target_params should be a vector of parameters to the target distribution
	//updated 20190301:
	//bounds_theta should be a vector of the lower and upper bounds
	//target_params and bounds_theta could be the same in an uniform distribution
	double theta_new, jp, alpha;
	double u;

	//Updated 20190308: Do not use reflect, it changes the proposed distribution
	//theta_new = r_random_walk_chain_reflect(theta0, var, bounds_theta[0], bounds_theta[1]);
	theta_new = r_random_walk_chain(theta0, var);
	//printf("theta new from dm step1 is: %g\n", theta_new);
	jp = jump_probability(theta0, theta_new, target, target_params);
	//printf("numerator from dm_step1 is %g\n", target(theta_new, target_params));
	///printf("denominator from dm_step1 is %g\n", target(theta0, target_params));
	//printf("jump probability is %g\n", jp);
	if (jp < 1.0)
		alpha = jp;
	else
		alpha = 1.0;
	u = runif(0.0, 1.0);
	if (u < alpha)
	{
		//printf("accepted\n");
		return theta_new;
	}
	else
	{
		//printf("not accepted\n");
		return theta0;
	}
}

double dm_step1_indep(double theta0, pdf target, double var, double *target_params, double *bounds_theta)
{
	//in double-metropolis step 1,;
	//target_pdf is the prior distribution of theta;
	//target_pdf should return log scale
	//target_params should be a vector of parameters to the target distribution
	//updated 20190301:
	//bounds_theta should be a vector of the lower and upper bounds
	//target_params and bounds_theta could be the same in an uniform distribution
	//updated 20190404:
	//this function uses independence proposal for new theta
	//only used in eta
	double theta_new, jp, alpha;
	double u;

	theta_new = runif(bounds_theta[0], bounds_theta[1]);
	//printf("theta new from dm step1 is: %g\n", theta_new);
	jp = jump_probability(theta0, theta_new, target, target_params);
	if (jp < 1.0)
		alpha = jp;
	else
		alpha = 1.0;
	u = runif(0.0, 1.0);
	if (u < alpha)
	{
		//printf("accepted\n");
		return theta_new;
	}
	else
	{
		//printf("not accepted\n");
		return theta0;
	}
}
int dm_step2(int size_x, double *x, double theta_new, double theta_current,
		double **neighbor, negp negpotential_theta,
		auxiliary auxiliary_y_gibbs_theta, double *par_neg,
		double *par_auxi)
{
	//target distribution here is the original
	//data distribution f(x|theta) in the Bayes rule,
	//without the intractable constant;
	//aka the negpotential function;
	//The negpotential function should return in log scale;
	//First generate a auxiliary variable y;
	double y[size_x];
	double r, alpha, u;
	double numerator, denominator;

	auxiliary_y_gibbs_theta(y, x, par_auxi, neighbor);
	numerator = negpotential_theta(y, theta_current, par_neg, neighbor)
		+ negpotential_theta(x, theta_new, par_neg, neighbor);

	denominator = negpotential_theta(x, theta_current, par_neg, neighbor)
		+ negpotential_theta(y, theta_new, par_neg, neighbor);

	r = exp(numerator-denominator);
	if (r < 1)
		alpha = r;
	else
		alpha = 1;

	u = runif(0.0, 1.0);
	if (u <= alpha)
		return 1;//accepted;
	else
		return 0;
}

void auxiliary_y_gibbs(int size_x, double *x, double *y, double **neighbor,
		double alpha, double eta, double tau2)
{
	//One Gibbs round for each element in x;
	//The generating pdf is the normal with mrf mu;
	int i = 0;
	double mu_i, y_new_i;
	double y_subtracted[size_x];

	for (int k = 0; k < size_x; ++k) {
		y[k] = x[k]; //starting value is x
		y_subtracted[k] = y[k] - alpha;
	}

	//this may be redundant since y is initialized outside of this function;
	for (i = 0; i < size_x; ++i) {
		mu_i = alpha + eta*vector_dot_product(size_x, neighbor[i], y_subtracted);
		//mu_i = alpha + eta*vector_multiplication(neighbor[i], y, size_x);
		//updated 20190228: 
		//second argument to rnorm should be sd!!
		//instead of variance!!
		y_new_i = rnorm(mu_i, sqrt(tau2));
		y[i] = y_new_i;
		y_subtracted[i] = y[i] - alpha;
	}
}

/**
 * *******************************************************************
 * Above are the general functions
 * Below are wrappers for specific parameters, iterations, etc
 * ******************************************************************
 **/

/**
 * For alpha:
 **/
double prior_alpha(double alpha_par, double *other_par)
{
	//other_par should be the lower and upper bound of uniform distribution
	//First lower than upper!!!
	//other_par might not necessarily be bounds
	//if the distribution is not uniform!
	double lb = other_par[0];
	double ub = other_par[1];

	return dunif(alpha_par, lb, ub, 1);
}

double negp_alpha_t(double *y, double alpha_par, double *other_par, double **neighbor)
{
	//other_par should have N as first param
	//eta[t] as second param
	//tau2[t] as third param
	int N = (int) other_par[0];
	double eta = other_par[1];
	double tau2 = other_par[2];
	double ystar[N];

	for (int i = 0; i < N; ++i)
		ystar[i] = 0.0;

	return negpotential(y, ystar, N, neighbor, alpha_par, eta, tau2);
}


void auxi_alpha_t(double *y, double *x, double *other_par, double **neighbor)
{
	int N = (int) other_par[0];
	double alpha_new = other_par[1];//!this is the new alpha generated;
	double eta = other_par[2];
	double tau2 = other_par[3];

	auxiliary_y_gibbs(N, x, y, neighbor, alpha_new, eta, tau2);
}

int dm_step2_t_alpha(int t, double **w_bycol, double *alpha, double *eta,
		double *tau2, int N, int T, double alpha_new, double **neighbor)
{
	double par_neg[3];
	double par_auxi[4];

	par_neg[0] = N;
	par_neg[1] = eta[t];
	par_neg[2] = tau2[t];
	par_auxi[0] = N;
	par_auxi[1] = alpha_new;
	par_auxi[2] = eta[t];
	par_auxi[3] = tau2[t];
	//updated 20190101
	//	edited from w_bycol[t] to w_bycol[t+1]
	//	because the first step of updating w_bycol has already been completed!
	return dm_step2(N, w_bycol[t+1], alpha_new, alpha[t], neighbor,
			negp_alpha_t, auxi_alpha_t, par_neg, par_auxi);
}

/**
 * For eta:
 **/
double prior_eta(double eta_par, double *other_par)
{
	//other_par should be the upper and lower bound of uniform distribution
	double lb = other_par[0];
	double ub = other_par[1];

	return dunif(eta_par, lb, ub, 1);
}

double negp_eta_t(double *y, double eta_par, double *other_par, double **neighbor)
{
	//other_par should have N as first param
	//alpha[t+1] as second param
	//tau2[t] as third param
	int N = (int) other_par[0];
	double ystar[N];
	double alpha = other_par[1];
	double tau2 = other_par[2];

	for (int i = 0; i < N; ++i)
		ystar[i] = 0.0;

	return negpotential(y, ystar, N, neighbor, alpha, eta_par, tau2);
}


void auxi_eta_t(double *y, double *x, double *other_par, double **neighbor)
{
	int N = (int) other_par[0];
	double alpha = other_par[1];
	double eta_new = other_par[2];//This is the new eta from dm_step1
	double tau2 = other_par[3];

	auxiliary_y_gibbs(N, x, y, neighbor, alpha, eta_new, tau2);

}

int dm_step2_t_eta(int t, double **w_bycol, double *alpha, double *eta,
		double *tau2, int N, int T, double eta_new, double **neighbor)
{
	double par_neg[3];
	double par_auxi[4];

	par_neg[0] = N;
	par_neg[1] = alpha[t+1];
	par_neg[2] = tau2[t];
	par_auxi[0] = N;
	par_auxi[1] = alpha[t+1];
	par_auxi[2] = eta_new;
	par_auxi[3] = tau2[t];

	return dm_step2(N, w_bycol[t+1], eta_new, eta[t], neighbor, negp_eta_t,
			auxi_eta_t, par_neg, par_auxi);
}


/**
 * For tau2
 **/
/**
 * updated 20190313:
 * 	changed tau2 prior to
 *	tau \sim unif(0,A)
 *	by change-of-variables
 *	f(tau2) = 1/(2*sqrt(tau2)*A)*indicator(0<tau2<A2)
 **/
double prior_tau2(double tau2_par, double *other_par)
{
	//other_par should be the upper and lower bound of uniform distribution
	double lb = other_par[0];
	double ub = other_par[1];
	double pdf;
	if (tau2_par<0) {//log(tau2_par) will give nan; 
		pdf = -1.0/0.0; //make it negative infinity;
	}
	else {
		pdf = -log(2)-log(tau2_par)/2+dunif(sqrt(tau2_par), lb, ub, 1);
	}
	return pdf;
}

double negp_tau2_t(double *y, double tau2_par, double *other_par, double **neighbor)
{
	//other_par should have N as first param
	//alpha[t+1] as second param
	//eta[t+1] as third param
	int N = (int) other_par[0];
	double ystar[N];
	double alpha = other_par[1];
	double eta = other_par[2];

	for (int i = 0; i < N; ++i)
		ystar[i] = 0.0;

	return negpotential(y, ystar, N, neighbor, alpha, eta, tau2_par);
}


void auxi_tau2_t(double *y, double *x, double *other_par, double **neighbor)
{
	int N = (int) other_par[0];
	double alpha = other_par[1];
	double eta = other_par[2];
	double tau2_new = other_par[3];//!this is the new tau2 from dm_step1

	auxiliary_y_gibbs(N, x, y, neighbor, alpha, eta, tau2_new);
}

int dm_step2_t_tau2(int t, double **w_bycol, double *alpha, double *eta,
		double *tau2, int N, int T, double tau2_new, double **neighbor)
{
	double par_neg[3];
	double par_auxi[4];

	par_neg[0] = N;
	par_neg[1] = alpha[t+1];
	par_neg[2] = eta[t+1];
	par_auxi[0] = N;
	par_auxi[1] = alpha[t+1];
	par_auxi[2] = eta[t+1];
	par_auxi[3] = tau2_new;

	return dm_step2(N, w_bycol[t+1], tau2_new, tau2[t], neighbor, negp_tau2_t,auxi_tau2_t, par_neg, par_auxi);
}
