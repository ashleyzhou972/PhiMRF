/**
 * @file ms_regular_metropolis.c
 * Regular Metropolis-Hastings algorithm for w
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <cblas.h>
#include "nmath.h"
#include "PhiMRF.h"

/**
 * Each of the parameters alpha, eta and tau2
 * has two types of posteriors:
 * 	- The multivariate normal -1
 * 	- The negpotential function -2
 * For the negpotential function type posterior,
 * there are two methods to compute the intractable constant
 * 	- The double- metropolis algorithm
 * 	- importance sampling -3
 * Only case 2 will be in this script
 **/

/**
 * Assuming random walk chain means:
 *  	- proposal distribution q is normal with mean 0
 *  	- q(current|proposed) == q(proposed|current)
 *  	- all q terms get cancelled in calculating jump probability
 **/


double r_random_walk_chain(double current_y, double var)
{
	//simulate (and return) a new value from the proposal
	//distribution (q), given current value;
	double z, y;

	z = rnorm(0.0, sqrt(var));
	y = current_y + z;

	return y;
}

double r_random_walk_chain_reflect(double current_y, double var, double lower, double upper)
{
	double z, y; 
	//printf("bounds in reflect %g, %g\n" , upper, lower);
	//printf("existing value %g\n", current_y);
	if (current_y > upper || current_y < lower)
		printf("invalid existing value %g\n", current_y);
	z = rnorm(0.0, sqrt(var));
	y = current_y + z;
	if (y > upper || y < lower)
		y = current_y - z;
	return y;
}

/**
 * ************************************
 * this function is not used under random walk chain
 **/
double d_random_walk_chain(double current_y, double proposed_y, double var)
{
	//calculate density value of the proposal distribution (q)
	//should remain the same if exchange the current_y and proposed_y
	//positions
	return dnorm(current_y-proposed_y, 0.0, sqrt(var), 1);
}
/**
 **************************************
 **/

/**
 * Assuming random walk chain means:
 *  	- proposal distribution q is normal with mean 0
 *  	- q(current|proposed) == q(proposed|current)
 *  	- all q terms get cancelled in calculating jump probability
 **/


double jump_probability(double current, double proposed, pdf target,
		double *param)
{
	//assuming random walk chain
	//Output of target should be log scaled
	double alpha, prob;
	double numerator, denominator;

	numerator = target(proposed, param);
	denominator = target(current, param);
	if (isnan(numerator)==1) {
		error("[ERROR] numerator for jump probability is NA for dm step1\n");
	}
	alpha = numerator-denominator;
	prob = exp(alpha);
	if (prob < 1.0)
		return prob;
	else
		return 1.0;
}



/**
 *
 * Another version!
 * Not multivariate in the statistical sense
 * But calculating mu as a vector
 *
 *
 *
 **/


/**
 * The log density of univariate y given w
 * updated 20180831
 * @param y is the full vector of all observations from all samples
 * @param m number of samples
 * @param N size of a single y vector (number of genes)
 * @param i gene index
 * @param w w[i]
 * This function drops the factorial!!
 * Because factorials overflow
 * (1/y!) here can be seen as a constant with regard to w in the posterior
 **/
double log_data_density_univar(double *y, int m, int N, int i, double w)
{
	double ysum = 0.0;
	for (int j = 0; j<m; j++) 
		ysum += y[i+N*j];
	return -m*exp(w) + w*ysum;
}



/**
 * Compute the vector mu using cblas in the Intel-mkl library
 *
**/
void mean_mu(int size_w, double *mu, double *w, double *neighbor_1d, double alpha, double eta)
{	
	double w_subtracted[size_w];
	for (int k = 0; k<size_w; ++k) 
		w_subtracted[k] = w[k] - alpha;
	double prod_vec[size_w]; //the product of matrix-vector multiplication;
	matrix_vector_multiplication(size_w, prod_vec, neighbor_1d, w_subtracted);
	scalar_vector_multiplication(size_w, eta, mu, prod_vec);
	scalar_vector_summation(size_w, alpha, mu);
}
/**
 * This function does NOT return the multivariate density value
 * Instead it simply takes a mu vector that's pre-computed for all i
 * and computes the UNIVARIATE density for a given i
 * The goal is to subsitute vector multiplication repeated (size_w) times
 * with matrix-vector multiplication, to speed up computation
 * @param size_w size of data vector
 * @param i position i
 * @param w_in current w observed
 * @param tau2 variance for the normal distribution
 * @param mu vector for all means of the normal distribution
**/
double log_mrf_density_vector_mu(int size_w, int i, double w_in, double tau2, double *mu)
{
	return dnorm(w_in, mu[i], sqrt(tau2), 1);
}

double log_sum_density_vector_mu(int size_w, int i, double *y, int m, double w_in, double *w, double tau2, double *mu_vec)
{
	double sum;
	sum = log_data_density_univar(y, m, size_w, i, w_in) + log_mrf_density_vector_mu(size_w, i, w_in, tau2, mu_vec);
	return sum;
}

int metropolis_for_w_vector_mu(int t, int N, double **w, double *y, int m, double var, double *neighbor_1d, double alpha, double eta, double tau2)
{
	// w is a N by T matrix;
	// y is a size N*m  vector
	// Calculating vector mu before iteration for each i;
	//initialize mu to 0 for each iteration;
	double mu[N];
	for (int k =0; k<N;k++) {
		mu[k] = 0.0;
	}
	mean_mu(N, mu, w[t], neighbor_1d, alpha, eta);	
	int jumps = 0;
	for (int i = 0; i < N; ++i) {
		double w_new_i, prob, jp;
		double numerator, denominator;
		double u = runif(0.0, 1.0);
		int jump_sig = 0;

		w_new_i = r_random_walk_chain(w[t][i], var);
		/**
		 * calculating jump probability
		 * not using the jump_probability function
		 * updated 20180523
		 **/
		numerator = log_sum_density_vector_mu(N, i, y, m,  w_new_i, w[t], tau2, mu);
		denominator = log_sum_density_vector_mu(N, i, y, m, w[t][i], w[t], tau2, mu);

		prob = exp(numerator - denominator );
		if (prob < 1.0)
			jp = prob;
		else
			jp = 1;
		//end of calculating jump probability;
		if (u < jp) {
			w[t+1][i] = w_new_i;
			jumps += 1;
			jump_sig = 1;
		}
		else 
			w[t+1][i] = w[t][i];
		/**
		if (i==999) {
			printf("mu[%d] is %g\n", i, mu[i]);
			printf("w_new[%d] is %g\n", i,  w_new_i);
			printf("w_current[%d] is %g\n", i,  w[t][i]);
			printf("numerator[%d] is %g\n", i,  numerator);
			printf("denominator[%d] is %g\n", i,  denominator);
			printf("jump alpha[%d] is %g\n", i,  jp);
			printf("jump[%d] is %d\n", i,  jump_sig);
		}
		**/
	}
	return jumps;
}

/**
 * Wrapper for matrix vector multiplication using MKL CBLAS
 * layout of matrix (rowmajor or columnmajor)
 * and triangle (upper or lower)
 * are irrelavant since matrix is symmetric
 * @param n size of matrix / vector
 * @param out output vector
 * @param matrix input matrix
 * @param vector input vector
 **/
void matrix_vector_multiplication(int n, double *out, double *matrix, double *vector)
{
      	int         lda;
      	double          alpha, beta;
      	//CBLAS_LAYOUT    layout;
      	//CBLAS_UPLO      uplo;
	alpha = 1.0;
	beta = 0.0;
	//layout = CblasRowMajor;	
	//uplo = CblasLower;
	lda = n ;
	int incx = 1.0;
	int incy = 1.0;
	//cast input matrix to double ;
	double * dmatrix;
	dmatrix = malloc(n*n*sizeof(double));
	if (dmatrix==NULL) {
		error("[ERROR] Failed to allocate memory in matrix vector multiplication\n");
	} else {
		for (int i = 0; i<n*n; ++i) {
			dmatrix[i] = (double) matrix[i];
		}
		//cblas_dsymv(layout, uplo, n, alpha, dmatrix, lda, vector, incx, beta, out, incy);
		cblas_dsymv(CblasRowMajor, CblasLower, n, alpha, dmatrix, lda, vector, incx, beta, out, incy);
		free(dmatrix); 
	}
}
void scalar_vector_multiplication(int n, double scalar, double *out, double *vector)
{
	int incx = 1.0;
	int incy = 1.0;
	cblas_daxpy(n, scalar, vector, incx, out, incy);
}

/**
 * scalar repeated n times, to form a identical vector of a;
 * plus another input vector
 * @param vector is updated;
 **/
void scalar_vector_summation(int n, double scalar, double *vector)
{

	double x[n];
	for (int i = 0; i<n ; ++i) 
		x[i] = 1.0;
	int incx = 1;
	int incy = 1;
	cblas_daxpy(n, scalar, x, incx, vector, incy);
} 
