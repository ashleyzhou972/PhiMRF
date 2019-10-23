/**
 * @file cblas_negpotential.c
 * Compute the negpotential function for the marginal distribution of w
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include <omp.h>
#include <cblas.h>
#include "nmath.h"
#include "PhiMRF.h"

/**
 * @file negpotential.c
 * @author Naihui Zhou (nzhou@iastate.edu)
 * Function to calculate negpotential function Q
 *
 **/

double negpotential(double *y, double *ystar, int size_y, double **neighbor, double alpha, double eta, double tau2)
{
	double result;

#ifdef OPENMP
#ifndef USE_NUM_THREADS
#	define USE_NUM_THREADS omp_get_max_threads()
#endif

	int num_threads = USE_NUM_THREADS;
	//int num_threads = 8;
	//printf("Number of threads used :%d\n", num_threads); 
	double sum_1[num_threads];

	for (int k = 0; k < num_threads; ++k)
		sum_1[k] = 0.0;

	omp_set_num_threads(num_threads); // Call needed threads
	# pragma omp parallel for // To run on all available threads
	for (int i = 0; i < size_y; ++i) {
		int id = omp_get_thread_num(); // Get id of thread

		sum_1[id] += H_i(y, i, ystar, size_y, neighbor, alpha, eta, tau2);
		double sum_2 = 0.0;
		for (int j = 0; j < size_y; ++j)
			sum_2 += H_ij(y, i, j, ystar, size_y, neighbor, alpha, eta, tau2);

		sum_1[id] += sum_2;
	}
	result = 0.0;
	for (int k = 0; k < num_threads; ++k)
		result += sum_1[k];

#else

	double summand_i = 0;
	double summand_ij = 0;

	for (int i = 0; i < size_y; ++i) {
		summand_i += H_i(y, i, ystar, size_y, neighbor, alpha, eta, tau2);
		for (int j = 0; j < size_y; ++j) {
			summand_ij += H_ij(y, i, j, ystar, size_y, neighbor,
					alpha, eta, tau2);
		}
	}
	result = summand_i + summand_ij;
	//Rprintf("Hi is %g\n", summand_i);
	//Rprintf("Hij is %g\n", summand_ij);
#endif
	//printf("Q is %g\n", result);
	return result;
}




double H_i(double *y, int i,  double *ystar, int size_y, double **neighbor,
		double alpha, double eta, double tau2)
{
	double mu, log_dnorm;
	double ystar_subtracted[size_y];

	for (int k = 0; k < size_y; ++k)
		ystar_subtracted[k] = ystar[k] - alpha;

	mu = alpha+eta*vector_dot_product(size_y, neighbor[i], ystar_subtracted);
	log_dnorm = dnorm(y[i], mu, sqrt(tau2), 1) - dnorm(ystar[i], mu, sqrt(tau2), 1);

	//1 for log = TRUE

	return log_dnorm;
}


double H_ij(double *y, int i, int j,  double *ystar, int size_y,
		double **neighbor, double alpha, double eta, double tau2)
{
	if (neighbor[i][j] == 0.0)
		return 0.0;
	else {
		double mu1, mu2, log_dnorm;
		double ystar_subtracted[size_y];

		for (int k = 0; k < size_y; ++k)
			ystar_subtracted[k] = ystar[k] - alpha;

		mu2 = alpha + eta * vector_dot_product(size_y, neighbor[i], ystar_subtracted);
		ystar_subtracted[j] = y[j] - alpha;
		mu1 = alpha + eta * vector_dot_product(size_y, neighbor[i], ystar_subtracted);

		double ldnorm11, ldnorm12, ldnorm21, ldnorm22;

		ldnorm11 = dnorm(y[i], mu1, sqrt(tau2), 1);
		ldnorm12 = dnorm(ystar[i], mu2, sqrt(tau2), 1);
		ldnorm21 = dnorm(ystar[i], mu1, sqrt(tau2), 1);
		ldnorm22 = dnorm(y[i], mu2, sqrt(tau2), 1);
		log_dnorm = ldnorm11 + ldnorm12-(ldnorm21+ldnorm22);

		return log_dnorm;
	}
}
double vector_dot_product(int n, double *vector1, double *vector2) {
	int incx = 1;
	int incy = 1;
	//convert int vector to double;
	double *dvector1;
	dvector1 = malloc(n*sizeof(double));
	for (int i = 0;i<n;++i) {
		dvector1[i] = (double) vector1[i];
	}
	double res;
	res = cblas_ddot(n, dvector1, incx, vector2, incy);
	free(dvector1);
	return res;
}
