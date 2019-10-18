#ifndef REGULAR_METROPOLIS_H
#define REGULAR_METROPOLIS_H

/**
 * @file regular_metropols.h
 * @author Naihui Zhou (nzhou@iastate.edu)
 *
 **/
typedef double (*pdf)(double, double *);
//pdf functions, fisrt  x, second is pointer to other parameters;


void get_max_min_eigenvalues(int dim, double *M_1d, double *out);

/**
 *
 * For regular metropolis-hastings within Gibbs algorithm
 **/
double r_random_walk_chain(double current_y, double var);
double d_random_walk_chain(double current_y, double proposed_y, double var);
double jump_probability(double current, double proposed, pdf targeti, double *param);
double r_random_walk_chain_reflect(double current_y, double var, double lower, double upper);

/**
 * Array helper functions
 **/
void matrix_vector_multiplication(int n, double *out, double *matrix, double *vector);
void scalar_vector_multiplication(int n, double scalar, double *out, double *vector);
void scalar_vector_summation(int n, double scalar, double *vector);

/**
 * regular metropolis for w
 **/

double log_data_density_univar(double *y, int m, int N, int i, double w);
double log_mrf_density_univar(int size_w, int i, double w_in, double *w, double **neighbor, double alpha, double eta, double tau2);
double log_sum_density_univar(int size_w, int i, double y, double w_in, double *w, double **neighbor, double alpha, double eta, double tau2);
int metropolis_for_w_univar(int t, int N, double **w, double *y, double var, double **neighbor, double alpha, double eta, double tau2);


void mean_mu(int size_w, double *mu, double *w, double *neighbor_1d, double alpha, double eta);
double log_mrf_density_vector_mu(int size_w, int i, double w_in, double tau2, double *mu);
double log_sum_density_vector_mu(int size_w, int i, double *y, int m,  double w_in, double *w, double tau2, double *mu_vec);
int metropolis_for_w_vector_mu(int t, int N, double **w, double *y, int m, double var, double *neighbor_1d, double alpha, double eta, double tau2);
#endif
