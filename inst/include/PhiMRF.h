#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
//#include "nmath.h"

void allocate_column(double *w, double **w_bycol, int N, int T);
void initialize_from_input(double *w_in, double alpha0, double eta0,
		double tau20, double *w, double *alpha, double *eta,
		double *tau2, int N, int T);
/**
 * SEXP functions
 **/
SEXP double_metropolis_cont(SEXP T_in, SEXP N_in, SEXP y_in, SEXP dim_in, SEXP val_in, SEXP row_ind_in, SEXP col_ptr_in, SEXP vars_in, SEXP bounds_alpha, SEXP bounds_eta, SEXP bounds_tau2, SEXP initials, SEXP wInitials);
SEXP preprocess(SEXP dim_in, SEXP val_in, SEXP row_ind_in, SEXP col_ptr_in);

/**
 * Mathlib functions
 **/
double dnorm(double x, double mu, double sigma, int give_log);
double dnorm4(double x, double mu, double sigma, int give_log);
double rnorm(double mu, double sigma);
double dunif(double x, double a, double b, int give_log);
double runif(double a, double b);

/**
 * typedef
 **/
typedef double (*pdf)(double, double *);
typedef double (*negp) (double *, double, double *, double **);
typedef void (*auxiliary)(double *, double *, double *, double **);
typedef struct sparse
{
	int dim;
	int nnz;
	double *val;
	int *row_ind;
	int *col_ptr;
} sparse;

/**
 * C functions
 */
double dm_step1(double theta0, pdf target, double var, double *target_params, double *bounds_theta);
double dm_step1_indep(double theta0, pdf target, double var, double *target_params, double *bounds_theta);
int dm_step2(int size_x, double *x, double theta_new, double theta_current,
		double **neighbor, negp negpotential,
		auxiliary auxiliary_y_gibbs_theta, double *par_neg,
		double *par_auxi);
void auxiliary_y_gibbs(int size_x, double *x, double *y, double **neighbor,
		double alpha, double eta, double tau2);
double vector_multiplication(int *array1, double *array2, int size);
double prior_alpha(double alpha_par, double *other_par);
double prior_eta(double eta_par, double *other_par);
double prior_tau2(double tau2_par, double *other_par);
double negp_alpha_t(double *y, double alpha_par, double *other_par, double **neighbor);
double negp_eta_t(double *y, double eta_par, double *other_par, double **neighbor);
double negp_tau2_t(double *y, double tau2_par, double *other_par, double **neighbor);
void auxi_alpha_t(double *y, double *x, double *other_par, double **neighbor);
void auxi_eta_t(double *y, double *x, double *other_par, double **neighbor);
void auxi_tau2_t(double *y, double *x, double *other_par, double **neighbor);
int dm_step2_t_alpha(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double alpha_new, double **neighbor);
int dm_step2_t_eta(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double eta_new, double **neighbor);
int dm_step2_t_tau2(int t, double **w_bycol, double *alpha, double *eta, double *tau2, int N, int T, double tau2_new, double **neighbor);
double negpotential(double *y, double *ystar, int size_y, double **neighbor,
		double alpha, double eta, double tau2);
double H_i(double *y, int i,  double *ystar, int size_y, double **neighbor,
		double alpha, double eta, double tau2);
double H_ij(double *y, int i, int j, double *ystar, int size_y,
		double **neighbor, double alpha, double eta, double tau2);
double vector_dot_product(int n, double *vector1, double *vector2) ;
void get_max_min_eigenvalues(int dim, double *M_1d, double *out);
double r_random_walk_chain(double current_y, double var);
double d_random_walk_chain(double current_y, double proposed_y, double var);
double jump_probability(double current, double proposed, pdf targeti, double *param);
double r_random_walk_chain_reflect(double current_y, double var, double lower, double upper);
void matrix_vector_multiplication(int n, double *out, double *matrix, double *vector);
void scalar_vector_multiplication(int n, double scalar, double *out, double *vector);
void scalar_vector_summation(int n, double scalar, double *vector);
double log_data_density_univar(double *y, int m, int N, int i, double w);
double log_mrf_density_univar(int size_w, int i, double w_in, double *w, double **neighbor, double alpha, double eta, double tau2);
double log_sum_density_univar(int size_w, int i, double y, double w_in, double *w, double **neighbor, double alpha, double eta, double tau2);
int metropolis_for_w_univar(int t, int N, double **w, double *y, double var, double **neighbor, double alpha, double eta, double tau2);
void mean_mu(int size_w, double *mu, double *w, double *neighbor_1d, double alpha, double eta);
double log_mrf_density_vector_mu(int size_w, int i, double w_in, double tau2, double *mu);
double log_sum_density_vector_mu(int size_w, int i, double *y, int m,  double w_in, double *w, double tau2, double *mu_vec);
int metropolis_for_w_vector_mu(int t, int N, double **w, double *y, int m, double var, double *neighbor_1d, double alpha, 
		double eta, double tau2);
void sparse_get_row(struct sparse M, double *output_row, int row_index);
void sparse_get_dense(struct sparse M, double **output_mat);
void sparse_get_dense_byrow(struct sparse M, double *output_mat);
void sparse_get_index(struct sparse M, int k, int *out_index);
void sparse_get_row_sum(struct sparse M, double *out_row_sums);

