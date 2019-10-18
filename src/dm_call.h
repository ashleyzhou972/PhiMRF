#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define MATHLIB_STANDALONE
#include <Rmath.h>
#include "cblas_double_metropolis.h"
#include "ms_regular_metropolis.h"
#include "cblas_negpotential.h"
#include "sparse.h"
#include "eigen.h"

void allocate_column(double *w, double **w_bycol, int N, int T);
void initialize_from_input(double *w_in, double alpha0, double eta0,
		double tau20, double *w, double *alpha, double *eta,
		double *tau2, int N, int T);
SEXP double_metropolis_cont(SEXP T_in, SEXP N_in, SEXP y_in, SEXP dim_in, SEXP val_in, SEXP row_ind_in, SEXP col_ptr_in, SEXP vars_in, SEXP bounds_alpha, SEXP bounds_eta, SEXP bounds_tau2, SEXP initials, SEXP wInitials);
