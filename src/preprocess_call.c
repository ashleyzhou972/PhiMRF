/**
 * @file preprocess_call.c
 * C interface with R's .Call() function for matrix transformation and eigenvalue computation
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

#include <Rinternals.h>
#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sparse.h"
#include "eigen.h"

/**
 * Goal 1: transform matrix by dividing nonzero values by rowsum + colsum
 * Goal 2: find out the bounds for eta by calculating eigen values
 **/
SEXP preprocess(SEXP dim_in, SEXP val_in, SEXP row_ind_in, SEXP col_ptr_in)
{
	//declaration of variables;
	//m is number of samples;
	int dim_neighbor;
	int *dim_p;
	int num_protected = 0;
	//Check incoming SEXP types and extract data
	//############ Begin reading neighborhood matrix ##################
	if (!isInteger(dim_in))
		error("[ERROR] matrix dimension must be integer");
	else {
		dim_p = INTEGER(dim_in);
		dim_neighbor = *dim_p;
	}
	// create sparse struct
	struct sparse neighbor_sparse;
	neighbor_sparse.dim = dim_neighbor;
	if (!isReal(val_in) || !isVector(val_in))
		error("[ERROR] Neighborhood argument 1: nonzero values must be double vector");
	else {
		neighbor_sparse.val = REAL(val_in);
		neighbor_sparse.nnz = length(val_in);
	}
	if (!isInteger(row_ind_in) || !isVector(row_ind_in))
		error("[ERROR] Neighborhood argument 2: row indices must be integer vector");
	else {
		neighbor_sparse.row_ind = INTEGER(row_ind_in);
	}
	if (!isInteger(col_ptr_in) || !isVector(col_ptr_in))
		error("[ERROR] Neighborhood argument 3: column pointers must be integer vector");
	else {
		neighbor_sparse.col_ptr = INTEGER(col_ptr_in);
	}
	//check if sparse matrix imported correctly;
	//Rprintf("val[1:2]: %g, %g\n", neighbor_sparse.val[0], neighbor_sparse.val[1]);
	//Rprintf("nnz: %d\n", neighbor_sparse.nnz);
	//Rprintf("col_ptr[1:3]: %d, %d, %d\n", neighbor_sparse.col_ptr[0], neighbor_sparse.col_ptr[1], neighbor_sparse.col_ptr[2]);
	//done creating sparse struct
		
	SEXP R_Return_List;
	SEXP R_b_eta;
	SEXP R_new_val;
	PROTECT(R_b_eta = allocVector(REALSXP, 2));
	PROTECT(R_new_val = allocVector(REALSXP, neighbor_sparse.nnz));
	PROTECT(R_Return_List = allocVector(VECSXP, 3));
	num_protected += 3;

	//################### Transforming non zero values #################;
	
	double *new_val = REAL(R_new_val);
	double rowsums[dim_neighbor];
	sparse_get_row_sum(neighbor_sparse, rowsums);
	for (int i = 0; i<neighbor_sparse.nnz; ++i)
	{
		int index[2];
		sparse_get_index(neighbor_sparse, i, index);
		//checks
		//printf("row index %d, row sum %g\n", index[0], rowsums[index[0]]);
		//printf("column index %d, column sum %g\n", index[1], rowsums[index[1]]);
		//This is assuming rowsums and colsums are the same
		//assuming symmetric matrix!!
		new_val[i] = neighbor_sparse.val[i]/(rowsums[index[0]]+rowsums[index[1]]);
		//printf("new val %g\n", new_val[i]);
	}
	SET_VECTOR_ELT(R_Return_List, 0, R_new_val);
	neighbor_sparse.val = new_val;
	//checks	
	//Rprintf("newval[1:2]: %g, %g\n", neighbor_sparse.val[0], neighbor_sparse.val[1]);
	//Rprintf("Done transforming the matrix\n");
	//################# End transforming non zero values ##############;	
	//################# Read transformed neighbor into dense format ###;
	double **neighbor_2d;
	//Rprintf("dim_neighbor is %d\n", dim_neighbor);
	neighbor_2d = malloc(dim_neighbor*dim_neighbor*sizeof(double *));
	if (neighbor_2d==NULL) {
		error("[ERROR] Failed to allocate memory for neighbor_2d\n");
	}
	//Rprintf("done allocating neighbor_2d\n");
	for (int i = 0; i < dim_neighbor; ++i) {
		neighbor_2d[i] = malloc(dim_neighbor*sizeof(double));
		if (neighbor_2d[i]==NULL) {
			error("[ERROR] Failed to allocated memory for neighbor_2d[%d]\n" , i);
		}
	}
	// get dense neighbor_2d
	sparse_get_dense(neighbor_sparse, neighbor_2d);
	//check if dense matrix converted correctly
	//for (int i = 0; i<10;++i) Rprintf("%g ", neighbor_2d[0][i]);
	//Rprintf("\n");
	//get dense neighbor_1d
	//double neighbor_1d[dim_neighbor*dim_neighbor];
	//sparse_get_dense_byrow(neighbor_sparse, neighbor_1d);
	double *neighbor_1d;
	neighbor_1d = malloc(dim_neighbor*dim_neighbor*sizeof(double));
	if (neighbor_1d ==NULL) {
		fprintf(stderr, "[ERROR] Failed to allocate memory for neighbor_1d\n");
	}
	for (int i = 0; i < dim_neighbor; ++i) {
		for (int j = 0; j < dim_neighbor; ++j)  {
			neighbor_1d[j*dim_neighbor+i] = neighbor_2d[i][j];
		}
	}
	//free neighbor_2d
	for (int i = 0; i< dim_neighbor; ++i) free(neighbor_2d[i]);
	free(neighbor_2d);

	//############## End processing dense eighborhood matrix ###################

	
	//###################### Calculating eigen values ###################
	double eigens[2];
	double *b_eta = REAL(R_b_eta);
	get_max_min_eigenvalues(dim_neighbor, neighbor_1d, eigens);
	if (eigens[0]<=0 && eigens[1]>=0) {
		b_eta[0] = 1/eigens[0];
		b_eta[1] = 1/eigens[1];
	}
	else {
		b_eta[0] = 1/eigens[1];
		b_eta[1] = 1/eigens[0];
	}
	SET_VECTOR_ELT(R_Return_List, 1, R_b_eta);
	Rprintf("Bounds for eta is (%g , %g)\n", b_eta[0], b_eta[1]);
	//Rprintf("Done calculating eigen values\n");
	//################### End calculating eigen values #################;
	Rprintf("num protected %d\n", num_protected);
	UNPROTECT(num_protected);
	//free(neighbor_1d);
	return R_Return_List;
}

