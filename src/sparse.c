#include<stdio.h>
#include<stdlib.h>


/**
 * Using structs to store CCS (compressed column storage) format sparse matrices
 * functions
 * 	get row
 * 	eigen (using MKL)
 *	inverse (solving linear system, MKL)
 *
 * Although for testing purposes M is stored in double pointer as well
 * in practice, we get the matrix directly from CCS format
 **/


typedef struct sparse
{
	int dim;
	int nnz;
	double *val;
	int *row_ind;
	int *col_ptr;
} sparse;

/**
 * @param row_index should be zero-based
 * @M sparse matrix stored in the sparse struct, make sure it's zero-based as well
 **/
void sparse_get_row(struct sparse M, double *output_row, int row_index)
{
	//iterate through each row_ind, if found  	
	if (row_index >= M.dim) 
	{
		fprintf(stderr, "[ERROR] desired row index is out of bounds\n");
		abort();
	}
	for (int i = 0; i<M.dim;++i) output_row[i] = 0;
	for (int i = 0; i<M.nnz; ++i)
	{
		if (row_index == M.row_ind[i]) 
		{
			//printf("%d\n" , i);
			for (int j = 0; j<M.dim; ++j) 
			{
				if ((i+1)> M.col_ptr[j] && (i+1)<=M.col_ptr[j+1])
				{
					output_row[j] = M.val[i];
					break;
				}
			}
		}
	}
}

/**
 * get back the dense matrix from the sparse format
 * @output_mat should be a DOUBLE pointer. Each row/column should be allocated memory;
 **/
void sparse_get_dense(struct sparse M, double **output_mat)
{
	for (int i = 0; i< M.dim; ++i)
	{
		sparse_get_row(M, output_mat[i], i);	
	}
}

