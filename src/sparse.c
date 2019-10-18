/**
 * @file sparse.c
 * sparse matrix related operations
 * @author Naihui Zhou {nzhou@iastate.edu}
 **/

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

/**
 * get back the dense matrix from the sparse format
 * @output_mat should be a single pointer, collapsed by row;
 **/
void sparse_get_dense_byrow(struct sparse M, double *output_mat)
{
	for (int i = 0; i< M.dim; ++i)
	{
		sparse_get_row(M, output_mat+i*M.dim, i);	
	}
}
/**
 * get the row/column index of the k^th non-zero value
 **/
void sparse_get_index(struct sparse M, int k, int *out_index)
{
	int out;
	for (int j = 0; j< M.dim; ++j)
	{
		if (M.col_ptr[j]<=k && M.col_ptr[j+1]>k) 
		{
			out = j;
			break;
		}
	}
	out_index[0] = M.row_ind[k];
	out_index[1] = out;
}

/**
 * get the rowsum of the matrix by only adding the nonzero values
 **/
void sparse_get_row_sum(struct sparse M, double *out_row_sums)
{
	for (int k = 0; k<M.dim; ++k)
	{
		out_row_sums[k] = 0;
		for (int i = 0; i< M.nnz; ++i)
		{
			if (M.row_ind[i]==k) out_row_sums[k] += M.val[i];
		}
	}
}


