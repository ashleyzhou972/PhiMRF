typedef struct sparse
{
	int dim;
	int nnz;
	double *val;
	int *row_ind;
	int *col_ptr;
} sparse;

void sparse_get_row(struct sparse M, double *output_row, int row_index);
void sparse_get_dense(struct sparse M, double **output_mat);


