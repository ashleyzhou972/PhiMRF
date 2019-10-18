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

void sparse_get_dense_byrow(struct sparse M, double *output_mat);
void sparse_get_index(struct sparse M, int k, int *out_index);

void sparse_get_row_sum(struct sparse M, double *out_row_sums);

