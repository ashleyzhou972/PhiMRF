#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "lapacke.h"

void get_max_min_eigenvalues(int dim, double *M_1d, double *out)
{
	//printf("dimension of neighborhood matrix %d\n", dim);
	lapack_int n, il, iu, m, lda, ldz, info;
	double abstol, vl, vu;
	double *w;
	//z and ifail are not referenced in the case jobz=='N';
	double z[2];
	lapack_int ifail[2];
	w = malloc(dim*sizeof(double));
	if (w==NULL) {
		fprintf(stderr, "[ERROR] Failed to allocated memory for eigenvalues in dsyevx\n");
		exit(1);
	}
	il = 1;
	iu = dim;
	abstol = -1;
	vl = 0;
	vu = 0;
	int layout = LAPACK_ROW_MAJOR;	
	char jobz = 'N';
	char range = 'I';
	char uplo = 'L';
	n = dim;
	lda = dim;
	ldz = dim;
	double *a;
	printf("Begin copying matrix\n");
	a = malloc(dim*dim*sizeof(double));
        if (a==NULL)
	{
		fprintf(stderr, "[ERROR] Failed to allocated memory for matrix A in dsyevx\n");
		exit(1);
	}
	/**
	for (int i = 0; i< dim*dim; ++i)
	{
		for (int j = 0; j< dim; ++j) 
		{
			//printf("%g ", M[i][j]);
			a[j] = M_1d[j];
		}
	}
	**/
	memcpy(a, M_1d, dim*dim);
	printf("Done copying matrix\n");
	//for (int k = 0; k<dim*dim; ++k) printf("%g ", a[k]);
	printf("beginning lapacke_dsyevx...\n");
	info = LAPACKE_dsyevx( layout, jobz, range, uplo, n, a, lda, vl, vu, il, iu, abstol, &m, w, z, ldz, ifail);
	if (info!=0) {
		fprintf(stderr, "[ERROR] Failed to calcuate eigenvalues using dsyevs\n");
		exit(1);
	}
	//for (int k = 0; k < dim; ++k) printf("%g ", w[k]);
	//printf("\n");
	out[0] = w[0];
	out[1] = w[dim-1];
	printf("lowest eigen value is %g.\n", out[0]);
	printf("highest eigen value is %g.\n", out[1]);
	free(a);
	free(w);
	//exit(0);
}
