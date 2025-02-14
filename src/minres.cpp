#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "minres.h"
#include "tiny_blas.h"
#include "array.h"

size_t minres_solve(const Matrix &A, const double *__restrict b,
				double *__restrict x, double *__restrict r,
				double *__restrict p, double *__restrict Ap,
				double *rel_error, double tol, int max_iter,
				bool inited)
{
	size_t N = A.rows;
	assert(A.rows == A.cols);

    TArray<double> p1(N,0.0), Ap1(N,0.0), 
    p2(N), Ap2(N);

	if (!inited) {
		/* r_0 = b - Ax_0 */
		A.mvp(x, r);
		blas_axpby(1, b, -1, r, N);
		/* p_0 = r_0 */
		blas_copy(r, p, N);
        A.mvp(p, Ap);
        blas_copy(p, p1.data, N);
        blas_copy(Ap, Ap1.data, N);
	}

	double r2 = blas_dot(r, r, N);
	*rel_error = r2;

	int iter = 0;
	
    while((iter < max_iter) && (*rel_error > tol*tol) && blas_dot(Ap1.data, Ap1.data, N) < 1e100){
        blas_copy(p1.data, p2.data, N);
        blas_copy(Ap1.data, Ap2.data, N);
        blas_copy(p, p1.data, N);
        blas_copy(Ap, Ap1.data, N);
    
        double alpha = blas_dot(r, Ap1.data, N)/blas_dot(Ap1.data, Ap1.data, N);
        blas_axpy(alpha, p1.data, x, N);
        blas_axpy(-alpha, Ap1.data, r, N);

        blas_copy(Ap1.data, p, N);
        A.mvp(Ap1.data, Ap);
        double beta = blas_dot(Ap, Ap1.data, N)/blas_dot(Ap1.data, Ap1.data, N);
        blas_axpby(-beta, p1.data, 1, p, N);
        blas_axpby(-beta, Ap1.data, 1, Ap, N);

        if (iter > 0){
            beta = blas_dot(Ap, Ap2.data, N)/blas_dot(Ap2.data, Ap2.data, N);
            blas_axpby(-beta, p2.data, 1, p, N);
            blas_axpby(-beta, Ap2.data, 1, Ap, N);
        }

        *rel_error = sqrt(blas_dot(r, r, N));
        iter++;
    }
    *rel_error = sqrt(*rel_error);
	return iter;
}


