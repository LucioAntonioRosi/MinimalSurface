#pragma once

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_OPENMP
	#include <omp.h>
#endif

void inline blas_copy(const double *__restrict src, double *__restrict dest,
		      size_t N)
{
	memcpy(dest, src, N * sizeof(double));
}

void inline blas_axpy(double a, const double *__restrict x,
		      double *__restrict y, size_t N)
{
#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for (size_t i = 0; i < N; ++i) {
		y[i] += a * x[i];
	}
}

void inline blas_axpby(double a, const double *__restrict x, double b,
		       double *__restrict y, size_t N)
{
#ifdef USE_OPENMP
	#pragma omp parallel for
#endif
	for (size_t i = 0; i < N; ++i) {
		y[i] = a * x[i] + b * y[i];
	}
}

double inline blas_dot_old(const double *x, const double *y, size_t N)
{
	double res = 0.0;
#ifdef USE_OPENMP
	#pragma omp parallel for reduction(+ : res)
#endif
	for (size_t i = 0; i < N; ++i) {
		res += x[i] * y[i];
	}
	return (res);
}

double inline blas_sum(const double *x, size_t N)
{
	double sum = 0;
#ifdef USE_OPENMP
	#pragma omp parallel for reduction(+ : sum)
#endif
	for (size_t i = 0; i < N; ++i) {
		sum += x[i];
	}
	return sum;
}

double inline blas_sum_in_place(double *x, size_t N)
{
	size_t step = 1;
	while (step < N) {
		size_t i = step;
		while (i < N) {
			x[i - step] += x[i];
			i += 2 * step;
		}
		step *= 2;
	}
	return x[0];
}

double inline blas_dot_new(const double *x, const double *y, size_t N)
{
	double *buf = (double *)malloc(N * sizeof(double));
	for (size_t i = 0; i < N; ++i) {
		buf[i] = x[i] * y[i];
	}
	double res = blas_sum_in_place(buf, N);
	free(buf);
	return res;
}

double inline blas_dot(const double *x, const double *y, size_t N)
{
	return blas_dot_old(x, y, N);
}
