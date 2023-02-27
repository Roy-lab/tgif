#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <list>
#include <math.h>
#include <queue>
#include <string>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include "initialization.h"

extern "C" {
	#include "random_svd/low_rank_svd_algorithms_gsl.h"
}

int init::nndsvd(gsl_matrix* A, gsl_matrix* W, gsl_matrix* H, int seed) {
	int m = A->size1;
	int n = A->size2;
	int k = W->size1;

	double a = 0; //avg 
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < n; j++) {
			a += A->data[i * A->tda + j];	
		}
	}
	a /= (m*n);
	
	
	gsl_matrix *U, *S, *V;
	
	randomized_low_rank_svd3(A, k, 3, 1, &U, &S, &V, seed);
	
	gsl_vector_view x = gsl_matrix_column(U, 0);
	gsl_vector_view y = gsl_matrix_column(V, 0);
	double sigma = sqrt(gsl_matrix_get(S, 0, 0));
	gsl_vector_scale(&x.vector, sigma);
	gsl_vector_scale(&y.vector, sigma);
	gsl_matrix_set_row(W, 0, &x.vector);
	gsl_matrix_set_row(H, 0, &y.vector);

	for (int i = 1; i < k; i++) {
		gsl_vector* xp = gsl_vector_calloc(m);
		gsl_vector* xn = gsl_vector_calloc(m);
		gsl_vector* yp = gsl_vector_calloc(n);
		gsl_vector* yn = gsl_vector_calloc(n);

		for (int j = 0; j < m; j++) {
			double val = gsl_matrix_get(U, j, i);
			if (val >= 0) {
				gsl_vector_set(xp, j, val);
			} else {	
				gsl_vector_set(xn, j, val*-1);
			}
		}		

		for (int j = 0; j < n; j++) {
			double val = gsl_matrix_get(V, j, i);
			if (val >= 0) {
				gsl_vector_set(yp, j, val);
			} else {	
				gsl_vector_set(yn, j, val*-1);
			}
		}		
	
		double xpNorm = gsl_blas_dnrm2(xp);
		double xnNorm = gsl_blas_dnrm2(xn);
		double ypNorm = gsl_blas_dnrm2(yp);
		double ynNorm = gsl_blas_dnrm2(yn);

		double mp = xpNorm * ypNorm;
		double mn = xnNorm * ynNorm;


		if (mp > mn) {
			gsl_vector_scale(xp, 1/xpNorm * sqrt(gsl_matrix_get(S,i,i)*mp));
			gsl_vector_scale(yp, 1/ypNorm * sqrt(gsl_matrix_get(S,i,i)*mp));	
			gsl_matrix_set_row(W, i, xp);
			gsl_matrix_set_row(H, i, yp);	
		} else {
			gsl_vector_scale(xn, 1/xnNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_vector_scale(yn, 1/ynNorm * sqrt(gsl_matrix_get(S,i,i)*mn));
			gsl_matrix_set_row(W, i, xn);
			gsl_matrix_set_row(H, i, yn);
		}

		gsl_vector_free(xp);
		gsl_vector_free(xn);
		gsl_vector_free(yp);
		gsl_vector_free(yn);
	}
	
	gsl_matrix_free(U);
	gsl_matrix_free(V);
	gsl_matrix_free(S);	
	
	for (int i=0; i < k; i++) {
		for (int j=0; j < m; j++) {
			double *val = &(W->data[i * W->tda + j]);
			if (*val <= 0) {
				*val = a;
			}			
		}
	}

	for (int i=0; i < k; i++) {
		for (int j=0; j < n; j++) {
			double *val = &(H->data[i * H->tda + j]);
			if (*val <= 0) {
				*val = a;
			}			
		}
	}

	return 0;
}

int init::random(gsl_matrix* W, gsl_matrix* H, int seed) {
	int k = W->size1;
	int rowCnt = W->size2;
	int colCnt = H->size2;

	const gsl_rng_type* T;
	gsl_rng* ri;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	ri = gsl_rng_alloc(T);
	gsl_rng_set(ri, seed);

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < rowCnt; j++) {
			double *val = &(W->data[i * W->tda + j]);
			*val = gsl_rng_uniform(ri);
		}
	}

	for (int i = 0; i < k; i++) {
		for (int j = 0; j < colCnt; j++) {
			double *val = &(H->data[i * H->tda + j]);
			*val = gsl_rng_uniform(ri);
		}
	}

	gsl_rng_free(ri);
	return 0;
}

