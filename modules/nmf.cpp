#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <cmath>
#include <iostream>
#include "initialization.h"
#include "nmf.h"
#include "utils.h"

NMF::NMF(int k, init_method initMethod, int maxIter, int seed, bool verb, double termTol) {
	n_components = k;
	init = initMethod;
	max_iter = maxIter;
	random_state = seed;
	verbose = verb;
	tol = termTol;
	list<double> reconstruction_err_;	
}

NMF::~NMF() {};

int NMF::update_kth_block_of_U(int k){
	gsl_vector_view u_k = gsl_matrix_row(U, k);
	gsl_vector_view v_k = gsl_matrix_row(V, k);
	double v_norm = pow(gsl_blas_dnrm2(&v_k.vector), 2);
	gsl_blas_dgemv(CblasNoTrans, 1/v_norm, R, &v_k.vector, 0, &u_k.vector);
	int n = R->size1;
	for (int i = 0; i < n; i++) {
		double *val = &((&u_k.vector)->data[i]);
		if (*val < 0) {
			*val = 0;
		}
	}
	if (gsl_vector_isnull(&u_k.vector)) {
		gsl_vector_set_all(&u_k.vector, 0.01);
	}
	return 0;
}

int NMF::update_kth_block_of_V(int k){
	gsl_vector_view u_k = gsl_matrix_row(U, k);
	gsl_vector_view v_k = gsl_matrix_row(V, k);
	double u_norm = pow(gsl_blas_dnrm2(&u_k.vector), 2);
	gsl_blas_dgemv(CblasTrans, 1/u_norm, R, &u_k.vector, 0, &v_k.vector);
	int m = R->size2;
	for (int i = 0; i < m; i++) {
		double *val = &((&v_k.vector)->data[i]);
		if (*val < 0) {
			*val = 0;
		}
	}
	if (gsl_vector_isnull(&v_k.vector)) {
		gsl_vector_set_all(&v_k.vector, 0.01);
	}
return 0;
}

int NMF::update() {
	for (int k = 0; k < n_components; k++) {
		gsl_vector_view u_k = gsl_matrix_row(U, k);
		gsl_vector_view v_k = gsl_matrix_row(V, k);
		gsl_blas_dger(1, &u_k.vector, &v_k.vector, R);
		update_kth_block_of_U(k);
		update_kth_block_of_V(k);
		gsl_blas_dger(-1, &u_k.vector, &v_k.vector, R);
	}
  normalize_and_scale();
	return 0;
}

int NMF::normalize_and_scale() {
	for (int k = 0; k < n_components; k++) {
		gsl_vector_view u_k = gsl_matrix_row(U, k);
		gsl_vector_view v_k = gsl_matrix_row(V, k);
		double norm = gsl_blas_dnrm2(&v_k.vector);
		gsl_vector_scale(&v_k.vector, 1/norm);
		gsl_vector_scale(&u_k.vector, norm);
	}
	return 0;
}

double NMF::calculate_objective() {
	gsl_matrix_memcpy(R, X);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, -1, U, V, 1, R);
	double error = utils::get_frobenius_norm(R);
	reconstruction_err_.push_back(error);
	return error;
}

int NMF::fit(gsl_matrix* inputmat, gsl_matrix* W, gsl_matrix* H) {
	X = inputmat;
	int n = X->size1;
	int m = X->size2;
	U = W;
	V = H;
	R = gsl_matrix_alloc(n,m);

	if ((U->size1 != n_components) || (V->size1 != n_components)) {
		cout << "The first dimension of U and V (i.e. their number of rows) should equal the number of components specified when instantiating NMF." << endl;
		return 1;
	} 

	if (verbose) {
		cout << "Initializing..." << endl;
	}	
	if (init == random_init) {
		init::random(U, V, random_state);
	} else {
		init::nndsvd(X, U, V, random_state);
	}

	double old_error = calculate_objective();
	double old_slope;
	for (int n_iter =0; n_iter < max_iter; n_iter++){
		update();
		double error = calculate_objective();
		double slope = old_error - error;
		if (verbose) {
			cout << "Itr " << n_iter+1 << " error = " << error << ", slope = " << slope << endl;
		}
		if (slope < tol) {
			if (verbose) {
				cout << "Converged at iteration " << n_iter+1 << endl;	
			}
			break;
		} else {
			old_error = error;
			old_slope = slope;
		}
	}
	gsl_matrix_free(R);

	/*
	for (int i = 0; i < n_components; i++) {
		for (int j = 0; j < 10; j++) {
			cout << U->data[i * U->tda + j] << ", ";
		}
		cout << endl;
	}
	*/

	return 0;
}

