#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <vector>
#include <cmath>
#include <iostream>
#include "node.h"
#include "branch.h"
#include "utils.h"

int Branch::update() {
	gsl_matrix_memcpy(V, parent->get_V());
	for (vector<Node*>::iterator itr=children.begin(); itr != children.end(); ++itr) {
		gsl_matrix_add(V, (*itr)->get_V());
	}	
	gsl_matrix_scale(V, 1.0/(children.size()+1));
	return 0;
}

int Branch::calculate_objective(double &frobError, double &regError) {
	frobError = 0;
	gsl_matrix* temp = gsl_matrix_alloc(V->size1, V->size2);
	gsl_matrix_memcpy(temp, V);
	gsl_matrix_sub(temp, parent->get_V());
	regError = alpha * utils::get_frobenius_norm(temp);
	//error /= ((V->size1) * (V->size2));
	gsl_matrix_free(temp);
	return 0;
}

