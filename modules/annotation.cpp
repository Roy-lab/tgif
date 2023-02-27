#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <vector>
#include "annotation.h"
#include "utils.h"

Annotations::Annotations(int nTasks, vector<string>& taskNames) {
	states = gsl_matrix_alloc((int) pow(2, nTasks), nTasks);
	utils::get_all_binary_strings_in_N_bits(states, nTasks);
	int nStates = states->size1;
	for (int i = 0; i < nStates; i++) {
		gsl_vector_view row = gsl_matrix_row(states, i);
		double sum = gsl_vector_sum(&row.vector);
		if (sum == 0) {
			stateNames.push_back("non-boundary");
		} else if (sum == 1) {
			int idx = gsl_vector_max_index(&row.vector);
			stateNames.push_back("boundary specific to "+taskNames[idx]); 
		} else if (sum == nTasks) {
			stateNames.push_back("common boundary");
		} else {
			string name("boundary in ");
			bool first = true;
			for (int j = 0; j < nTasks; j++) {
				if (gsl_matrix_get(states, i, j) == 1) {
					if (first) {
						name += taskNames[j];
						first = false;
					} else {
						name += "+" + taskNames[j];
					}
				}
			}
			stateNames.push_back(name);
		}
	}
}

Annotations::~Annotations() {
	gsl_matrix_free(states);
}

int Annotations::assign_states(gsl_matrix* X) {
	int N = X->size1; 
	int S = states ->size1;
	for (int i = 0; i < N; i++) {
		gsl_vector_view row = gsl_matrix_row(X, i);
		for (int j = 0; j < S; j++) {
			gsl_vector_view state = gsl_matrix_row(states, j);
			if (gsl_vector_equal(&row.vector, &state.vector)) {
				stateAssignments.push_back(stateNames[j]);
			}
		}
	}
	return 0;
}

int Annotations::print_annotations(string outputFile, string chro, vector<int>& coord, int binSize, int N, int map[]) {
	ofstream ofs;
	ofs.open(outputFile.c_str());
	string defaultStateName = "N/A (not enough counts in region)";
	string stateName = defaultStateName;
	for (int i = 0; i < N; i++) {
		if (map[i] >= 0) {
			stateName = stateAssignments[map[i]]; 
		} else {
			stateName = defaultStateName;
		}
		ofs << chro << "\t" << coord[i] << "\t" << coord[i] + binSize << "\t" << stateName << endl;
	}
	ofs.close();
	return 0;
}
