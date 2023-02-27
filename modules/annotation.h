#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_randist.h>
#include <math.h>
#include <iostream>
#include <float.h>
#include <vector>
#ifndef _annot
#define _annot_
using namespace std;

class Annotations {
	public:
		Annotations(int, vector<string>&); 
		~Annotations();

		vector<string> get_state_names() {return stateNames;};
		gsl_matrix* get_states() {return states;};
		int assign_states(gsl_matrix*);
		int print_annotations(string, string, vector<int>&, int, int, int[]);

	private:
		gsl_matrix* states;
		vector<string> stateNames;
		vector<string> stateAssignments;
};
#endif
