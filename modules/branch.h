#include <gsl/gsl_matrix.h>
#include <vector>
#include <string>
#include <fstream>
#include "node.h"
#ifndef _branch_
#define _branch_
using namespace std;

class Branch : public Node {
	public:
		Branch(int k, double alpha, string alias): Node(k, alpha, alias) {isleaf=false;};
		Branch(gsl_matrix* initV, int k, double alpha, string alias) 
			: Node(initV, k, alpha, alias) {isleaf=false;};
		~Branch() {};

		int update();
		int calculate_objective(double&, double&);

};
#endif
