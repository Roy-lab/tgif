#include <gsl/gsl_matrix.h>
#include <string>
#include <vector>
#include <fstream>
#include "node.h"
#ifndef _root_
#define _root_
using namespace std;

class Root : public Node {
	public:
		Root(int k, double alpha, string alias): Node(k, alpha, alias) {isleaf=false;};
		Root(gsl_matrix* initV, int k, double alpha, string alias) 
			: Node(initV, k, alpha, alias) {isleaf=false;};
		~Root() {};

		int update();
		int calculate_objective(double&, double&);

		Node* get_parent() {return NULL;};
		void set_parent(Node* p) {
			cout << "Root node can't have a parent." << endl;
		};
};
#endif
