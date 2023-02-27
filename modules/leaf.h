#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <iostream>
#include <string>
#include "io.h"
#include "node.h"
#ifndef _leaf_
#define _leaf_
using namespace std;

class Leaf : public Node {
	public:
		Leaf(int k, double alpha, string alias): Node(k, alpha, alias) {isleaf=true;};
		Leaf(gsl_matrix* inputMat, gsl_matrix* initU, gsl_matrix* initV, 
				int k , double alpha, string alias) 
				: Node(initV, k, alpha, alias) {
			isleaf=true;
			X = inputMat;
			U = initU;
			int n = X->size1;
			int m = X->size2;
			R = gsl_matrix_alloc(n, m);
		};
		~Leaf() {gsl_matrix_free(R);};

		int update();
		int calculate_objective(double&, double&);

		gsl_matrix* get_X() {return X;};
		void set_X(gsl_matrix* input) {
			X = input;
			int n = X->size1;
			int m = X->size2;
			R = gsl_matrix_alloc(n, m);
		};
		gsl_matrix* get_U() {return U;};
		void set_initial_U(gsl_matrix* initU) {U = initU;};
		vector<Node*>* get_children() {return NULL;};
		void add_child(Node* c) {cout << "Leaf node can't have children." << endl;};

		int write_factors_to_file(string prefix) {
			gsl_matrix* Vtrans = gsl_matrix_alloc(V->size2, V->size1);
			gsl_matrix_transpose_memcpy(Vtrans, V);
			io::write_dense_matrix(prefix+alias+"_V.txt", Vtrans);
			gsl_matrix_free(Vtrans);
			gsl_matrix* Utrans = gsl_matrix_alloc(U->size2, U->size1);
			gsl_matrix_transpose_memcpy(Utrans, U);
			io::write_dense_matrix(prefix+alias+"_U.txt", Utrans);
			gsl_matrix_free(Utrans);
			return 0;
		};

	private:
		gsl_matrix* X;
		gsl_matrix* U;
		gsl_matrix* R;
		//gsl_matrix* L;
		//gsl_matrix* Linv;
		//double lambda;
		int update_kth_block_of_U(int);
		int update_kth_block_of_V(int);
		int normalize_and_scale();
};
#endif
