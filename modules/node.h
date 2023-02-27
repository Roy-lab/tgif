#include <gsl/gsl_matrix.h>
#include <iostream>
#include <list>
#include <string>
#include <fstream>
#include "io.h"
#ifndef _node_
#define _node_
using namespace std;

class Node {
	public:
		Node(int k, double a, string name): n_components(k), alpha(a), alias(name) {};
		Node(gsl_matrix* initV, int k, double a, string name)
			: n_components(k),alpha(a),alias(name) 
		{
			V = initV;
		};
		~Node() {};

		virtual int update() = 0;
		virtual int calculate_objective(double&, double&) = 0;

		int get_n_components () {return n_components;};
		void set_n_components(int k) {n_components = k;};
		double get_alpha() {return alpha;};
		void set_alpha(double a) {alpha = a;};
 
		bool is_leaf(){return isleaf;};
		string get_alias(){return alias;};

		gsl_matrix* get_V() {return V;};
		void set_initial_V(gsl_matrix* initV) {V = initV;};
	
		virtual Node* get_parent() {return parent;};
		virtual void set_parent(Node* p) {parent=p;};
		virtual vector<Node*>* get_children() {return &children;};
		virtual void add_child(Node* child) {children.push_back(child);}; 

		virtual int write_factors_to_file(string prefix) {
			gsl_matrix* Vtrans = gsl_matrix_alloc(V->size2, V->size1);
			gsl_matrix_transpose_memcpy(Vtrans, V);
			io::write_dense_matrix(prefix+alias+"_V.txt", Vtrans);
			gsl_matrix_free(Vtrans);
			return 0;
		};

	protected:
		string alias;
		int n_components;
		int alpha;
		Node* parent;
		vector<Node*> children;
		gsl_matrix* V;
		bool isleaf;
};
#endif
