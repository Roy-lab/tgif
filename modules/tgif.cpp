#include <gsl/gsl_matrix.h>
#include <list>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "node.h"
#include "nmf.h"
#include "leaf.h"
#include "branch.h"
#include "root.h"
#include "io.h"
#include "utils.h"
#include "tgif.h"

int TGIF::make_tree(vector<int>& parentIds, vector<string>& aliases, int nLeaves) { 
	int nNodes = parentIds.size();
	for (int i = 0; i < nLeaves; i++) {
		Leaf* leaf = new Leaf(n_components, alpha, aliases[i]);
		tree.push_back(leaf); 
	}
	for (int i = nLeaves; i < (nNodes-1); i++) {
		Branch* branch = new Branch(n_components, alpha, aliases[i]);
		tree.push_back(branch);
	}
	Root* root = new Root(n_components, alpha, aliases[nNodes-1]);
	tree.push_back(root);
	for (int i =0; i < nNodes; i++) {
		int parentId = parentIds[i]-1;
		if (parentId > 0) {
			Node* child = tree[i];
			Node* parent = tree[parentId]; 
			child->set_parent(parent);	
			parent->add_child(child);
		}
	}
	return 0;
} 

int TGIF::init(vector<gsl_matrix*>& input, int initTaskIdx) {
	int nLeaves = input.size();
	int nNodes = tree.size();
	int n = input[initTaskIdx] -> size1;
	int m = input[initTaskIdx] -> size2;
	gsl_matrix* U1 = gsl_matrix_calloc(n_components, n);
	gsl_matrix* V1 = gsl_matrix_calloc(n_components, m);
	NMF nmf = NMF(n_components, nndsvd_init, 100, random_state, false, 1);
	nmf.fit(input[initTaskIdx], U1, V1); 
	for (int i = 0; i < nNodes; i++) {
		gsl_matrix* V = gsl_matrix_calloc(n_components, m);
		gsl_matrix_memcpy(V, V1);
		tree[i]->set_initial_V(V);
		if (i < nLeaves) {
			((Leaf*)tree[i])->set_X(input[i]);
			gsl_matrix* U = gsl_matrix_calloc(n_components, n);
			gsl_matrix_memcpy(U, U1);
			((Leaf*)tree[i])->set_initial_U(U);
		}
	}
	gsl_matrix_free(U1);
	gsl_matrix_free(V1);
	return 0;
}

int TGIF::fit(vector<gsl_matrix*>& input) {
	TGIF::fit(input, 0);
	return 0;
}

int TGIF::fit(vector<gsl_matrix*>& input, int initTaskIdx) {
	init(input, initTaskIdx);

	double old_error = 0;
	for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			double frobError, regError;
			(*node)->calculate_objective(frobError, regError);
			old_error += (frobError + regError);	
	}
	double old_slope;

	ofstream ofs;
	if (debug) {
		string fname = output_prefix + "track.log";
		ofs.open(fname.c_str());
		ofs << "Itr\t";
		for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			if ((*node)->get_parent() != NULL) {
				ofs << (*node)->get_alias() << "\t";
				if((*node)->is_leaf()) {
					ofs << (*node)->get_alias() << "-parent\t";
				}
			}
		}
		ofs << "total\tslope" << endl;
	}
	
	for (int n_iter =0; n_iter < max_iter; n_iter++){
		for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			(*node)->update();
		}
		double error = 0;
		if (debug) {
			ofs << n_iter+1 << "\t";
		}
		for(vector<Node*>::iterator node=tree.begin(); node!=tree.end(); ++node) {
			double frobError, regError;
			(*node)->calculate_objective(frobError, regError);
			error += (frobError + regError);
			if (debug && ((*node)->get_parent() != NULL)) {
				if ((*node)->is_leaf()) {
					ofs << frobError << "\t" << regError << "\t";
				} else {
					ofs << regError << "\t";
				}
			}
		}
		reconstruction_err_.push_back(error);
		double slope = old_error - error;
		if (verbose) {
			cout << "Itr " << n_iter +1 << " error = " << error << ", slope = " << slope << endl;
		}
		if (debug) {
			ofs << error << "\t" << slope << endl;
		}
		if (slope < tol) {
			if (verbose) {
				cout << "Converged at iteration " << n_iter << endl;	
			}
			break;
		} else {
			old_error = error;
			old_slope = slope;
		}
	}
	if (debug) {
		ofs.close();
	}
	return 0;
}

int TGIF::print_factors() {
	for (vector<Node*>::iterator itr=tree.begin(); itr!=tree.end(); ++itr) {
		(*itr)->write_factors_to_file(output_prefix);
	}
	return 0;
}

int TGIF::clean() {
	for (vector<Node*>::iterator itr=tree.begin(); itr!=tree.end(); ++itr) {
		gsl_matrix_free((*itr)->get_V());
		if ((*itr)->is_leaf()) {
			gsl_matrix_free(((Leaf*)(*itr))->get_U());
		}
	}
	int l = tree.size();
	for (int i = 0; i < l; i++) {
		delete tree[i];
	}
	return 0;
}
