#include <gsl/gsl_matrix.h>
#include <list>
#include <vector>
#include <string>
#include "node.h"
#ifndef _tgif_
#define _tgif_
using namespace std;

class TGIF {
	public:
		TGIF(double a, int k, int nitr, double t, int seed, bool verb, bool deb, string outputPrefix) : 
				alpha(a),
				n_components(k),
				max_iter(nitr),
				tol(t),
				random_state(seed),
				verbose(verb),
				debug(deb),
				output_prefix(outputPrefix) {};

		~TGIF() {};

		int make_tree(vector<int>&, vector<string>&, int); //parent IDs, aliases, numLeaves
		int fit(vector<gsl_matrix*>&);
		int fit(vector<gsl_matrix*>&, int);
		int print_factors();
		int clean(); //frees up all data structures!

		vector<Node*>& get_tree(){return tree;};
		list<double>& get_errors(){return reconstruction_err_;};

	private:
		int n_components;
		int max_iter;
		int random_state;
		bool verbose;
		bool debug;
		double tol;
		double alpha;
		vector<Node*> tree;
		list<double> reconstruction_err_;
		string output_prefix;
		int init(vector<gsl_matrix*>&, int);
};
#endif
