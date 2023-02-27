#include <gsl/gsl_matrix.h>
#include <list>
#ifndef _nmf_
#define _nmf_
using namespace std;

enum init_method {nndsvd_init, random_init};

class NMF {
	public:
		NMF(int, init_method, int, int, bool, double);
		~NMF();
		int fit(gsl_matrix*, gsl_matrix*, gsl_matrix*);
		int n_components;
		init_method init;
		int max_iter;
		int random_state;
		bool verbose;
		double tol;
		list<double> reconstruction_err_;

	private:
		gsl_matrix* X;
		gsl_matrix* U;
		gsl_matrix* V;
		gsl_matrix* R;
		int update_kth_block_of_U(int);
		int update_kth_block_of_V(int);
		int update();
		int normalize_and_scale();
		double calculate_objective();
};
#endif
