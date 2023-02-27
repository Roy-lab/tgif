#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <vector>
#ifndef _utils_
#define _utils_
using namespace std;

namespace utils
{
	double get_frobenius_norm(gsl_matrix*);
	int get_neighborhood_graph(gsl_matrix*, int);
	int get_laplacian(gsl_matrix*, gsl_matrix*);
	int get_inverse(gsl_matrix*);
	int sample_by_distance(gsl_matrix*, gsl_matrix*); //X, submatrix
	int shuffle_by_distance(gsl_matrix*);
	int get_elementwise_mean(vector<gsl_matrix*>, gsl_matrix*); //list of matrices, mean matrix
	int get_empirical_pval(gsl_vector*, gsl_vector*, gsl_vector*); //x, sorted_background, pval
	int correct_for_fdr(gsl_vector*, double, gsl_vector*, gsl_vector*); // pval, FDR(alpha), qval, significance	
	int get_cosine_distance_to_upstream_neighbor(gsl_matrix*, gsl_vector*); //U, dist	
	double logsumexp_vector(gsl_vector*);
	double logsumexp_matrix(gsl_matrix*); 
	int mvn_logpdf_per_row(gsl_matrix*, gsl_vector*, gsl_matrix*, gsl_vector*);//X, mu, Covars, pdfs
	int get_all_binary_strings_in_N_bits(gsl_matrix*, int, int, int&, double[]);
	int get_all_binary_strings_in_N_bits(gsl_matrix*, int);
	int print_progress_bar(double);
	int get_csr_submatrix(gsl_spmatrix*, gsl_matrix*, int, int);
	int get_local_pval_from_neighbors(gsl_vector*, gsl_vector*, gsl_vector*, int);
	int pick_summit(gsl_vector*, gsl_vector*, gsl_vector*);
};
#endif
