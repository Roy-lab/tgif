#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_matrix.h>
#include <list>
#include <vector>
#include <map>
#include <string>
#ifndef _io_
#define _io_
using namespace std;
namespace io
{
	int read_sparse_matrix(string, gsl_matrix*);
	int map_nonzero_rows(int, int[], vector<string>, int&, int, double);
	int read_sparse_matrix_nonzero_rows_only(string, gsl_matrix*, int[]);
	int read_dense_matrix(string, gsl_matrix*, bool, string&);
	int write_dense_matrix(string, gsl_matrix*);
	int write_dense_matrix_with_header(string, string, gsl_matrix*);
	int write_dense_matrix_with_header_and_map(string, string, gsl_matrix*, int[], int);
	int write_cluster_label(string, gsl_vector*);
	int write_vector(string, gsl_vector*);
	int read_vector(string, gsl_vector*);
	int write_list(string, list<double>&);
	int read_tree(string, vector<int>&, vector<string>&, vector<string>&);
	int read_metadata(string, int&, int&, vector<int>&, string&);
	int read_into_coo_matrix_within_distance(string, int, gsl_spmatrix*, int[]);
	int print_usage(string);
	int write_compartments_to_file(string,string,vector<int>&,int,string, gsl_matrix*, map<int,char>);
};
#endif
