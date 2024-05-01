#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_permute.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_statistics_double.h>
#include <list>
#include <vector>
#include <map>
#include <math.h>
#include <iostream>
#include "utils.h"
#include "io.h"

int utils::get_laplacian(gsl_matrix* W, gsl_matrix* L) {
	int n = W->size1;
	gsl_vector* D = gsl_vector_alloc(n);
	for (int i = 0; i < n; i++) {
		double *val = &(D->data[i]);	
		gsl_vector_view row = gsl_matrix_row(W, i);
		*val = gsl_blas_dasum(&row.vector);
	}
	gsl_vector_view diag = gsl_matrix_diagonal(L);
	gsl_vector_memcpy(&diag.vector, D);
	gsl_matrix_sub(L, W);
	gsl_vector_free(D);
	return 0;
}

int utils::get_neighborhood_graph(gsl_matrix* W, int neighborhood_radius) {
	gsl_vector_view diag = gsl_matrix_diagonal(W);
	gsl_vector_set_all(&diag.vector, 1);
	for (int i = 1; i <= neighborhood_radius; i++) {
		gsl_vector_view lower = gsl_matrix_subdiagonal(W, i);
		gsl_vector_set_all(&lower.vector, 1);
		gsl_vector_view upper = gsl_matrix_superdiagonal(W, i);
		gsl_vector_set_all(&upper.vector, 1);
	}
	return 0;
}

double utils::get_frobenius_norm(gsl_matrix* X) {
	int n = X->size1;
	int m = X->size2;
	double sum = 0;
	gsl_vector* row_copy = gsl_vector_alloc(m);
	for (int i = 0; i < n; i++) {
		gsl_vector_view row = gsl_matrix_row(X, i);
		gsl_vector_memcpy(row_copy, &row.vector); 
		gsl_vector_mul(row_copy, &row.vector);
		double rowsum = gsl_blas_dasum(row_copy);
		sum += rowsum;
	}
	gsl_vector_free(row_copy);
	return sum;
}

int utils::get_inverse(gsl_matrix* A) {
	gsl_linalg_cholesky_decomp1(A);
	gsl_linalg_cholesky_invert(A);
	return 0;
}

int utils::get_elementwise_mean(vector<gsl_matrix*> matlist, gsl_matrix* M) {
	int n = matlist.size();
	gsl_matrix_set_zero(M);
	for(vector<gsl_matrix*>::iterator mat=matlist.begin(); mat != matlist.end(); ++mat) {
		gsl_matrix_add(M, *mat);
	}
	gsl_matrix_scale(M, 1.0/n);
	return 0;
}

int utils::sample_by_distance(gsl_matrix* X, gsl_matrix* S) {
	gsl_rng_env_setup();
	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	int n = X->size1;
	int m = S->size1; //size of submatrix to be filled with samples from X
	for (int k = 0; k < m; k++) {
		//shuffle indices of entries in the given diagonal of X
		gsl_permutation* p = gsl_permutation_alloc(n-k);
		gsl_permutation_init(p);
		gsl_ran_shuffle(r, p->data, n-k, sizeof(size_t));
		//fill the digonal of S with random samples from X
		for (int i = 0; i < (m-k); i++) {
			int sample_i = gsl_permutation_get(p, i);
			double sample_val = gsl_matrix_get(X, sample_i, sample_i+k);
			int j = i + k;;
			gsl_matrix_set(S, i, j, sample_val);
			gsl_matrix_set(S, j, i, sample_val);
		}
		gsl_permutation_free(p);
	}
	gsl_rng_free(r);
	return 0;
}

int utils::shuffle_by_distance(gsl_matrix* X) {
	gsl_rng_env_setup();
	const gsl_rng_type* T = gsl_rng_default;
	gsl_rng* r = gsl_rng_alloc(T);

	int n = X->size1;
	//shuffle entries in the given diagonal of X
	for (int k = 0; k < n-1; k++) {
		gsl_permutation* p = gsl_permutation_alloc(n-k);
		gsl_permutation_init(p);
		gsl_ran_shuffle(r, p->data, n-k, sizeof(size_t));
		if (k == 0) {
			gsl_vector_view diag = gsl_matrix_diagonal(X);	
			gsl_permute_vector(p, &diag.vector);
		} else {
			gsl_vector_view subdiag = gsl_matrix_subdiagonal(X, k);
			gsl_vector_view supdiag = gsl_matrix_superdiagonal(X, k);
			gsl_permute_vector(p, &subdiag.vector);
			gsl_permute_vector(p, &supdiag.vector);
		}
		gsl_permutation_free(p);
	}
	gsl_rng_free(r);
	return 0;
}


int utils::get_empirical_pval(gsl_vector* x, gsl_vector* sorted_background, gsl_vector* pval) {
	int n = x->size;
	for (int i = 0; i < n; i++) {
		double val = gsl_vector_get(x, i);
		for (int j = 0; j < n; j++) {
			if (gsl_vector_get(sorted_background,j) >= val) {
				gsl_vector_set(pval, i, (n-j)*1.0/n);
				break;
			}
			if ((j+1)==n) {
				gsl_vector_set(pval, i, 0);
			}
		} 
	}
	return 0;
}

int utils::correct_for_fdr(gsl_vector* pval, double fdr, gsl_vector* qval, gsl_vector* reject_null) {
	size_t n = pval->size;
	//get rank of each p-value after sorting from smallest to largest
	//(rank is the ivnerse of the index permutation)
	gsl_permutation* perm = gsl_permutation_alloc(n);
	gsl_permutation* rank = gsl_permutation_alloc(n);
	gsl_sort_vector_index(perm, pval);
	gsl_permutation_inverse(rank, perm);
	double mpval = 0; // largest p-value that is smaller than its q-value
	for (int i =0; i<n; i++) {
		double q = (rank->data[i]+1)*fdr/n; //q-value = (rank/n) * FDR
		gsl_vector_set(qval, i, q);
		double p = gsl_vector_get(pval, i);
		//cout << p << "\t" << rank->data[i] << "\t" <<  q << endl;
		if ((p < q) && (p > mpval)) {
			mpval = p;
		}
	}
	//cout << mpval << endl;
	for (int i =0; i<n; i++) {
		if (gsl_vector_get(pval, i) <= mpval) {
			gsl_vector_set(reject_null, i, 1);
		} else {
			gsl_vector_set(reject_null, i, 0);
		}
	}
	gsl_permutation_free(perm);
	gsl_permutation_free(rank);
	return 0;
}

int utils::get_cosine_distance_to_upstream_neighbor(gsl_matrix* U, gsl_vector* dist) {
	int n = dist->size;
	gsl_vector_set(dist, 0, 0);
	gsl_vector_view last_col = gsl_matrix_column(U, 0);
	double upstream_neighbor_norm = gsl_blas_dnrm2(&last_col.vector);
	for (int i = 1; i<n; i++) {
		gsl_vector_view curr_col = gsl_matrix_column(U, i);
		double curr_col_norm = gsl_blas_dnrm2(&curr_col.vector);
		double denom;
		gsl_blas_ddot(&last_col.vector, &curr_col.vector, &denom);
		double num = upstream_neighbor_norm * curr_col_norm;
		if (num == 0) {
			gsl_vector_set(dist, i, 0);
		} else {
			gsl_vector_set(dist, i, 1 - (denom/num));
		}
		last_col = curr_col;
		upstream_neighbor_norm = curr_col_norm;
	} 
	return 0;
}

double utils::logsumexp_vector(gsl_vector* x) {
	size_t n = x->size;
	double max_exp = gsl_vector_max(x), sum = 0.0;
 	for (size_t i = 0; i < n; i++) {
		sum += exp(x->data[i * x->stride] - max_exp);
	}
	return log(sum)+max_exp;
}

double utils::logsumexp_matrix(gsl_matrix* X) {
	size_t n = X->size1, m = X->size2, mtda = X->tda;
	double max_exp = gsl_matrix_max(X), sum = 0.0;
 	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < m; j++) {
			sum += exp(X->data[i * mtda + j] - max_exp);
		}
	}
	return log(sum)+max_exp;
}

int utils::mvn_logpdf_per_row(gsl_matrix* X, gsl_vector* mu, gsl_matrix* Covars, gsl_vector* pdf){
	size_t n = X->size1;
	size_t m = X->size2;
	size_t stride = pdf->stride;
	gsl_vector* work = gsl_vector_alloc(m);
	gsl_matrix* L = gsl_matrix_alloc(m,m);
	gsl_matrix_memcpy(L, Covars);
	if (gsl_matrix_isnull(L)) {
		gsl_vector_view diag = gsl_matrix_diagonal(L);
		gsl_vector_set_all(&diag.vector, 0.01);	
	}
	for (int i = 0; i < m; i++) {
		for (int j = i; j < m; j++) {
			double val = gsl_matrix_get(L, i, j);
			if (i == j) {
				if (val <= 0) {
					gsl_matrix_set(L, i, j, 0.01);
				}
			} else {
				if (val < 0 ) {
					gsl_matrix_set(L, i, j, 0);
					gsl_matrix_set(L, j, i, 0);
				}
			}
		}
	}
	io::write_dense_matrix("covars.txt",L);
	gsl_linalg_cholesky_decomp1(L);
	double logpdf = 0.0;
	for (int i = 0; i < n; i++) {
		gsl_vector_view x = gsl_matrix_row(X, i);
		gsl_ran_multivariate_gaussian_log_pdf(&x.vector, mu, L, &(pdf->data[i*stride]), work);	
	}
	gsl_matrix_free(L);
	gsl_vector_free(work);
	return 0;
}

int utils::get_all_binary_strings_in_N_bits(gsl_matrix* S, int N, int colId, int& rowId, double arr[]) {
	if (colId == N) {
		gsl_vector_view a = gsl_vector_view_array(arr, N);
		gsl_matrix_set_row(S, rowId, &a.vector);
		rowId++;		
		return 0;
	}
	arr[colId] = 0;
	get_all_binary_strings_in_N_bits(S, N, colId +1, rowId, arr);

	arr[colId] = 1;
	get_all_binary_strings_in_N_bits(S, N, colId +1, rowId, arr);

	return 0;
}

int utils::get_all_binary_strings_in_N_bits(gsl_matrix* S, int N) {
	int rowId = 0;
	double arr[N];
	get_all_binary_strings_in_N_bits(S, N, 0, rowId, arr);
	return 0;
}

int utils::print_progress_bar(double progress) {
	int barWidth = 70;
	cout << "[";
	int pos = barWidth * progress;
	for (int i = 0; i < barWidth; ++i) {
		if (i < pos) cout << "=";
		else if (i == pos) cout << ">";
		else cout << " ";
	}
	cout << "] " << int(progress * 100.0) << "%" << endl;
	return 0;
}

int utils::get_csr_submatrix(gsl_spmatrix* X, gsl_matrix* block, int start, int end) {
	//assumes X is a CSR matrix
	//block is a symmetric matrix smaller than X
	//start and end are the min/max row(col) indices of the submatrix to extract
	gsl_matrix_set_zero(block);
	for (int row = start; row < end; row++) {
		int ptr1 = X->p[row];
		int ptr2 = X->p[row+1];
		for (int idx = ptr1; idx < ptr2; idx++) {
			int col = X->i[idx];
			if ((col >= start) && (col < end)) {
				gsl_matrix_set(block, row - start, col - start, X->data[idx]);
			}
		}
	}
	return 0;
}

int utils::get_local_pval_from_neighbors(gsl_vector* x, gsl_vector* background, gsl_vector* pval, int num_neighbors) {
	// assumes x, background, pval are all the same length
	int n = x->size;
	for (int i = 0; i < n; i++) {
		double val = gsl_vector_get(x,i);	
		int start = i - num_neighbors;
		int end = i + num_neighbors + 1;
		if (start < 0) {
			end = min(end - start, n);
			start = 0;
		}
		if (end > n) {
			start = max(start - (end-n), 0);
			end = n;
		}
		//cout << start << "\t" << end << endl;
		int more = 0;
		int denom = 0;
		for (int j = start; j < end; j++) {
			denom++;
			if (val <= gsl_vector_get(background,j)) {
				more = more+1;
			}
		}
		gsl_vector_set(pval, i, more * 1.0 / denom);
	}
	return 0;
}

int utils::pick_summit(gsl_vector* rejectNull, gsl_vector* score, gsl_vector* summit_only) {
	int n = rejectNull->size;
	double curr_summit_score = gsl_vector_get(score, 0);
	int curr_summit_idx = 0;
	double last_val = gsl_vector_get(rejectNull, 0);
	for (int i = 1; i < n; i++) {
		double curr_val = gsl_vector_get(rejectNull, i);
		double curr_score = gsl_vector_get(score, i);
		if (curr_val == 1) { 
			if (last_val == 1) {
				if (curr_score > curr_summit_score) {
					// if in the middle of a stretch of significant boundaries and we find a new summit,
					// update summit tracker
					curr_summit_score = curr_score;
					curr_summit_idx = i;
				}
			} else  { //if a new stretch of significant boundaries started, reset summit tracker
				curr_summit_score = curr_score;
				curr_summit_idx = i;
			}
		} else { //if a stretch of significant boundaries ended, set the summit/peak
				if (last_val == 1) {
					gsl_vector_set(summit_only, curr_summit_idx, 1);
				}
		}
		last_val = curr_val; 
	}
	gsl_vector_set(summit_only, curr_summit_idx, 1);
	return 0;
}

int utils::get_rowwise_correlation(gsl_matrix* X, gsl_matrix* C) {
	int n = X->size1;
	int m = X->size2;
	for (int i = 0; i < n; i++) {
		gsl_vector_view x_i = gsl_matrix_row(X, i);
		bool i_all_zero = gsl_vector_isnull(&x_i.vector);
		if (i_all_zero == true) {
			gsl_vector_view c_i_ = gsl_matrix_row(C, i);
			gsl_vector_set_zero(&c_i_.vector);
			gsl_vector_view c__i = gsl_matrix_column(C, i);
			gsl_vector_set_zero(&c__i.vector);
		} else {
			gsl_matrix_set(C, i, i, 1);
			for (int j = i+1; j < n; j++) {
				gsl_vector_view x_j = gsl_matrix_row(X, j);
				double c = 0;
				bool j_all_zero = gsl_vector_isnull(&x_j.vector);
				if (j_all_zero == false) {
					c = gsl_stats_correlation(x_i.vector.data,1,x_j.vector.data,1,n);
				}
				gsl_matrix_set(C, i, j, c);
				gsl_matrix_set(C, j, i, c);
			}
		} 
	}
	return 0;
}

int utils::get_columnwise_max_idx(gsl_matrix* X, gsl_vector* idx) {
	int m = X->size2;
	for (int i = 0; i < m; i++) {
		gsl_vector_view x_i = gsl_matrix_column(X, i);
		gsl_vector_set(idx, i, gsl_vector_max_index(&x_i.vector));
	}
	return 0;
}

int utils::distancewise_normalization(gsl_matrix* X) {
	int n = X->size1;
	for (int i = 0; i < n; i++) {
		gsl_vector_view diag = gsl_matrix_superdiagonal(X, i);
		//double mean = gsl_stats_mean(diag.vector.data, 1, n-i);
		double mean = 0;
		for (int j = 0; j < (n-i); j++) {
			mean += gsl_vector_get(&diag.vector, j);
		}
		mean /= (n-i);
		if (mean > 0) {
			gsl_vector_scale(&diag.vector, 1.0/mean);
			if (i > 0) {
				gsl_vector_view subdiag = gsl_matrix_subdiagonal(X, i);
				gsl_vector_memcpy(&subdiag.vector, &diag.vector);
			}
		}
	}
	return 0;
}

int utils::get_rowwise_mean(gsl_matrix* X, gsl_vector* mean) {
	size_t n = X->size1, m=X->size2, mtda=X->tda;
	for (size_t i = 0; i < n; i++) {
		double sum = 0;
		for (size_t j = 0; j < m; j++) {
			sum += X->data[i * mtda + j];
		}
		gsl_vector_set(mean, i, sum/((double) m));
	}
	return 0;
}


map<int,char> utils::map_cluster_to_compartment(gsl_matrix* cluster_id, gsl_matrix* mean_oe) {
	//cluster whose regions have higher mean O/E value (stored in M) 
	//is assigned to compartment A
	map<int,double> cluster_mean;
	map<int,int> cluster_num_regions;
	map<int,char> cid_to_ab;
	int n = cluster_id->size1;
	int T = cluster_id->size2;
	for (int i = 0; i < n; i++) {
		for (int t = 0; t < T; t++) {
			int cid = (int) gsl_matrix_get(cluster_id, i, t);
			cluster_mean[cid] += gsl_matrix_get(mean_oe, i, t);
			cluster_num_regions[cid] += 1;
		}
	}
	for (int i = 0; i < 2; i++) {
		cluster_mean[i] /= cluster_num_regions[i];
	}
	if (cluster_mean[0] > cluster_mean[1]) {
		cid_to_ab[0] = 'A';
		cid_to_ab[1] = 'B';
	} else {
		cid_to_ab[0] = 'B';
		cid_to_ab[1] = 'A';
	}
	return cid_to_ab;
}

int utils::get_norm_by_row(gsl_matrix* X, gsl_vector* norm) {
	int n = X->size2;
	for (int i=0; i < n; i++) {
		gsl_vector_view x = gsl_matrix_column(X, i);
		gsl_vector_set(norm, i, gsl_blas_dnrm2(&x.vector));
	}
	return 0;
}

int utils::get_cosine_distance_by_row(gsl_matrix* X, gsl_vector* x_norm, gsl_matrix* Y, gsl_vector* y_norm, gsl_vector* dist) {
	int n = X->size2;
	for (int i=0; i < n; i++) {
		gsl_vector_view x = gsl_matrix_column(X, i);
		gsl_vector_view y = gsl_matrix_column(Y, i);
		double num = 0;
		gsl_blas_ddot(&x.vector, &y.vector, &num);
		double denom = gsl_vector_get(x_norm,i) * gsl_vector_get(y_norm,i);
		gsl_vector_set(dist, i, 1-num/denom); 
	}
	return 0;
}
