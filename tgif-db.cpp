#include <iostream>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_spmatrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_statistics_double.h>
#include <list>
#include <vector>
#include <math.h>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "modules/io.h"
#include "modules/nmf.h" 
#include "modules/utils.h"
#include "modules/node.h"
#include "modules/leaf.h"
#include "modules/initialization.h"
#include "modules/root.h"
#include "modules/branch.h"
#include "modules/tgif.h"
#include "modules/annotation.h"

int main(int argc, char **argv)
{
	struct timeval beginTime;
	gettimeofday(&beginTime,NULL);

	struct rusage bUsage;
	getrusage(RUSAGE_SELF,&bUsage);

	const char* treeFile;
	const char* bedFile;
	bool diff_compartment = false; //by default, do differential TAD analysis
	int original_N = -1;
	int N = -1;
	int windowSize = 2000000;
	int stepSize = 1000000;
	int minK = 2;
	int maxK = 8;

	double FDR = 0.05;
	string outputPrefix = string("");
	int randomState = 1010;
	bool verbose = true;
	bool log = false;
	int maxIter = 300; // used in TGIF, will keep fixed
	double tol = 1; // used in TGIF, will keep fixed
	double alpha = 1000000;
	double threshold = 0.5; // sparsity threshold for filtering rows/cols
	bool outputBoundaryScore = false;
	bool outputPval = false;
	bool outputQval = false;
	string usage = string("usage.txt");

	int c;
	while((c = getopt(argc, argv, "w:s:m:x:o:r:a:f:t:bpqlh")) != -1)
		switch (c) {
			case 'w':
				windowSize = atoi(optarg);
				break;
			case 's':
				stepSize = atoi(optarg);
				break;
			case 'm':
				minK = atoi(optarg);
				break;
			case 'x':
				maxK = atoi(optarg);
				break;
			case 'o':
				outputPrefix = string(optarg);
				break;
			case 'r':
				randomState = atoi(optarg);
				break;
			case 'a':
				alpha = atof(optarg);
				break;
			case 't':
				threshold = atof(optarg);
				break;
			case 'b':
				outputBoundaryScore = true;
				break;
			case 'f':
				FDR = atof(optarg);
				break;
			case 'p':
				outputPval = true;
				break;
			case 'q':
				outputQval = true;
				break;
			case 'l':
				log = true;
				break;
			case 'h':
				io::print_usage(usage);
				return 0;
			case '?':
				io::print_usage(usage);
				return 0;
			default:
				io::print_usage(usage);
				return 0;
		}	

	if ((argc - optind) != 2) {
		io::print_usage(usage);
		return 1;
	} else {
		treeFile = argv[optind];
		bedFile = argv[optind+1];
	}

	string treeFileName(treeFile);
	string bedFileName(bedFile);
	string chro;
	vector<int> parentIds, coords;
	vector<string> aliases, fileNames;

	cout << "[TGIF in differential boundary analysis mode]" << endl;
	cout << "Reading in input files..." << endl;	
	int binSize;
	io::read_metadata(bedFileName, binSize, original_N, coords, chro);
	io::read_tree(treeFileName, parentIds, aliases, fileNames);
	int nTasks = fileNames.size();
	int nNodes = parentIds.size();

	// sliding window size for generating submatrices
	windowSize /= binSize;
	stepSize /= binSize;
	if (windowSize < 100) {
		cout << "The sliding window size was reset to required minimum of 100 bins/regions and the stride to 50." << endl;
		windowSize = 100;
		stepSize = 50;
	}

	// remove sprase rows
	int map[original_N];
	io::map_nonzero_rows(original_N, map, fileNames, N, windowSize, threshold);
	cout <<"(Filtered out " << original_N-N << "/" << original_N << " regions with fewer than half of neighbors with non-zero interactions within the specified window size in at least one of the input datasets.)" << endl;

	// read full matrices
	vector<gsl_spmatrix*> fullMatrices;
	int nzmax = N + (windowSize-1) * (2*N - windowSize); 
	for (int i = 0; i < nTasks; i++) {
		gsl_spmatrix* tmp = gsl_spmatrix_alloc_nzmax(N, N, nzmax, GSL_SPMATRIX_COO);
		gsl_spmatrix* M = gsl_spmatrix_alloc_nzmax(N, N, nzmax, GSL_SPMATRIX_CSR);
		io::read_into_coo_matrix_within_distance(fileNames[i], windowSize, tmp, map);
		gsl_spmatrix_csr(M, tmp);
		gsl_spmatrix_free(tmp);
		fullMatrices.push_back(M);		
	}

	// storage for boundary score calcuation
	gsl_matrix* boundaryScore = gsl_matrix_calloc(N, nTasks);
	gsl_vector* numFactorizations = gsl_vector_calloc(N);
	gsl_vector* backgroundScore = gsl_vector_calloc(N);

	// placeholder submatrices to copy into
	gsl_matrix* bg_subX = gsl_matrix_alloc(windowSize, windowSize);
	vector<gsl_matrix*> submatrices;
	for (int t = 0; t < nTasks; t++) {
		gsl_matrix* subX = gsl_matrix_alloc(windowSize, windowSize);
		submatrices.push_back(subX);
	}

	// for printing progress bar
	double checkpoint_arr[] = {0.01, 0.05, 0.09, 0.15, 0.31, 0.47, 0.63, 0.80, 0.96};
	vector<double> checkpoints (checkpoint_arr, checkpoint_arr + sizeof(checkpoint_arr)/sizeof(double));
	vector<double>::iterator progress = checkpoints.begin();
	int cntFactorizations = 0;
	double totalFactorizations = N/stepSize * (maxK - minK + 1);

	bool lastMatDone = false; // to break out of loop below if we reach the end
	gsl_vector* tmp = gsl_vector_calloc(windowSize); // used for calculating boundary score

	cout << "Performing TGIF..." << endl;
	for (int i = 0; i < N; i += stepSize) { 
		// get start and end index of this window
		int upper_corner = i;
		if ((upper_corner + windowSize) > N) {
			if (lastMatDone) {
				break;
			} else{
				upper_corner = N - windowSize;
				lastMatDone = true;
			}
		}

		// get submatrices
		for (int t = 0; t < nTasks; t++) {
			utils::get_csr_submatrix(fullMatrices[t], submatrices[t], upper_corner, upper_corner+windowSize);	
		}

		//making background signal matrix
		utils::get_elementwise_mean(submatrices, bg_subX);
		utils::shuffle_by_distance(bg_subX);

		for (int k = minK; k < (maxK+1); k++) {
			// do TGIF in this window
			TGIF tgif = TGIF(alpha, k, maxIter, tol, randomState, false, false, outputPrefix);
			tgif.make_tree(parentIds, aliases, nTasks);
			tgif.fit(submatrices, cntFactorizations % nTasks);

			// get boundary scores in this window
			vector<Node*>& tree = tgif.get_tree(); 
			for (int t = 0; t < nTasks; t++) {
				gsl_matrix* U = ((Leaf*) tree[t])->get_U();
				utils::get_cosine_distance_to_upstream_neighbor(U, tmp);
				gsl_vector_view b = gsl_matrix_subcolumn(boundaryScore, t, upper_corner, windowSize);
				gsl_vector_add(&b.vector, tmp);
			}
			tgif.clean();

			// Do NMF on background matrix
			gsl_matrix* bg_U = gsl_matrix_calloc(k, windowSize);
			gsl_matrix* bg_V = gsl_matrix_calloc(k, windowSize);
			NMF nmf = NMF(k, nndsvd_init, 200, randomState, false, 1);
			nmf.fit(bg_subX, bg_U, bg_V);
			utils::get_cosine_distance_to_upstream_neighbor(bg_U, tmp);
			gsl_matrix_free(bg_U);
			gsl_matrix_free(bg_V);

			// get background boundary scores in this window
			gsl_vector_view b = gsl_vector_subvector(backgroundScore, upper_corner, windowSize);
			gsl_vector_add(&b.vector, tmp);

			// increment number of factorizations done for regions in this window
			gsl_vector_view stretch = gsl_vector_subvector(numFactorizations, upper_corner, windowSize);
			gsl_vector_add_constant(&stretch.vector, 1);

			// occasionally print progress bar	
			cntFactorizations++;
			if ((cntFactorizations/totalFactorizations) >= *progress && progress != checkpoints.end()) {
				utils::print_progress_bar(*progress);
				++progress;
			} 
		}
	}		

	// free storage no longer needed
	gsl_vector_free(tmp); 
	gsl_matrix_free(bg_subX);
	for (int t = 0; t < submatrices.size(); t++) {
		gsl_matrix_free(submatrices[t]);
	} 	
	for (int i = 0; i < fullMatrices.size(); i++) {
		gsl_spmatrix_free(fullMatrices[i]);
	}

	// get mean boundary score for each region
	gsl_vector_div(backgroundScore, numFactorizations);
	for (int i = 0; i < N; i++) {
		int denom = gsl_vector_get(numFactorizations, i);
		if (denom == 0) {
			gsl_vector_set(numFactorizations, i, 1.0);
		} else {
			gsl_vector_set(numFactorizations, i, 1.0/denom);
		}
	}
	gsl_matrix_scale_rows(boundaryScore, numFactorizations);
	gsl_vector_free(numFactorizations);

	// optionally print boundary scores:
	string header = aliases[0];
	for (int i =1; i < nTasks; i++) {
		header += "\t" + aliases[i];
	}
	if (outputBoundaryScore) {
		io::write_dense_matrix_with_bed_header_and_map(outputPrefix + "boundary_score.txt", chro, coords, binSize, header, boundaryScore, map, original_N);
		io::write_vector(outputPrefix + "background_score.txt",backgroundScore);
	}

	// get empirical p-val and do FDR correction
	cout << "Testing for significant boundaries..." << endl;
	gsl_matrix* Pvals = gsl_matrix_alloc(N, nTasks);
	gsl_matrix* Qvals = gsl_matrix_alloc(N, nTasks);
	gsl_matrix* rejectNull = gsl_matrix_alloc(N, nTasks);
	gsl_matrix* summit_only = gsl_matrix_calloc(N,nTasks);
	gsl_sort_vector(backgroundScore);
	for (int i = 0; i < nTasks; i++ ) {
		gsl_vector_view p = gsl_matrix_column(Pvals, i);
		gsl_vector_view q = gsl_matrix_column(Qvals, i);
		gsl_vector_view r = gsl_matrix_column(rejectNull, i);
		gsl_vector_view s = gsl_matrix_column(boundaryScore, i);
		gsl_vector_view so = gsl_matrix_column(summit_only, i);
		utils::get_empirical_pval(&s.vector, backgroundScore, &p.vector);	
		//utils::get_local_pval_from_neighbors(&s.vector, backgroundScore, &p.vector, windowSize);
		utils::correct_for_fdr(&p.vector, FDR, &q.vector, &r.vector);
		utils::pick_summit(&r.vector, &s.vector, &so.vector);
	}
	gsl_vector_free(backgroundScore);
	io::write_dense_matrix_with_bed_header_and_map(outputPrefix + "significant_boundaries.txt", chro, coords, binSize, header, rejectNull, map, original_N);
	io::write_dense_matrix_with_bed_header_and_map(outputPrefix + "significant_boundaries_summit_only.txt", chro, coords, binSize, header, summit_only, map, original_N);

	if (outputPval) {
		io::write_dense_matrix_with_bed_header_and_map(outputPrefix + "boundary_pval.txt",chro, coords, binSize, header, Pvals, map, original_N);
		io::write_dense_matrix_with_bed_header_and_map(outputPrefix + "boundary_adj_pval.txt",chro, coords, binSize, header, Qvals, map, original_N);
	}
	gsl_matrix_free(Pvals);
	gsl_matrix_free(Qvals);

	// Annotate each region for context specificity
	//Annotations annot = Annotations(nTasks, aliases);
	//annot.assign_states(rejectNull);
	//annot.assign_states(summit_only);
	//annot.print_annotations(outputPrefix + "annotations.txt", chro, coords, binSize, original_N, map);

	/***
 	*Pairwise significantly differential boundaries
 	***/
	gsl_vector* dist = gsl_vector_alloc(N);
	int loss[N];
	gsl_vector* pval = gsl_vector_alloc(N);
	gsl_vector* qval = gsl_vector_alloc(N);
	gsl_vector* reject_null = gsl_vector_alloc(N);

	for (int i=0; i < nTasks; i++) {
		for (int j=i+1; j < nTasks; j++) {
			// get abs difference in boundary score between timepoint i and j for region idx
			for (int idx =0; idx < N; idx++) {
				double diff = gsl_matrix_get(boundaryScore, idx, i) - gsl_matrix_get(boundaryScore, idx, j);
				if (diff < 0) {
					diff  = diff * -1;
					loss[idx] = i;	
				} else {
					loss[idx] = j;
				}
				gsl_vector_set(dist, idx, diff);
			}		
			// get mean and std of differences in consistently non-boundary regions
			double distances_in_consistent_regions[N];
			int count_consistent_regions = 0;
			for (int idx=0;idx < N; idx++) {
				if ((gsl_matrix_get(summit_only, idx, i) == 0) && (gsl_matrix_get(summit_only, idx, j)==0)) { 
					distances_in_consistent_regions[count_consistent_regions] = gsl_vector_get(dist,idx);
					count_consistent_regions++;
				}
			}
			double mean = gsl_stats_mean(distances_in_consistent_regions, 1, count_consistent_regions);
			double std = gsl_stats_sd_m(distances_in_consistent_regions, 1, count_consistent_regions, mean);

			// get p-value (gaussian CDF upper tail) based on given mean and std
			for (int idx=0; idx < N; idx++) {
				double x = gsl_vector_get(dist, idx) - mean;
				gsl_vector_set(pval, idx, gsl_cdf_gaussian_Q(x, std));
			}

			// FDR correction
			utils::correct_for_fdr(pval, 0.05, qval, reject_null);
			// output significantly differential boundary regions
			string pairFileName = outputPrefix + aliases[i] + "_vs_" + aliases[j] + "_significantly_differential_boundary_regions"; 
			bool writeToFile[N];
			for (int idx = 0; idx < N; idx++) {
				bool changed = false;
				if ((gsl_matrix_get(summit_only, idx, i) - gsl_matrix_get(summit_only, idx, j)) != 0) {
					changed = true;
				}
				writeToFile[idx] = changed;
			} 
			//io::write_significant_regions(pairFileName + ".txt", chro, coords, binSize, writeToFile, dist, pval, qval, reject_null, map, original_N);	
			io::write_sigDB(pairFileName + ".txt", chro, coords, binSize, writeToFile, dist, pval, qval, reject_null, map, original_N, aliases, loss);	
			//io::write_significant_regions_for_debug(pairFileName + "_DEBUG.txt", chro, coords, binSize, writeToFile, dist, pval, qval, reject_null, map, original_N);
		}
	}

	gsl_matrix_free(summit_only);
	gsl_matrix_free(rejectNull);
	gsl_matrix_free(boundaryScore);
	gsl_vector_free(pval);
	gsl_vector_free(qval);
	gsl_vector_free(reject_null);
	gsl_vector_free(dist);

	// Output parameters info log optionallly
	if (log) {
		ofstream ofs;
		ofs.open((outputPrefix + "parameters.log").c_str());
		ofs << "alpha = " << alpha << endl;
		ofs << "FDR = " << FDR << endl;
		ofs << "window size = " << windowSize * binSize << endl;
		ofs << "step size = " << stepSize *binSize << endl;
		ofs << "min k = " << minK << endl;
		ofs << "max k = " << maxK << endl;
		ofs.close();
	}


	struct timeval endTime;
	gettimeofday(&endTime,NULL);

	struct rusage eUsage;
	getrusage(RUSAGE_SELF,&eUsage);

	unsigned long int bt = beginTime.tv_sec;
	unsigned long int et = endTime.tv_sec;

	cout << "Total time elapsed: " << et - bt << " seconds" << endl;

	unsigned long int bu = bUsage.ru_maxrss;
	unsigned long int eu = eUsage.ru_maxrss;
	
	cout << "Memory usage: " << (eu - bu)/1000 << "MB" << endl;
	
	return 0;
}
