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
	int original_N = -1;
	int N = -1;
	int k = 2;
	string outputPrefix = string("");
	int randomState = 1010;
	bool verbose = true;
	bool log = false;
	int maxIter = 300; // used in TGIF, will keep fixed
	double tol = 1; // used in TGIF, will keep fixed
	double alpha = 10000;
	bool outputUandV = false; 
	bool inputAlreadyOE = false; // if input matrices are already O/E counts
	string usage = string("usage.txt");

	int c;
	while((c = getopt(argc, argv, "k:o:r:a:eulh")) != -1)
		switch (c) {
			case 'k':
				k = atoi(optarg);
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
			case 'e':
				inputAlreadyOE = true;
				break;
			case 'u':
				outputUandV = true;
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

	cout << "[TGIF in differential compartment analysis mode]" << endl;
	cout << "Reading in input files..." << endl;	
	int binSize;
	io::read_metadata(bedFileName, binSize, original_N, coords, chro);
	io::read_tree(treeFileName, parentIds, aliases, fileNames);
	int nTasks = fileNames.size();
	int nNodes = parentIds.size();

	int idx_map[original_N];
	io::map_nonzero_rows(original_N, idx_map, fileNames, N, original_N, 0);

	// get correlation matrices 
	vector<gsl_matrix*> corrMatrices;
	gsl_matrix* OE = gsl_matrix_calloc(N, N);	
	gsl_matrix* mean_OE_per_region = gsl_matrix_alloc(N, nTasks);
	for (int i = 0; i < nTasks; i++) {
		//io::read_sparse_matrix(fileNames[i], OE);
		io::read_sparse_matrix_nonzero_rows_only(fileNames[i], OE, idx_map);
		// 1. get O/E matrix, if given raw counts:
		if (!inputAlreadyOE) {
			utils::distancewise_normalization(OE);
		}
		// used later to determine A/B compartment from cluster IDs:
		gsl_vector_view mean_vector = gsl_matrix_column(mean_OE_per_region, i);
		utils::get_rowwise_mean(OE, &mean_vector.vector);
		// 2. get correlation of O/E matrix and upshift by 1:
		gsl_matrix* C = gsl_matrix_alloc(N, N);
		utils::get_rowwise_correlation(OE, C);
		gsl_matrix_add_constant(C, 1);
		corrMatrices.push_back(C);	
		//io::write_dense_matrix(outputPrefix+aliases[i]+"_corr.txt",C);
		gsl_matrix_set_zero(OE);	
	}
	gsl_matrix_free(OE);

	// run tgif
	TGIF tgif = TGIF(alpha, k, maxIter, tol, randomState, false, false, outputPrefix);
	tgif.make_tree(parentIds, aliases, nTasks);
	tgif.fit(corrMatrices, 0);

	// print factors & leaf cluster assignments
	if (outputUandV) {
		tgif.print_factors();
	}
	gsl_matrix* maxDimClusters = gsl_matrix_alloc(N, nTasks);
	vector<Node*>& tree = tgif.get_tree(); 
	for (int i = 0; i < nTasks; i++) {
		gsl_vector_view c = gsl_matrix_column(maxDimClusters, i);
		utils::get_columnwise_max_idx(((Leaf*) tree[i])->get_U(),&c.vector);
	}
	string header = aliases[0];
	for (int i =1; i < nTasks; i++) {
		header += "\t" + aliases[i];
	}
	io::write_dense_matrix_with_bed_header_and_map(outputPrefix+"cluster_assignment.txt",chro, coords, binSize,header,maxDimClusters,idx_map,original_N);

	// if the mean O/E counts in cluster 0 is higher than that of cluster 1, cluster 0 = compartment A
	if (k == 2) {
		//gsl_vector_view cid = gsl_matrix_column(maxDimClusters,0);
		//gsl_vector_view mean_oe0 = gsl_matrix_column(mean_OE_per_region, 0);
		//map<int,char> cid_to_ab = utils::map_cluster_to_compartment(&cid.vector, &mean_oe0.vector);
		map<int,char> cid_to_ab = utils::map_cluster_to_compartment(maxDimClusters, mean_OE_per_region);
		io::write_compartments_to_file(outputPrefix+"compartment_assignment.txt",chro,coords,binSize,header,maxDimClusters,cid_to_ab, idx_map, original_N);
	}

	gsl_matrix_free(mean_OE_per_region);
	// free all correlation matrices
	for (int i = 0; i < corrMatrices.size(); i++) {
		gsl_matrix_free(corrMatrices[i]);
	}

	vector<gsl_matrix*> factors;
	vector<gsl_vector*> norms;

	for (int t=0; t < nTasks; t++) {
		gsl_matrix* factor = ((Leaf*) tree[t])->get_U(); 
		factors.push_back(factor);
		gsl_vector* norm = gsl_vector_alloc(N);
		utils::get_norm_by_row(factor,norm);
		norms.push_back(norm);
	} 

	gsl_vector* dist = gsl_vector_alloc(N);
	gsl_vector* pval = gsl_vector_alloc(N);
	gsl_vector* qval = gsl_vector_alloc(N);
	gsl_vector* reject_null = gsl_vector_alloc(N);

	for (int i=0; i < nTasks; i++) {
		for (int j=i+1; j < nTasks; j++) {
			utils::get_cosine_distance_by_row(factors[i],norms[i],factors[j],norms[j], dist);

			//get mean & std for distances of 'consistent' regions
			double distances_in_consistent_regions[N];
			int count_consistent_regions = 0;
			for (int idx=0;idx < N; idx++) {
				if (gsl_matrix_get(maxDimClusters, idx, i) == gsl_matrix_get(maxDimClusters, idx, j)) { 
					distances_in_consistent_regions[count_consistent_regions] = gsl_vector_get(dist,idx);
					count_consistent_regions++;
				}
			}
			double mean = gsl_stats_mean(distances_in_consistent_regions, 1, count_consistent_regions);
			double std = gsl_stats_sd_m(distances_in_consistent_regions, 1, count_consistent_regions, mean);
			
			//get p-value (gaussian CDF upper tail) for each dist with given mean and std
			for (int idx=0; idx < N; idx++) {
				double x = gsl_vector_get(dist, idx) - mean;
				gsl_vector_set(pval, idx, gsl_cdf_gaussian_Q(x, std));
			}
			//FDR correction
			utils::correct_for_fdr(pval, 0.05, qval, reject_null);	
			//output sig diff regions
			string pairFileName = outputPrefix + aliases[i] + "_vs_" + aliases[j] + "_significantly_differential_compartmental_regions"; 
			bool writeToFile[N];
			for (int idx = 0; idx < N; idx++) {
				bool changed = false;
				if (gsl_matrix_get(maxDimClusters, idx, i) != gsl_matrix_get(maxDimClusters, idx, j)) {
					changed = true;
				}
				writeToFile[idx] = changed;
			} 
			io::write_significant_regions(pairFileName + ".txt", chro, coords, binSize, writeToFile, dist, pval, qval, reject_null, idx_map, original_N);	
			//io::write_significant_regions_for_debug(pairFileName + "_DEBUG.txt", chro, coords, binSize, writeToFile, dist, pval, qval, reject_null, idx_map, original_N);	
		}
	}

	gsl_matrix_free(maxDimClusters);
	gsl_vector_free(pval);
	gsl_vector_free(qval);
	gsl_vector_free(reject_null);
	gsl_vector_free(dist);

	for (int t=0; t < nTasks;t++) {
		gsl_vector_free(norms[t]);
	}

	// Output parameters info log optionallly
	if (log) {
		ofstream ofs;
		ofs.open((outputPrefix + "parameters.log").c_str());
		if (inputAlreadyOE) {
			ofs << "input matrices = O/E counts" << endl;
		} else {
			ofs << "input matrices = raw counts" << endl;
		}
		ofs << "alpha = " << alpha << endl;
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
