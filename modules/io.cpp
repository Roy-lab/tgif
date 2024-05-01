#include <iostream>
#include <fstream>
#include <gsl/gsl_matrix.h>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <map>
#include <math.h>
#include <cmath>
#include <sstream>
#include <vector>
#include "io.h"
#include <gsl/gsl_spmatrix.h>

int io::print_usage(string inputFile) {
	ifstream f(inputFile.c_str());
	string line;
	while(getline(f, line)) {
		cout << line << endl;
	}
	f.close();
	return 0;
}

int io::write_dense_matrix(string outputFile, gsl_matrix* X) {
	int rowNum = X-> size1;
	int colNum = X-> size2;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			ofs << X->data[i * X->tda + j] << "\t";
		}	
		ofs << endl;
	}
	ofs.close();
	return 0;
}

int io::write_dense_matrix_with_header(string outputFile, string header, gsl_matrix* X) {
	int rowNum = X-> size1;
	int colNum = X-> size2;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	ofs << header << endl;
	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			ofs << X->data[i * X->tda + j] << "\t";
		}	
		ofs << endl;
	}
	ofs.close();
	return 0;
}

int io::write_dense_matrix_with_header_and_map(string outputFile, string header, gsl_matrix* X, int map[], int n) {
	int rowNum = n; 
	int colNum = X-> size2;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	ofs << header << endl;
	for (int i = 0; i < rowNum; i++) {
		for (int j = 0; j < colNum; j++) {
			if (map[i] >= 0) {
				ofs << X->data[map[i] * X->tda + j] << "\t";
			} else {
				ofs << "N/A\t";
			}
		}	
		ofs << endl;
	}
	ofs.close();
	return 0;
}


int io::write_cluster_label(string outputFile, gsl_vector* x) {
	int n = x->size;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	for (int i = 0; i < n; i++) {
		ofs << ((int) x->data[i * x->stride]) << endl;
	}
	ofs.close();
	return 0;
}

int io::write_vector(string outputFile, gsl_vector* x) {
	int n = x->size;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	for (int i = 0; i < n; i++) {
		ofs << (x->data[i * x->stride]) << endl;
	}
	ofs.close();
	return 0;
}

int io::read_vector(string inputFile, gsl_vector* x) {
	ifstream input(inputFile.c_str());
	double val;
	int i = 0; 
	while (input >> val) {
		gsl_vector_set(x, i, val);
		i++;
	}
	input.close();
	return 0;
}

int io::map_nonzero_rows(int n, int map[], vector<string> filenames, int &N, int radius, double threshold) {
	int numTasks = filenames.size();
	gsl_matrix* nonzero = gsl_matrix_calloc(n, numTasks);	
	for (int t = 0; t < numTasks; t++) {
		ifstream input(filenames[t].c_str());
		int i, j;
		double cnt;
		while (input >> i >> j >> cnt) {
			if ((abs(j-i)<=radius) && (cnt >0)) {
				(nonzero->data[i * nonzero->tda + t])  += 1;
				(nonzero->data[j * nonzero->tda + t])  += 1;
			}
		}
		input.close();
	}
	int new_idx = 0;
	//double threshold = 0.5;
	for (int i = 0; i < n; i++) {
		bool pass = true;
		int numNeighbors = min(radius, n-i) + min(i, radius) +1;
		for (int t = 0; t < numTasks; t++) {
			// if in at least one of the tasks the row is too sparse skip it
			if (gsl_matrix_get(nonzero, i, t)/numNeighbors <= threshold) {
				pass = false;
				map[i] = -1;
				break;
			} 
		}
		if (pass) {
			map[i] = new_idx;
			new_idx++;
		}
	}
	N = new_idx;
	return 0;
}

int io::read_sparse_matrix_nonzero_rows_only(string inputFile, gsl_matrix* X, int map[]) {
	int rowNum = X->size1;
	int colNum = X->size2;
	ifstream input(inputFile.c_str());
	int i, j;
	double val;
	while (input >> i >> j >> val) {
		//val = log2(val+1);
		int new_i = map[i];
		int new_j = map[j];
		if ((new_i >= 0) && (new_j >= 0)) {
			gsl_matrix_set(X, new_i, new_j, val);
			gsl_matrix_set(X, new_j, new_i, val);
		}
	}
	input.close();
	return 0;
}

int io::read_into_coo_matrix_within_distance(string inputFile, int max_dist, gsl_spmatrix* X, int map[]) {
	int n = X->size1; 
	if ( ((X->sptype) != GSL_SPMATRIX_COO) || (n != (X->size2))) {
		cout << "Allocate a symmetric COO matrix to read the file into." << endl;
		return 1;
	}
	ifstream input(inputFile.c_str());
	int i, j;
	double val;
	while (input >> i >> j >> val) {
		if (abs(j -i) < max_dist) {
			//cout << val << endl;
			int new_i = map[i];
			int new_j = map[j];
			if ((new_i >= 0) && (new_j >= 0)) {
				gsl_spmatrix_set(X, new_i, new_j, val);
				gsl_spmatrix_set(X, new_j, new_i, val);
			}
		}
	}
	input.close();
	return 0;
}

 
int io::read_sparse_matrix(string inputFile, gsl_matrix* X) {
	int rowNum = X->size1;
	int colNum = X->size2;
	ifstream input(inputFile.c_str());
	int i, j;
	double val;
	while (input >> i >> j >> val) {
		//val = log2(val+1);
		gsl_matrix_set(X, i, j, val);
		gsl_matrix_set(X, j, i, val);
	}
	input.close();
	return 0;
}

int io::read_dense_matrix(string inputFile, gsl_matrix* X, bool skipHeader, string& header) {
	int rowNum = X->size1;
	int colNum = X->size2;
	ifstream input(inputFile.c_str());
	char* buff=NULL;
	string buffstr;
	int bufflen=0;
	int i, j;
	double val;
	int rowid=0;
	if (skipHeader) {
		getline(input, header);
		cout << header << endl;
	}
	while(input.good())
	{
		getline(input, buffstr);
		if(buffstr.length()>=bufflen)
		{
			if(buff!=NULL)
			{
				delete[] buff;
			}
			bufflen=buffstr.length()+1;
			buff=new char[bufflen];
		}
		strcpy(buff,buffstr.c_str());
		int colid=0;
		char* tok=strtok(buff,"\t");
		while(tok!=NULL)
		{
			val=atof(tok);
			gsl_matrix_set(X, rowid, colid, val);
			tok=strtok(NULL,"\t");
			colid++;
		}
		rowid++;
	}
	input.close();
	return 0;
}

int io::write_list(string outputFile, list<double>& err) {
	ofstream ofs;
	ofs.open(outputFile.c_str());
	for (list<double>::iterator itr = err.begin(); itr != err.end(); ++itr) {
		ofs << *itr << endl;
	}
	ofs.close();
	return 0;
}

int io::read_tree(string inputFile, 
	vector<int>& parentIds, vector<string>& aliases, vector<string>& fileNames) {
	ifstream input(inputFile.c_str());
	int id, pid;
	string alias, filename;
	while (input >> id >> pid >> alias >> filename) {
		parentIds.push_back(pid);
		aliases.push_back(alias);
		if (filename != "N/A") {
			fileNames.push_back(filename);
		}
	}
	input.close();
	return 0;
}

int io::read_metadata(string inputFile, int &binsize, int &n, vector<int>& coords, string& chro) {
	ifstream f(inputFile.c_str());
	int x,y,idx;
	string s;
	while(f >> s >> x >> y >> idx) {
		if (idx == 0) {
			chro = s;
			binsize = y - x;
		}
		coords.push_back(x);
		n = idx;
	}
	n = n + 1; // zero-indexed
	f.close();
	return 0;
}

int io::write_compartments_to_file(string outputFile,string chro, vector<int>& coord, int binSize,string header, gsl_matrix* X, map<int,char> cid_to_ab, int map[], int n) {
	int rowNum = X-> size1;
	int colNum = X-> size2;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	ofs << "#chro\tstart\tend\t" << header << endl;
	for (int i = 0; i < n; i++) {
		ofs << chro << "\t" << coord[i] << "\t" << coord[i] + binSize;
		int new_i = map[i];
		for (int j = 0; j < colNum; j++) {
			if (new_i < 0) {
				ofs << "\tN/A";
			} else {
				int cid = (int) X->data[new_i*X->tda + j];
				ofs << "\t" << cid_to_ab[cid];
			}
		}
		ofs << endl;
	}
	ofs.close();
	return 0;
}

int io::write_dense_matrix_with_bed_header_and_map(string outputFile, string chro, vector<int>& coord, int binSize, string header, gsl_matrix* X, int map[], int n) {
	int rowNum = n; 
	int colNum = X-> size2;
	ofstream ofs;
	ofs.open(outputFile.c_str());
	ofs << "#chro\tstart\end\t" << header << endl;
	for (int i = 0; i < rowNum; i++) {
		ofs << chro << "\t" << coord[i] << "\t" << coord[i]+binSize << "\t";
		for (int j = 0; j < colNum; j++) {
			if (map[i] >= 0) {
				ofs << X->data[map[i] * X->tda + j] << "\t";
			} else {
				ofs << "N/A\t";
			}
		}	
		ofs << endl;
	}
	ofs.close();
	return 0;
}

int io::write_significant_regions(string outputFile, string chro, vector<int>& coord, int binSize, bool writeToFile[], gsl_vector* metric, gsl_vector* pval, gsl_vector* qval, gsl_vector* reject_null, int map[], int n) {
	ofstream ofs;
	ofs.open(outputFile.c_str());
	ofs << "#chro\tstart\tend\t|diff|\tpval\tpadj" << endl;
	for (int i=0; i <n; i++) {
		int new_i = map[i];
		if (new_i >= 0) {
			if ((writeToFile[new_i]) && (gsl_vector_get(reject_null, new_i) == 1)) { 
				ofs << chro << "\t" << coord[i] << "\t" << coord[i]+binSize << "\t" << gsl_vector_get(metric, new_i) << "\t" << gsl_vector_get(pval, new_i) << "\t" << gsl_vector_get(qval, new_i) << "\t" <<  endl;
			}
		}
	}
	ofs.close();
	return 0;
}

int io::write_significant_regions_for_debug(string outputFile, string chro, vector<int>& coord, int binSize, bool writeToFile[], gsl_vector* metric, gsl_vector* pval, gsl_vector* qval, gsl_vector* reject_null, int map[], int n) {
	ofstream ofs;
	ofs.open(outputFile.c_str());
	for (int i=0; i <n; i++) {
		int new_i = map[i];
		ofs << chro << "\t" << coord[i] << "\t" << coord[i]+binSize << "\t";
		if (new_i >= 0) {
			ofs << gsl_vector_get(metric, new_i) << "\t" << gsl_vector_get(pval, new_i) << "\t" << gsl_vector_get(qval, new_i) << "\t" << writeToFile[new_i] << "\t" << (bool) gsl_vector_get(reject_null, new_i) << endl;
		} else {
			ofs << "N/A\tN/A\tN/A\tN/A\tN/A" << endl;
		}
	}
	ofs.close();
	return 0;
}
