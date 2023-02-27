### Tree-Guided Integrated Factorization (TGIF) to identify differential TAD boundaries

### [Step 1] Install 

Installation instructions below were tested in Linux Centos 7 distribution. [GSL (GNU Scientific Library) 2.7](https://www.gnu.org/software/gsl/doc/html/index.html) is used to handle matrix- and vector-related operations. Older versions of GSL (even 2.6) will not work since we use newly added functions. 

1. __If you already have GSL 2.7 installed__, edit the first few lines of the Makefile to point to the correct include and shared library directory, then jump to step 3.
```
#CHANGE PATHS AS NEEDED:
INCLUDE_PATH = ${CONDA_PREFIX}/include
LIBRARY_PATH = ${CONDA_PREFIX}/lib
```
2. __If you do not have GSL 2.7 installed, or you are not sure__, one way to get it installed is to use [conda](https://anaconda.org/conda-forge/gsl/):
```
conda install -c conda-forge gsl
```
3. Make sure to add the location of the installed shared library to where the compiler/linker will be looking. If you used conda to install GSL to the default location in step 2, run the following command after activating the appropriate environment (or add the appropriate path if you already have GSL installed):
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${CONDA_PREFIX}/lib
```
4. And let's install! In the same directory you downloaded the code/Makefile (either by cloning the repository or by downloading a release), run:
```
make tgif
```
5. If all went well, you won't get any alarming messages, and you will see an executable named `run_tgif` created in the same directory. A quick test below will print the manual for running TGIF:
```
./run_tgif -h
```

Note: in order to implement NNDSVD initialization of factors, a fast randomized SVD algorithm from [RSVDPACK](https://github.com/sergeyvoronin/LowRankMatrixDecompositionCodes) was used. A minor modification to allow random seed specification was made to the original code from [RSVDPACK](https://github.com/sergeyvoronin/LowRankMatrixDecompositionCodes/tree/main/single_core_gsl_code). This updated code is included under modules/random_svd directory. Compilation of this code is part of the included Makefile; no additional step is necessary for installation.

### [Step 2] Run

#### Basic usage
```
./run_tgif input/tree.txt input/coordinates.bed -o output/
```
- `input/tree.txt` specifies the tree file, which contains file locations to individual task matrices (paths are relative to location of run_tgif executable location; more details below). 
- `input/coordinates.bed` is a bed file specifying the start and end coordinates of each bin whose index is referred to in the sparse-format matrix files (more details below).
-	Optional: `-o output/` will put all output files to output/ directory. By default the output files will be saved to the current directory. Check out the example output directory in the repo.


#### Input tree file format
- This file specifies the task hierarchy, where each task represents some biological condition and has its own input matrix.
- See example in [input/tree.txt](https://github.com/Roy-lab/tgif/blob/main/input/tree.txt).
```
1	5	leaf1	input/matrix1.txt
2	5	leaf2	input/matrix2.txt
3	6	leaf3	input/matrix3.txt
4	6	leaf4	input/matrix4.txt
5	7	branch12	N/A
6	7	branch34	N/A
7	-1	root	N/A
```
All 4 columns are required. If you want a flat tree or do not care about the task hierarchy, point all nodes to have the root node as the parent.
- Column 1: **node ID**; start from 1 and move up.
- Column 2: **parent node ID**; right now the implementation will only work correctly if you ID all children nodes before a parent node (so start from the lowest depth of tree, move to next level, till you hit the root, which should be the last node ID.)
- Column 3: **node alias**, used to refer to each input matrix/task in the output files.
- Column 4: **location of input matrix file for leaf nodes**, relative to where the run_tgif executable is. Set as N/A for non-leaf nodes.

#### Input matrix file format
- Each of the matrix files referred to in the tree file.
- 0-indexed, tab-delimited, sparse-matrix format, no header
- See example in [input/matrix1.txt](https://github.com/Roy-lab/tgif/blob/main/input/matrix1.txt).
```
0	0	1000.2
0	1	1201.78
10	1	200.7
...
```

#### Input coordinates file format
- This file specifies the starting and ending coordinates of each bin whose index is referred to in the matrix files.
- The matrices are 0-indexed (i.e. first entry index is [0,0]) and they must share a common coordinates file.
- This coordinates bed file is tab-delimited.
- See exmaple in [input/coordinates.bed](https://github.com/Roy-lab/tgif/blob/main/input/coordinates.bed).
```
chr5	0	40000	0
chr5	40000	80000	1
chr5	80000	120000	2
chr5	120000	160000	3
...
```

#### Optional parameters
| Position | Parameter                    | Description                                                                                                                                                                                | Default Value/Behavior                    |
|:--------:|------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------|
|     1    | <input tree file>            | This file specifies the task hierarchy, where each task represents some biological condition and has its own input matrix. See  documentation above for more details, format, and example. | N/A                                       |
|     2    | <input coordinates bed file> | This is a bed file annotating each of the genomic regions with its task-specific boundary status. See documentation above for more  details, format, and example.                          | N/A                                       |
| optional | -o <output path and prefix>  | Output file path and prefix. Note: will NOT create a directory if the specified directory does not already exist.                                                                          | Output files saved to current  directory. |
| optional | -l                           | Generate a file called parameters.log under the specified output path (or current directory by default) that lists the parameter settings used.                                            | No log file of parameters saved.          |
| optional | -b                           | Generate files boundary_score.txt and background_score.txt under the  specified output path (or current directory by default).                                                             | No boundary-score related files  saved.   |
| optional | -p                           | Generate files boundary_pval.txt and boundary_adj_pval.txt under the  specified output path (or current directory by default).                                                             | No p-value related files saved.           |
| optional | -w <window size>             | Specify, in number of basepairs (e.g. 2000000), the size of the "window" or the dimension of the submatrices where factorization will be applied.                                          | 2000000 (i.e., 2MB)                       |
| optional | -s <step size>               | Specify, in number of basepairs (e.g. 1000000), the size of the step size for sliding down the block diagonal of the input matrix to generate each  set of submatrices.                    | 1000000 (i.e., 1MB)                       |
| optional | -a <alpha>                   | Specify an integer value used for the strength of the tree regularization.                                                                                                                 | 1000000                                   |

#### While running...
The program will print its execution progress while it runs and print the total time (in seconds) and max memory (in MB) consumption at the end. 
```
Reading in input files...
Performing TGIF...
[>                                                                     ] 1%
[===>                                                                  ] 5%
[======>                                                               ] 9%
[==========>                                                           ] 15%
[=====================>                                                ] 31%
[================================>                                     ] 47%
[============================================>                         ] 63%
[========================================================>             ] 80%
[===================================================================>  ] 96%
Testing for significant boundaries...
Total time elapsed: 213 seconds
Memory usage: 241MB
```

#### Output file 1: significant boundaries
- Each of the input matrix/task will have a column in the `significant_boundaries_summit_only.txt` file.
- Each row corresponds to a genomic region (i.e. 1st row below the header = 1st genomic region/bin in the coordinates bed file).
- 0 = non-boundary; 1 = significant boundary
- See example in [output/significant\_boundaries_summit_only.txt](https://github.com/Roy-lab/tgif/blob/main/output/significant_boundaries_summit_only.txt).
- NOTE: `significant_boundaries_summit_only.txt` marks only the region with the highest boundary score if a consecutive stretch of regions is flagged as significant boundaries. For the full set of such regions, see `significant_boundaries.txt` file which is also generated by default.
```
leaf1	leaf2	leaf3	leaf4
0	0	0	0	
0	0	0	0	
0	0	0	0	
0	0	0	1
...
```
In the example output shown above, the 4th bin/region in leaf4 task has a significant boundary.

#### Output file 2: annotations
- This is a bed file annotating each of the genomic regions with its task-specific boundary status.
- See example in [output/annotations.txt](https://github.com/Roy-lab/tgif/blob/main/output/annotations.txt).
```
chr5	0	40000	non-boundary
chr5	40000	80000	non-boundary
chr5	80000	120000	non-boundary
chr5	120000	160000	boundary specific to leaf4
```
In the example shown above, the genomic region chr5:120000-160000 contains a boundary specific to leaf4 task.

#### Optional output file: boundary scores
Using the `-b` flag, the user can print the boundary score for each region and for each task to `boundary_score.txt` file. Also, the scores used as the null distribution to test each boundary score's significance will be saved to `background_score.txt`. The format is the same as `significant_boundaries_summit_only.txt` file, but each entry is a boundary score instead.

#### Optional output file: boundary p-values
Using the `-p` flag, the user can print the p-values of each boundary score to `boundary_pval.txt` file. The adjusted p-values after FDR correction is printed to `boundary_adj_pval.txt`. The format is the same as `significant_boundaries_summit_only.txt` file, but each entry is a p-value instead.

#### Optional output file: parameter settings
Using the `-l` flag, the user can save the list of parameter values used in this run of TGIF to `parameters.log` file.
