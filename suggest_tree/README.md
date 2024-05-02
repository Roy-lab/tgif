### Tree structure based on input matrix similarities

TGIF requires a tree file specifying the relationship among the input matrices and their corresponding conditions in a hierarchy, e.g. https://github.com/Roy-lab/tgif/blob/main/input/tgif-db/chr19_tree.txt. However, the structure of the tree or the relationship among the conditions may not be known in some situations, for instance among pseudo-bulk interaction profiles from single-cell HiC. This script builds a tree based on the similarity among the input matrix structures and outputs a file that can be used as direct input to TGIF. The similarity of each pair of input HiC matrix is measured with stratum-adjusted correlation coefficient (SCC), originally developed for HiCRep to measure the reproducibility of HiC count matrices. SCC is converted to distance (1-SCC), then this distance matrix is used to build a hierarch using average linkage, just like one of the intermediate steps in hierarchical clustering. Instead of using the hierarchy to cluster, we simply use the structure as is as input to TGIF.

#### Usage
"""
python suggestTree.py

#### Input 1: table or list of input matrices and conditions

#### Input 2: 

