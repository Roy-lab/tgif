### Tree structure based on input matrix similarities

TGIF requires a tree file specifying the relationship among the input matrices and their corresponding conditions in a hierarchy, e.g. https://github.com/Roy-lab/tgif/blob/main/input/tgif-db/chr19_tree.txt. However, the structure of the tree or the relationship among the conditions may not be known in some situations, for instance among pseudo-bulk interaction profiles from single-cell HiC. 

This script builds a tree based on the similarity among the input matrix structures and outputs a file that can be used as direct input to TGIF. The similarity of each pair of input HiC matrix is measured with stratum-adjusted correlation coefficient (SCC), originally developed for HiCRep ([Yang et al., 2017](https://doi.org/10.1101/gr.220640.117)) to measure the reproducibility of HiC count matrices. SCC is converted to distance (1-SCC), then this distance matrix is used to build a hierarch using average linkage, just like one of the intermediate steps in hierarchical clustering. Instead of using the hierarchy to cluster, we simply use the structure as is as input to TGIF. Original SCC implementation for HiCRep was adopted for this script.

#### Usage
```
python suggestTree.py [file containing list of input matrices] [coordinates bed file] [output file prefix]
python suggestTree.py input/chr19_25kb_mat_files.txt input/chr19_25kb.bed output/chr19_
```

#### Input 1: file containing list of input matrices/conditions
- Column 1: ID of each condition/input matrix; start at 1
- Column 2: Alias for each condition (also used for output later)
- Column 3: path to each input matrix (it is a relative path here, but use absolute path)
```
1	ES	input/ES_chr19_25kb.txt
2	NPC	input/NPC_chr19_25kb.txt
3	CN	input/CN_chr19_25kb.txt
```

#### Input 2: bed file containing genomic coordinates
- See https://github.com/Roy-lab/tgif/tree/main#input-coordinates-file-format.

#### Final parameter: output file prefix
- Specify the path and any prefix to the output files

#### Output 1: tree file
- The constructed tree, written to file in a format required by TGIF
```
1	5	ES	input/ES_chr19_25kb.txt
2	4	NPC	input/NPC_chr19_25kb.txt
3	4	CN	input/CN_chr19_25kb.txt
4	5	branch1	N/A
5	-1	root	N/A
```

#### Output 2: a dendrogram plot
- A very simple plot of the constructed tree structure in pdf.

#### References

Yang T, Zhang F, Yardımcı GG, Song F, Hardison RC, Noble WS, Yue F, Li Q. HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient. Genome Res. 2017 Nov;27(11):1939-1949. doi: 10.1101/gr.220640.117. Epub 2017 Aug 30. PMID: 28855260; PMCID: PMC5668950.

