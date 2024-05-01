import pandas as pd
import numpy as np
import sys
from scipy.cluster import hierarchy
from scipy.sparse import coo_matrix
from hicrep import scc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

def readSparseMatrix(fname,n):
	df = pd.read_csv(fname, sep="\t",index_col=False, header=None, names=['i','j','cnt'])
	temp1 = df.loc[df.i != df.j]
	temp2 = pd.DataFrame({'i':temp1.j,'j':temp1.i,'cnt':temp1.cnt})
	df = pd.concat((df, temp2))
	X = coo_matrix((df.cnt,(df.i,df.j)),shape=(n,n))
	return X

def main(args):
	tableFile = args[1]
	bedFile = args[2]
	outputPrefix = args[3]

	bed = pd.read_csv(bedFile, sep="\t", index_col=False, header=None, names=['chrom','start','end','idx'])
	n = bed.shape[0]
	resolution = int(bed.end.values[0] - bed.start.values[0])

	df = pd.read_csv(tableFile, sep="\t", index_col=False, header=None, names=['id','alias','fname'])
	matrices = []
	for idx, row in df.iterrows():
		X = readSparseMatrix(row.fname, n)
		matrices.append(X)

	maxDiag = max(100, int(2000000/resolution))
	maxDiag = min(n, maxDiag)
	T = len(matrices)
	dist = []
	for i in range(T):
		for j in range(i+1,T):
			d = scc.sccByDiag(matrices[i],matrices[j],maxDiag) 
			dist.append(d)
	dist = 1 - np.array(dist)
	link = hierarchy.linkage(dist, method='average')
	root, nodes = hierarchy.to_tree(link, rd=True)
	hierarchy.dendrogram(link)
	plt.savefig(outputPrefix + 'dendrogram.pdf')

	nNodes = len(nodes)
	parent = np.full(nNodes, -2, dtype=int)
	for node in nodes:
		if (not node.is_leaf()):
			idx = node.get_id()
			left = node.get_left().get_id()
			right = node.get_right().get_id()
			parent[left] = idx
			parent[right] = idx
	
	aliases = df.alias.values.tolist()
	for i in range(nNodes-T-1):
		aliases.append('branch{:d}'.format(i+1))
	aliases.append('root')

	fileloc = df.fname.values.tolist()
	for i in range(nNodes-T):
		fileloc.append("N/A")

	dftree = pd.DataFrame({'id':np.arange(nNodes,dtype=int),'parent':parent,'alias':aliases,'fileloc':fileloc})
	dftree.id += 1
	dftree.parent += 1
	dftree.to_csv(outputPrefix + 'tree.txt',sep="\t",float_format="%d",index=False,header=False)
 
if __name__=="__main__":
	main(sys.argv)
