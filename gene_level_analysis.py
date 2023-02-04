# Gene-level analysis (DEGs, Enrichment, GRN)

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

def get_rank_array(adata,key):
	rank_genes_groups = adata.uns['rank_genes_groups']
	return np.array(rank_genes_groups[key].tolist()).flatten()

def deg(adata,groupby,refKey,pvals_adj=1,pvals=1,log2FC=0,od='.'):
	adata.uns['log1p']["base"] = None
	sc.tl.rank_genes_groups(adata,groupby=groupby,reference=refKey,method='wilcoxon',use_raw=True)
	sc.pl.rank_genes_groups(adata,show=False)
	plt.savefig(od+'/rank_genes_groups.pdf')

def deg_df(adata,pvals_adj=1,pvals=1,log2FC=0):
	degs = pd.DataFrame({
		'names':get_rank_array(adata,'names'),
		'scores':get_rank_array(adata,'scores'),
		'pvals':get_rank_array(adata,'pvals'),
		'pvals_adj':get_rank_array(adata,'pvals_adj'),
		'logfoldchanges':get_rank_array(adata,'logfoldchanges')
	})
	filter = (degs.pvals_adj.values<pvals_adj) & (np.abs(degs.logfoldchanges.values)>log2FC) & (degs.pvals.values<pvals)
	return degs.loc[filter,:]