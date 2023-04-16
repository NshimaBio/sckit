# Gene-level analysis (DEGs, Enrichment, GRN)
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt

import pyscenic.aucell
from pyscenic.utils import GeneSignature

def aucell(adata:anndata.AnnData,score_name:str,gene_list:list,use_raw:bool=True):
    if use_raw:
        adata_raw = adata.raw.to_adata()
        exp_mtx=pd.DataFrame(adata_raw.X.toarray(),columns=adata_raw.var_names,index=adata_raw.obs_names)
    else:
        exp_mtx=pd.DataFrame(adata.X.toarray(),columns=adata.var_names,index=adata.obs_names)
    
    pd.DataFrame({score_name:gene_list}).to_csv("temp.csv",index=False)
    g=GeneSignature.from_grp("temp.csv",name=score_name)
    Path("temp.csv").unlink()
    auc_mtx = pyscenic.aucell.aucell(exp_mtx=exp_mtx, signatures=[g],noweights=True)
    adata.obs = pd.merge(adata.obs,auc_mtx,left_index=True,right_index=True)


def get_rank_array(adata,key,rank_name="rank_genes_groups"):
	deg_df = adata.uns[rank_name]
	return np.array(deg_df[key].tolist()).flatten()

def deg(adata,groupby,refKey="rest",pts=True,rank_name="rank_genes_groups",use_raw=True,method='wilcoxon'):
	adata.uns['log1p']["base"] = None
	sc.tl.rank_genes_groups(adata,groupby=groupby,pts=pts,reference=refKey,method=method,use_raw=use_raw,key_added=rank_name)

def deg_df(adata:anndata.AnnData,rank_name:str="rank_genes_groups",pts=True,use_raw=True):
	identy = list(adata.uns[rank_name]['names'].dtype.names)
	degs = pd.DataFrame({
		'Score':get_rank_array(adata,'scores',rank_name=rank_name),
		'Pvalue':get_rank_array(adata,'pvals',rank_name=rank_name),
		'Padj':get_rank_array(adata,'pvals_adj',rank_name=rank_name),
		'LogFC':get_rank_array(adata,'logfoldchanges',rank_name=rank_name),
		'Gene':get_rank_array(adata,'names',rank_name=rank_name)
		}
		)
	if use_raw:
		degs.insert(loc=0,value= identy * adata.raw.n_vars,column="Identy")
	else:
		degs.insert(loc=0,value= identy * adata.n_vars,column="Identy")
	if pts:
		pts_rest = adata.uns[rank_name]['pts_rest'].melt(ignore_index=False,value_name="PTS_Rest")
		pts_rest.drop(columns="variable",inplace=True)
		pts = adata.uns[rank_name]['pts'].melt(ignore_index=False,value_name="PTS")
		pts = pd.concat([pts,pts_rest],axis=1)
		pts.insert(loc=0,column="TEMP",value=pts.index.values + "XXX" + pts.variable.values)
		pts.reset_index(drop=True,inplace=True)

		degs.insert(loc=0,column="TEMP",value=degs.Gene.values + "XXX" + degs.Identy.values)
		degs = pd.merge(pts,degs,on="TEMP")
		degs.drop(columns=["variable","TEMP"],inplace=True)
	
	degs.set_index("Gene",drop=True,inplace=True)
	degs.sort_values(by="Identy",inplace=True)
	return degs