from io import StringIO
import numpy as np
import pandas as pd
import scanpy as sc
import anndata
from typing import Optional
import matplotlib.pyplot as plt
import bioquest as bq

def labeled(
	adata: anndata.AnnData, 
	cluster_names: str, 
	reference_key: str, 
	cell_type_key: str = 'CellType', 
	inplace: bool = True
	) -> Optional[anndata.AnnData]:
	'''
	labeled(hdata,cluster_names=new_cluster_names,reference_key='leiden',cell_type_key='CellType')
	'''

	_adata = adata if inplace else adata.copy()
	_ref_df = _adata.obs.loc[:, [reference_key]]
	_annot_df = pd.read_csv(StringIO(cluster_names), header=None, dtype='object')
	_adata.obs[cell_type_key] = pd.merge(
		_ref_df, _annot_df, left_on=reference_key, right_on=0, how='left')[1].values
	return None if inplace else _adata


def label_helper(number_of_cluster: int):
	'''
	number_of_cluster: 最后一个cluster的数字
	'''
	_s1 = ",\n".join([str(i) for i in range(number_of_cluster+1)])
	_s2 = "\nnew_cluster_names ='''\n" + _s1 + ",\n'''\n"
	print(_s2)

def anno_heatmap(adata,marker_df,reference_key="Cluster",figsize=(18,6),return_score=False,save_fig=False):
    obs = adata.obs
    markers_dict = {x:np.intersect1d(marker_df.loc[:,x].dropna().values,adata.raw.to_adata().var_names) for x in  marker_df.columns}
    for x in markers_dict.keys():
        sc.tl.score_genes(adata,gene_list=markers_dict[x],score_name=f"{x}_score")
    dt = bq.tl.select(adata.obs,columns=[reference_key],pattern="_score$")
    adata.obs = obs
    dt=dt.loc[adata.obs.loc[:,reference_key].sort_values().index,:]
    a=dt.groupby(by=reference_key).apply(np.median,axis=0)
    dt2 = pd.DataFrame({x:y for x,y in enumerate(a)},index=bq.st.removes(string=dt.columns[1:],pattern=r"_score$"))
    import seaborn as sns
    sns.clustermap(dt2,method='complete',standard_scale=True,cmap="viridis",figsize=figsize);
    if return_score:
        return dt2
    if save_fig:
        plt.savefig(f"{save_fig}/anno_heatmap.pdf",bbox_inches='tight')