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

def auc_heatmap(adata,marker,out_prefix,ref_key="Cluster",figsize=(12,6),use_raw=True):
    import decoupler
    net=marker.melt(var_name="source",value_name="target").dropna()
    decoupler.run_aucell(adata,net,source="source",target="target",min_n=1,seed=1314,use_raw=use_raw)
    dt2=adata.obsm["aucell_estimate"].groupby(by=adata.obs.loc[:,ref_key]).agg(np.mean)
    import seaborn
    seaborn.clustermap(dt2.T,method='complete',z_score=0,cmap="viridis",figsize=figsize);
    plt.savefig(f"{out_prefix}.pdf",bbox_inches='tight')
    dt2.index.name="CellType"
    dt2.to_csv(f"{out_prefix}_score.csv.gz")

def score_heatmap(adata,marker_df,reference_key="Cluster",figsize=(9,6),return_score=False,save_fig=False):
    obs = adata.obs
    markers_dict = {x:np.intersect1d(marker_df.loc[:,x].dropna(),adata.raw.var_names) for x in  marker_df.columns}
    for x in markers_dict.keys():
        sc.tl.score_genes(adata,gene_list=markers_dict[x],score_name=f"{x}_Marker_Score")
    dt = bq.tl.select(adata.obs,columns=[reference_key],pattern="_Marker_Score$")
    adata.obs = obs
    a=dt.groupby(by=reference_key).apply(np.mean,axis=0)
    a.columns = bq.st.removes(string=a.columns,pattern=r"_Marker_Score$")
    import seaborn as sns
    sns.clustermap(a.T,method='complete',standard_scale=0,cmap="viridis",figsize=figsize);
    if return_score:
        return dt
    if save_fig:
        plt.savefig(f"{save_fig}/anno_heatmap.pdf",bbox_inches='tight')
    a.to_csv(f"{out_prefix}_score.csv.gz")