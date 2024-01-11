import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, Union, Optional
import pandas as pd
import scanpy as sc
import bioquest
from .utils import subset, mkdir

def cell_ratio(
    adata, 
    x:str, 
    y:str,
    *,
    palette = None,
    normalize = True,
    od=None,
    legend=True,
    figsize=(6, 3)
    ):

    df = adata.obs.loc[:,[x,y]]
    x_items = sorted(df[x].unique().tolist())
    y_items = sorted(df[y].unique().tolist())

    if palette is None:
        palette = dict(zip(y_items, adata.uns[f'{y}_colors']))

    heights = []
    for x_item in x_items:
        tmp_result = []
        x_item_counter = df[df[x]==x_item][y].value_counts().to_dict()
        for y_item in y_items:
            tmp_result.append(x_item_counter.get(y_item, 0))
        heights.append(tmp_result)
    heights = np.asarray(heights)

    if normalize:
        heights = heights/np.sum(heights, axis=0)
    heights = (heights.T/np.sum(heights, axis=1)).T

    plt.figure(figsize=figsize)
    _last = np.matrix([0.]* heights.shape[0])
    for i, y_item in enumerate(y_items):
        p = plt.bar(range(0, heights.shape[0]), heights[:, i],
                    bottom=np.asarray(_last)[0],
                    color=palette.get(y_item, 'b'),
                    label=y_item
                    )
        _last = _last + np.matrix(heights[:, i])
    plt.xticks(range(0, len(x_items)),labels=x_items,rotation=90)
    plt.ylim((0, 1))
    if legend:
        plt.legend()
        ax = plt.gca()
        ax.legend(bbox_to_anchor=(1.05, 1),loc='upper left', borderaxespad=0.)
    if od is not None:
        plt.savefig(od + f'/{x}_{y}_cell_ratio.pdf', dpi=300,bbox_inches="tight")
    return pd.DataFrame(heights,columns=y_items,index=x_items)

def plot_batch_effect(adata:sc.AnnData,
    cluster_key:str,
    batch_key:str = "SampleID",
    resolution: float = 1.0,
    use_rep:Optional[str] = None,
    output_dir:str = './',
    legend_loc:str = "right margin", # right margin, on data, 
    legend_fontsize:str = "small" , # [‘xx-small’, ‘x-small’, ‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’]
    n_pcs:int = 30,
    pca_use_hvg:bool = True,
    n_neighbors:int = 15,
    suffix:str = '',
    dpi:int = 300,
    formats:tuple = ('pdf','png')
    ) -> Optional[sc.AnnData] :
    mkdir(output_dir)
    _export = bioquest.tl.export(formats=formats,od=output_dir,suffix=suffix,dpi=dpi)
    sc.tl.pca(adata, svd_solver='arpack',n_comps=50,use_highly_variable=pca_use_hvg)
    sc.pl.pca_variance_ratio(adata,n_pcs=50,show=False)
    _export(f"PCA_variance_ratio")

    sc.pp.neighbors(adata, n_neighbors=n_neighbors,use_rep=use_rep, n_pcs=n_pcs)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color=batch_key,show=False,legend_loc=legend_loc,legend_fontsize=legend_fontsize)
    _export(f"UMAP_{batch_key}")
    
    sc.tl.leiden(adata, key_added=cluster_key, resolution=resolution)
    sc.pl.umap(adata, color=cluster_key,legend_loc="on data",size=20,show=False,legend_fontoutline=3);
    _export(f"UMAP_Cluster");