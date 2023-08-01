import scanpy as sc
from .utils import subset, mkdir
from bioquest.tl import export

def scanorama(
    adata,
    batch_key:str,
    *,
    output_dir:str='./',
    n_neighbors:int=15,
    n_pcs:int =50,
    n_jobs:int=1,
    prefix:str='',
    suffix:str='',
    dpi:int=300,
    formats:str = 'pdf',
    legend_loc:str = "right margin", # right margin, on data, 
    legend_fontsize:str="x-small" , # [‘xx-small’, ‘x-small’, ‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’]
    inplace:bool = True
    ):
    """
    scanorama 
    adata = sk.scanorama(adata,batch_key="Sample",output_dir="output",legend_fontsize='x-small')
    """
    mkdir(output_dir)
    _adata = adata if inplace else adata.copy()
    _export = export(formats=formats,od=output_dir,prefix=prefix,suffix=suffix,dpi=dpi)
    sc.external.pp.scanorama_integrate(_adata,key=batch_key,adjusted_basis="X_scanorama")
    sc.pp.neighbors(_adata,n_neighbors=n_neighbors,use_rep="X_scanorama",n_pcs=n_pcs)
    sc.tl.umap(_adata)
    sc.pl.umap(_adata, color=batch_key,legend_loc=legend_loc,legend_fontsize=legend_fontsize,show=False)
    _export("UMAP_after_scanorama")
    sc.pl.embedding(_adata,color=batch_key,basis ='X_scanorama',legend_loc=legend_loc,legend_fontsize=legend_fontsize,show=False)
    _export("PCA_after_scanorama")
    sc.tl.tsne(_adata, use_rep="X_scanorama",n_jobs=n_jobs)
    sc.pl.tsne(_adata, color=batch_key,show=False);
    _export("TSNE_after_scanorama")
    return None if inplace else _adata