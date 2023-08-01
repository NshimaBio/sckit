from typing import Tuple, Union, Optional,List
import scanpy as sc
from .utils import subset, mkdir
import bioquest

def normalise(adata:sc.AnnData,
    batch_key:str = "SampleID",
    output_dir:str = './',
    n_jobs:int = 1,
    flavor:str = 'Seurat',
    n_top_genes:int = 3000,
    target_sum:int = 10000,
    scale:bool = False,
    regress_out:Optional[List[str]] = ['pct_counts_Mito'],
    prefix:str = '',
    suffix:str = '',
    dpi:int = 300,
    formats:tuple = ('pdf','png'),
    inplace:bool = False,
    ) -> Optional[sc.AnnData]:
    """
    normalise
    """
    mkdir(output_dir)
    _export = bioquest.tl.export(formats=formats,od=output_dir,prefix=prefix,suffix=suffix,dpi=dpi)
    _adata = adata if inplace else adata.copy()
    _adata.layers["counts"] = _adata.X.copy() # preserve counts
    sc.pp.normalize_total(_adata, target_sum=target_sum, inplace=True)
    sc.pp.log1p(_adata,base=None)
    if flavor == 'Seurat_v3':
        span = 0.3
        error_status = True
        while error_status and span<0.8:
            try:
                sc.pp.highly_variable_genes(_adata, n_top_genes = n_top_genes,
                flavor="seurat_v3", subset=False, batch_key=batch_key,span=span,layer='counts',inplace=True)
                error_status = False
            except:
                span = span + 0.1
                print(f"highly_variable_genes erorr: add 0.1 to span, span={span} now.")
    
    if flavor == 'Seurat':
        sc.pp.highly_variable_genes(_adata,n_top_genes=n_top_genes,flavor="seurat", batch_key=batch_key, subset=False,inplace=True)

    sc.pl.highly_variable_genes(_adata,show=False, log=True)
    _export("highly_variable_genes_" + batch_key)
    
    _adata.raw = _adata
    _adata = _adata[:, _adata.var.highly_variable]
    
    if regress_out:
        for i in regress_out:
            sc.pp.regress_out(_adata, keys=regress_out, n_jobs=n_jobs)
    if scale:
        sc.pp.scale(_adata, max_value=10)

    return None if inplace else _adata