from typing import Tuple, Union, Optional
import scanpy as sc
from bioquest.tl import export

def normalise(adata:sc.AnnData,
	batch_key:str = "Sample",
	od:str = './',
	n_top_genes:int = 3000,
	target_sum:int = 10000,
	regress_out:str = ('pct_counts_Mito'),
	legend_loc:str = "on data", # right margin, on data, 
	legend_fontsize:str = "small" , # [‘xx-small’, ‘x-small’, ‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’]
	n_pcs:int = 30,
	n_neighbors:int = 15,
	n_jobs:int = 8,
	prefix:str = '',
	suffix:str = '',
	dpi:int = 300,
	flavor = 'seurat_v3',
	formats:tuple = ('pdf','png'),
	inplace:bool = True,
	) -> Optional[sc.AnnData] :
	"""
	normalise
	"""
	_export = export(formats=formats,od=od,prefix=prefix,suffix=suffix,dpi=dpi)

	_adata = adata if inplace else adata.copy()
	_adata.layers["counts"] = _adata.X.copy() # preserve counts

	sc.pp.normalize_total(_adata, target_sum=target_sum)
	sc.pp.log1p(_adata,base=None)

	if flavor == 'seurat_v3':
		span = 0.3
		error_status = True
		while error_status and span<0.8:
			try:
				sc.pp.highly_variable_genes(_adata, n_top_genes = n_top_genes,
				flavor="seurat_v3", subset=False, batch_key=batch_key,span=span,layer='counts')
				error_status = False
			except:
				span = span + 0.1
				print(f"highly_variable_genes erorr: add 0.1 to span, span={span} now.")
	
	if flavor == 'seurat':
		sc.pp.highly_variable_genes(_adata,
				flavor="seurat", subset=False)

	sc.pl.highly_variable_genes(_adata,show=False, log=True)
	_export("normalise" + "_highly_variable_genes_" + batch_key)
	
	_adata.raw = _adata
	_adata = _adata[:, _adata.var.highly_variable]
	
	if regress_out:
		sc.pp.regress_out(_adata, keys=regress_out, n_jobs=n_jobs)

	sc.pp.scale(_adata, max_value=10)

	sc.tl.pca(_adata, svd_solver='arpack',n_comps=50,use_highly_variable=True)
	sc.pl.pca(_adata, color=batch_key,legend_loc=legend_loc,legend_fontsize=legend_fontsize,show=False)

	_export("normalise" + "_PCA_" + batch_key)
	sc.pl.pca_variance_ratio(_adata,n_pcs=50,show=False)
	_export("normalise" + "_PCA_pca_variance_ratio_" + batch_key)

	sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
	sc.tl.umap(_adata)
	sc.pl.umap(_adata, color=batch_key,show=False,legend_loc=legend_loc,legend_fontsize=legend_fontsize)
	_export("normalise" + "_UMAP_" + batch_key)

	return None if inplace else _adata
