from typing import Tuple, Union, Optional
import scanpy as sc
import anndata
import bioquest

def normalise(adata:anndata.AnnData,
	batch_key:str = "Sample",
	output_dir:str = './',
	n_jobs:int = 1,
	flavor = 'Seurat',
	n_top_genes:int = 3000,
	target_sum:int = 10000,
	scale:bool = False,
	regress_out:str = ('pct_counts_Mito'),
	legend_loc:str = "right margin", # right margin, on data, 
	legend_fontsize:str = "small" , # [‘xx-small’, ‘x-small’, ‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’]
	n_pcs:int = 30,
	pca_use_hvg:bool = False,
	n_neighbors:int = 15,
	prefix:str = '',
	suffix:str = '',
	dpi:int = 300,
	formats:tuple = ('pdf','png'),
	inplace:bool = True,
	) -> Optional[anndata.AnnData] :
	"""
	normalise
	"""
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
		sc.pp.regress_out(_adata, keys=regress_out, n_jobs=n_jobs)

	if scale:
		sc.pp.scale(_adata, max_value=10)

	sc.tl.pca(_adata, svd_solver='arpack',n_comps=50,use_highly_variable=pca_use_hvg)
	sc.pl.pca(_adata, color=batch_key,legend_loc=legend_loc,legend_fontsize=legend_fontsize,show=False)
	_export("PCA_" + batch_key)
	sc.pl.pca_variance_ratio(_adata,n_pcs=50,show=False)
	_export("PCA_variance_ratio_" + batch_key)

	sc.pp.neighbors(_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
	sc.tl.umap(_adata)
	sc.pl.umap(_adata, color=batch_key,show=False,legend_loc=legend_loc,legend_fontsize=legend_fontsize)
	_export("UMAP_" + batch_key)

	return None if inplace else _adata
