import scanpy as sc
from bioquest.tl import export

def harmony(
	adata,
	batch_key:str,
	*,
	od:str='./',
	max_iter_harmony:int =20,
	n_pcs:int =20,
	n_neighbors:int=15,
	prefix:str='',
	suffix:str='',
	dpi:int=300,
	formats:str = 'pdf',
	legend_loc:str = "right margin", # right margin, on data, 
	legend_fontsize:str="x-small" , # [‘xx-small’, ‘x-small’, ‘small’, ‘medium’, ‘large’, ‘x-large’, ‘xx-large’]
	inplace:bool = True
	):
	"""
	harmony 
	adata = sk.pp.harmony(adata,batch_key="sample",od="results",legend_fontsize='x-small')
	"""
	_adata = adata if inplace else adata.copy()
	_export = export(formats=formats,od=od,prefix=prefix,suffix=suffix,dpi=dpi)
	sc.external.pp.harmony_integrate(_adata,key=batch_key,max_iter_harmony=max_iter_harmony,adjusted_basis='X_pca')
	sc.pp.neighbors(_adata,n_neighbors=n_neighbors)
	sc.tl.umap(_adata)
	sc.pl.umap(_adata, color=batch_key,legend_loc=legend_loc,legend_fontsize=legend_fontsize,show=False)
	_export("UMAP_after_harmony")
	sc.pl.embedding(_adata,color=batch_key,basis ='X_pca',legend_loc=legend_loc,legend_fontsize=legend_fontsize,show=False)
	_export("PCA_after_harmony")

	return None if inplace else _adata