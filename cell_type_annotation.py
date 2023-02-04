from io import StringIO
import pandas as pd
import scanpy as sc
from typing import Optional

def labeled(
	adata: sc.AnnData, 
	cluster_names: str, 
	reference_key: str, 
	cell_type_key: str = 'CellType', 
	inplace: bool = True
	) -> Optional[sc.AnnData]:
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
