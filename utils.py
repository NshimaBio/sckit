import scanpy as sc
from typing import Optional, Union

def subset(
        adata: sc.AnnData,
        subsets: dict,
        inplace: bool = False
) -> Optional[sc.AnnData]:
    """
    filter/subset a AnnData according to subsets conditions
    """
    _a = adata if inplace else adata.copy()
    for k in subsets:
        v = subsets.get(k)
        if isinstance(v, list):
            _lg = _a.obs[k].isin(v)
            _a = _a[_lg, :]
        else:
            _lg = _a.obs[k].apply(lambda x: eval(v))
            _a = _a[_lg, :]
    return None if inplace else _a
