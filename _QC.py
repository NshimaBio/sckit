import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from bioquest.tl import export
from typing import Tuple, Union, Optional

def qcMetrics(adata:sc.AnnData,
    *,
    batch_key:Optional[str] = None,
    od:str='./',
    min_genes:int=200,
    min_cells:int=3,
    prefix='',
    suffix='',
    dpi = 300,
    inplace:bool = True,
    formats:Union[str,Tuple[str,...]] = ('pdf','png')
    ) -> Optional[sc.AnnData] :
    """
    QC plot: violin and scatter plot
    - doublet_method: one of None, scrublet, scvi, doubletdetection
    QCPlot(adata,min_cells=10,od='./')
    """
    _export = export(formats=formats,od=od,prefix=prefix,suffix=suffix,dpi=dpi)

    _adata = adata if inplace else adata.copy()

    if min_cells:
        sc.pp.filter_cells(_adata,min_genes=min_genes)
    if min_cells:
        sc.pp.filter_genes(_adata,min_cells=min_cells)

    _adata.var['HB'] = _adata.var_names.isin(("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ"))
    _adata.var['Mito'] = _adata.var_names.str.startswith(r'MT-',r"mt")
    _adata.var['Ribo'] = _adata.var_names.str.startswith(r'^RP[SL0-9]')
    # ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    # ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
    # _adata.var['Ribo'] = _adata.var_names.isin(tuple(ribo_genes[0].values))
    
    sc.pp.calculate_qc_metrics(_adata, qc_vars=['Mito',"Ribo","HB"], percent_top=None, log1p=False, inplace=True)

    ks = ('total_counts','n_genes_by_counts','pct_counts_Mito','pct_counts_Ribo','pct_counts_HB')

     # /* violin plot */
    _, axes = plt.subplots(3, 2, figsize=(5, 5))
    axes = axes.flatten()
    for a,k in enumerate(ks):
        sc.pl.violin(_adata, keys=k,jitter=False,show=False,ax=axes[a])
    plt.subplots_adjust(wspace = 0.3,hspace=0.3) 
    _export("QC_Violin_whole")
    
    # /* violin plot */
    _, axes = plt.subplots(2, 2, figsize=(8, 8))
    axes = axes.flatten()
    for a,k in enumerate(ks[1:]):
        sc.pl.scatter(_adata, x='total_counts', y=k,color=batch_key,ax=axes[a],legend_loc="none",show=False)
    plt.subplots_adjust(wspace = 0.3,hspace=0.3) 
    _export("QC_Scatter")

    # /* violin plot */
    _, axes = plt.subplots(3, 2, figsize=(30, 20))
    axes = axes.flatten()
    for a,k in enumerate(ks):
        sc.pl.violin(_adata,keys=k,groupby=batch_key,rotation=45,jitter=False,show=False,stripplot=False,ax=axes[a])
    plt.subplots_adjust(wspace = 0.2,hspace=0.3) 
    _export("QC_Violin_single")

    return None if inplace else _adata


def doubletMetrics(
    adata:sc.AnnData,
    *,
    od:str='./',
    prefix:str='',
    suffix:str='',
    dpi:int = 300,
    formats:Union[str,Tuple[str,...]] = ('pdf','png'),
    inplace:bool = True,
    ) -> Optional[sc.AnnData]:
    """
    """
    import doubletdetection
    _adata = adata if inplace else adata.copy()
    _export = functools.partial(sckit.tl.export,formats=formats,od=od,prefix=prefix,suffix=suffix,dpi=dpi)
    clf = doubletdetection.BoostClassifier()
    # raw_counts is a cells by genes count matrix
    _adata.obs['doublet'] = clf.fit(_adata.X.copy()).predict()
    # higher means more likely to be doublet
    _adata.obs['doublet_score'] = clf.doublet_score()
    doubletdetection.plot.convergence(clf, 
        save=od+'/QC_convergence.pdf', 
        show=False,
        p_thresh=1e-16, 
        voter_thresh=0.5)
    doubletdetection.plot.threshold(clf, save=od+'/QC_DoubletDetection_threshold.pdf', show=False, p_step=6)
    sc.pl.violin(_adata,"doublet_score",show=False)
    _export("QC_doublet_score")

    return None if inplace else _adata

    # # /* doublet */
    # if doublet_method == 'scrublet':
    #     sc.external.pp.scrublet(_adata, expected_doublet_rate = 0.05, threshold = 0.25,batch_key=batch_key)
    #     sc.external.pl.scrublet_score_distribution(_adata,show=False)
    #     _export("QC_scrublet_score_distribution")