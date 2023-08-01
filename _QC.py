import itertools
import re
import numpy as np
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import functools

import bioquest
from .utils import subset, mkdir
from typing import Tuple, Union, Optional


def qc_hist(
    adata: sc.AnnData,
    batch_key,
    n_genes_by_counts=(0, 4000),
    total_counts=(0, 4000),
    output_dir="./",
    prefix="",
    suffix="",
    dpi=300,
    formats: Union[str, Tuple[str, ...]] = ("pdf",),
) -> Optional[sc.AnnData]:
    """ """
    mkdir(output_dir)
    _export = bioquest.tl.export(
        formats=formats, od=output_dir, prefix=prefix, suffix=suffix, dpi=dpi
    )
    _adata = adata
    for x in np.unique(_adata.obs.loc[:, batch_key]):
        fig, axs = plt.subplots(1, 4, figsize=(12, 3))
        _adata2 = _adata[_adata.obs.loc[:, batch_key] == x]
        fig.suptitle(f"Covariates for filtering: {x}")
        sns.histplot(_adata2.obs.total_counts, kde=False, ax=axs[0])
        _a1 = bioquest.tl.subset(
            _adata2.obs,
            subsets={
                "n_genes_by_counts": f"{n_genes_by_counts[0]}<x<{n_genes_by_counts[1]}"
            },
        )
        sns.histplot(_a1.total_counts, kde=False, ax=axs[1])

        sns.histplot(_adata2.obs.n_genes_by_counts, kde=False, ax=axs[2])
        _a2 = bioquest.tl.subset(
            _adata2.obs,
            subsets={"total_counts": f"{total_counts[0]}<x<{total_counts[1]}"},
        )
        sns.histplot(_a2.n_genes_by_counts, kde=False, ax=axs[3])
        plt.subplots_adjust(wspace=0.5)
        _export(f"hist_for_{x}")


def qcMetrics(
    adata: sc.AnnData,
    *,
    batch_key: Optional[str] = None,
    output_dir: str = "./",
    min_genes: int = 0,
    min_cells: int = 0,
    mitochondrion: bool = False,
    hemoglobin: bool = False,
    ribosome: bool = False,
    percent_top=(50,),
    log1p: bool = True,
    prefix="",
    suffix="",
    dpi=300,
    inplace: bool = True,
    formats: Union[str, Tuple[str, ...]] = ("pdf", "png"),
) -> Optional[sc.AnnData]:
    """
    qcMetrics
    sk.qcMetrics(adata_spatial,batch_key="Sample",output_dir=OUTPUT_DIR,mitochondrion=True)
    """
    mkdir(output_dir)
    _export = bioquest.tl.export(
        formats=formats, od=output_dir, prefix=prefix, suffix=suffix, dpi=dpi
    )

    _adata = adata if inplace else adata.copy()

    if min_cells:
        sc.pp.filter_cells(_adata, min_genes=min_genes)
    if min_cells:
        sc.pp.filter_genes(_adata, min_cells=min_cells)
    if hemoglobin:
        _adata.var.loc[:, "Hb"] = _adata.var_names.isin(
            ("HBA1", "HBA2", "HBB", "HBD", "HBE1", "HBG1", "HBG2", "HBM", "HBQ1", "HBZ")
        )
    if mitochondrion:
        _adata.var.loc[:, "Mito"] = bioquest.st.detects(
            string=_adata.var_names, pattern=r"^mt-", flags=re.I
        )
    if ribosome:
        _adata.var.loc[:, "Ribo"] = _adata.var_names.str.startswith(r"^RP[SL0-9]")
    # ribo_url = "http://software.broadinstitute.org/gsea/msigdb/download_geneset.jsp?geneSetName=KEGG_RIBOSOME&fileType=txt"
    # ribo_genes = pd.read_table(ribo_url, skiprows=2, header = None)
    # _adata.var['Ribo'] = _adata.var_names.isin(tuple(ribo_genes[0].values))

    qc_vars = np.array(["Mito", "Ribo", "Hb"])[[mitochondrion, ribosome, hemoglobin]]
    if qc_vars:
        sc.pp.calculate_qc_metrics(
            _adata, qc_vars=qc_vars, percent_top=percent_top, log1p=log1p, inplace=True
        )

    keys = np.array(
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_Mito",
            "pct_counts_Ribo",
            "pct_counts_HB",
        ]
    )
    ks = keys[[True, True, mitochondrion, ribosome, hemoglobin]]
    _n = len(ks)
    # /* violin plot */
    _, axes = plt.subplots(1, _n, figsize=(3 * _n, 3))
    for a, k in enumerate(ks):
        sc.pl.violin(_adata, keys=k, jitter=False, show=False, ax=axes[a])
    plt.subplots_adjust(wspace=0.5)
    _export("QC_Violin_whole")

    # /* scatter plot */
    if qc_vars:
        iters = list(
            itertools.product(
                ks[[0, 1]], keys[[False, False, mitochondrion, ribosome, hemoglobin]]
            )
        )
        _, axes = plt.subplots(1, _n, figsize=(3 * _n, 3))
        for i, (x, y) in enumerate(iters):
            sc.pl.scatter(
                _adata,
                x=x,
                y=y,
                color=batch_key,
                ax=axes[i],
                legend_loc="none",
                show=False,
            )
        sc.pl.scatter(
            _adata,
            x=ks[0],
            y=ks[1],
            color=batch_key,
            ax=axes[-1],
            legend_loc="none",
            show=False,
        )
        plt.subplots_adjust(wspace=0.5)
        _export("QC_Scatter")
    else:
        sc.pl.scatter(
            _adata, x=ks[0], y=ks[1], color=batch_key, legend_loc="none", show=False
        )
        _export("QC_Scatter")

    # /* violin plot */
    _, axes = plt.subplots(1, _n, figsize=(6 * _n, 6))
    for a, k in enumerate(ks):
        sc.pl.violin(
            _adata,
            keys=k,
            groupby=batch_key,
            rotation=90,
            jitter=False,
            show=False,
            stripplot=False,
            ax=axes[a],
        )
    plt.subplots_adjust(wspace=0.5)
    _export("QC_Violin_single")

    return None if inplace else _adata


def single_qc_plot(
    adata,
    batch_key: str,
    output_dir="./",
    prefix="",
    suffix="",
    dpi=300,
    formats: Union[str, Tuple[str, ...]] = ("pdf",),
) -> Optional[sc.AnnData]:
    """ """
    mkdir(output_dir)
    _export = bioquest.tl.export(
        formats=formats, od=output_dir, prefix=prefix, suffix=suffix, dpi=dpi
    )
    ks = np.array(
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_Mito",
            "pct_counts_Ribo",
            "pct_counts_HB",
        ]
    )
    ks = np.intersect1d(ks, adata.obs.columns)
    _n = len(ks)
    for id in np.unique(adata.obs.loc[:, batch_key]):
        _adata = adata[adata.obs.loc[:, batch_key] == id]
        # /* violin plot */
        _, axes = plt.subplots(1, _n, figsize=(3 * _n, 3))
        for a, k in enumerate(ks):
            sc.pl.violin(_adata, keys=k, jitter=False, show=False, ax=axes[a])
        plt.subplots_adjust(wspace=0.5)
        _export(f"violin_for_{id}")

        # /* scatter plot */
        if _n > 2:
            iters = list(itertools.product(ks[[0, 1]], ks[2:]))
            _, axes = plt.subplots(1, _n, figsize=(3 * _n, 3))
            for i, (x, y) in enumerate(iters):
                sc.pl.scatter(
                    _adata,
                    x=x,
                    y=y,
                    color=batch_key,
                    ax=axes[i],
                    legend_loc="none",
                    show=False,
                )
            sc.pl.scatter(
                _adata,
                x=ks[0],
                y=ks[1],
                color=batch_key,
                ax=axes[-1],
                legend_loc="none",
                show=False,
            )
            plt.subplots_adjust(wspace=0.5)
            _export(f"scatter_for_{id}")
        else:
            sc.pl.scatter(
                _adata, x=ks[0], y=ks[1], color=batch_key, legend_loc="none", show=False
            )
            _export(f"scatter_for_{id}")


# def doubletMetrics(
#     adata:sc.AnnData,
#     *,
#     od:str='./',
#     prefix:str='',
#     suffix:str='',
#     dpi:int = 300,
#     formats:Union[str,Tuple[str,...]] = ('pdf','png'),
#     inplace:bool = True,
#     ) -> Optional[sc.AnnData]:
#     """
#     """
#     import doubletdetection
#     _adata = adata if inplace else adata.copy()
#     _export = functools.partial(sckit.tl.export,formats=formats,od=od,prefix=prefix,suffix=suffix,dpi=dpi)
#     clf = doubletdetection.BoostClassifier()
#     # raw_counts is a cells by genes count matrix
#     _adata.obs['doublet'] = clf.fit(_adata.X.copy()).predict()
#     # higher means more likely to be doublet
#     _adata.obs['doublet_score'] = clf.doublet_score()
#     doubletdetection.plot.convergence(clf,
#         save=od+'/QC_convergence.pdf',
#         show=False,
#         p_thresh=1e-16,
#         voter_thresh=0.5)
#     doubletdetection.plot.threshold(clf, save=od+'/QC_DoubletDetection_threshold.pdf', show=False, p_step=6)
#     sc.pl.violin(_adata,"doublet_score",show=False)
#     _export("QC_doublet_score")

#     return None if inplace else _adata

# # /* doublet */
# if doublet_method == 'scrublet':
#     sc.external.pp.scrublet(_adata, expected_doublet_rate = 0.05, threshold = 0.25,batch_key=batch_key)
#     sc.external.pl.scrublet_score_distribution(_adata,show=False)
#     _export("QC_scrublet_score_distribution")


def batch_subset(
    adata: sc.AnnData,
    batch_key: str,
    afilters: dict,
) -> Optional[sc.AnnData]:
    """ """
    _adata = adata
    temp_list = []
    for x in np.unique(_adata.obs.loc[:, batch_key]):
        _adata2 = _adata[_adata.obs.loc[:, batch_key] == x]
        temp_list.append(subset(_adata2, afilters.get(x), inplace=False))
    return sc.concat(temp_list)


def percent_subset(
    adata: sc.AnnData,
    batch_key: str,
    n_genes_by_counts: tuple = (2, 98),
    total_counts: tuple = (2, 98),
    pct_counts_Mito: float = 10.0,
) -> sc.AnnData:
    """
    percent_subset
    """
    _adata = adata
    temp_list = []
    for x in np.unique(_adata.obs.loc[:, batch_key]):
        _adata2 = _adata[_adata.obs.loc[:, batch_key] == x]
        gene_bottom, gene_up = np.percentile(
            _adata2.obs.n_genes_by_counts.values, n_genes_by_counts
        )
        count_bottom, count_up = np.percentile(
            _adata2.obs.total_counts.values, total_counts
        )
        afilter = {
            "pct_counts_Mito": f"x<{pct_counts_Mito}",
            "n_genes_by_counts": f"{gene_bottom}<x<{gene_up}",
            "total_counts": f"{count_bottom}<x<{count_up}",
        }
        temp_list.append(subset(_adata2, afilter, inplace=False))
    return sc.concat(temp_list)


def mad_filter(adata: sc.AnnData, *metric_nmad: tuple[str, int], **kwds) -> sc.AnnData:
    """
    use median ± n * mad as indicators to filter out outliers
    """
    batch_key = kwds.get("batch_key", None)

    def is_outlier(
        adata, metric: str, nmads: int, batch_key: Optional[str] = batch_key
    ):
        def helper(m: pd.Series, nmads: int = nmads):
            from scipy.stats import median_abs_deviation

            median_ = m.median()
            n_mad = median_abs_deviation(m) * nmads
            return np.logical_or(m < median_ - n_mad, m > median_ + n_mad)

        if batch_key:
            batches = adata.obs.loc[:, batch_key]
            return adata.obs.loc[:, metric].groupby(batches).apply(helper).droplevel(0)
        return helper(adata.obs.loc[:, metric])

    outliers = functools.reduce(
        np.logical_or, [[is_outlier(adata, x, y)] for x, y in metric_nmad]
    )
    return adata[~outliers]
