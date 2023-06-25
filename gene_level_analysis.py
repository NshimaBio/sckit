# Gene-level analysis (DEGs, Enrichment, GRN)
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt


def aucell(adata: sc.AnnData, gene_set: dict, use_raw: bool = True, seed: int = 1314):
    """
    AUCell uses the Area Under the Curve (AUC) to calculate whether a
    set of targets is enriched within the molecular readouts of each sample.
    To do so, AUCell first ranks the molecular features of each sample from
    highest to lowest value, resolving ties randomly.
    Then, an AUC can be calculated using by default the top 5% molecular
    features in the ranking. Therefore, this metric, aucell_estimate, represents
    the proportion of abundant molecular features in the target set, and their
    relative abundance value compared to the other features within the sample.
    """
    import decoupler

    net = (
        pd.DataFrame(dict([(k, pd.Series(v)) for k, v in gene_set.items()]))
        .melt(var_name="source", value_name="target")
        .dropna()
    )
    decoupler.run_aucell(adata, net=net, seed=seed, use_raw=use_raw)


def get_rank_array(adata, key, rank_name="rank_genes_groups"):
    deg_df = adata.uns[rank_name]
    return np.array(deg_df[key].tolist()).flatten()


def deg(
    adata,
    groupby,
    refKey="rest",
    pts=True,
    rank_name="rank_genes_groups",
    use_raw=True,
    method="wilcoxon",
    log_base=None,
):
    adata.uns["log1p"]["base"] = log_base
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        pts=pts,
        reference=refKey,
        method=method,
        use_raw=use_raw,
        key_added=rank_name,
    )


def deg_df(
    adata: sc.AnnData, rank_name: str = "rank_genes_groups", pts=True, use_raw=True
):
    identy = list(adata.uns[rank_name]["names"].dtype.names)
    degs = pd.DataFrame(
        {
            "Score": get_rank_array(adata, "scores", rank_name=rank_name),
            "Pvalue": get_rank_array(adata, "pvals", rank_name=rank_name),
            "Padj": get_rank_array(adata, "pvals_adj", rank_name=rank_name),
            "LogFC": get_rank_array(adata, "logfoldchanges", rank_name=rank_name),
            "Gene": get_rank_array(adata, "names", rank_name=rank_name),
        }
    )
    if use_raw:
        degs.insert(loc=0, value=identy * adata.raw.n_vars, column="Identy")
    else:
        degs.insert(loc=0, value=identy * adata.n_vars, column="Identy")
    if pts:
        pts_rest = adata.uns[rank_name]["pts_rest"].melt(
            ignore_index=False, value_name="PTS_Rest"
        )
        pts_rest.drop(columns="variable", inplace=True)
        pts = adata.uns[rank_name]["pts"].melt(ignore_index=False, value_name="PTS")
        pts = pd.concat([pts, pts_rest], axis=1)
        pts.insert(
            loc=0, column="TEMP", value=pts.index.values + "XXX" + pts.variable.values
        )
        pts.reset_index(drop=True, inplace=True)

        degs.insert(
            loc=0, column="TEMP", value=degs.Gene.values + "XXX" + degs.Identy.values
        )
        degs = pd.merge(pts, degs, on="TEMP")
        degs.drop(columns=["variable", "TEMP"], inplace=True)

    degs.set_index("Gene", drop=True, inplace=True)
    degs.sort_values(by="Identy", inplace=True)
    return degs


def deseq(
    adata: sc.AnnData, ref_level: str, design_factors: str, n_jobs: int = 16
) -> pd.DataFrame:
    """
    differential expression analysis (DEA) with scRNA-seq data
    adata: AnnData object with counts X matrix and obs
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    genes_to_keep = adata.X.sum(axis=0) >= 10
    adata = adata[:, ~genes_to_keep]

    count_df = adata.to_df()
    group_df = adata.obs
    count_df.index.name = None

    # 构建DeseqDataSet 对象
    dds = DeseqDataSet(
        # adata=adata,
        counts=count_df,
        clinical=group_df,
        ref_level=[design_factors, ref_level],
        design_factors=design_factors,
        refit_cooks=True,
        n_cpus=n_jobs,
    )
    # 离散度和log fold-change评估.
    dds.deseq2()
    # 差异表达统计检验分析
    stat_res = DeseqStats(
        dds, alpha=0.05, cooks_filter=True, independent_filter=True, n_cpus=n_jobs
    )
    stat_res.summary()
    return stat_res.results_df.rename(columns={"log2FoldChange": "LFC"})
