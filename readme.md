self-defining function for scRNA-seq, snRNA-seq and Spatial Transcriptomics.

```python
import sckit as sk

# quality control
sk.qcMetrics(adata,batch_key="library_id",output_dir=OUTPUT_DIR,mitochondrion=True,min_genes=300, min_cells=10)
sk.qc_hist(adata,batch_key="library_id",n_genes_by_counts=(4000,60000),total_counts=(40000,500000),output_dir=OUTPUT_DIR)
sk.single_qc_plot(adata,batch_key="library_id",output_dir=OUTPUT_DIR)
adata = sk.subset(adata,afilter,inplace=False)

# normalise
adata = sk.normalise(adata,batch_key='library_id',flavor="Seurat",n_top_genes=3000,output_dir=OUTPUT_DIR,pca_use_hvg=True,n_jobs=24, inplace=False)

# integrate
sk.harmony(adata,batch_key='library_id',output_dir=OUTPUT_DIR)

# cell type annotation
sk.label_helper(36)
sk.labeled(adata,cluster_names=new_cluster_names,reference_key='Cluster',cell_type_key='CellType')

# get differential expression gene dataframe
deg_df = sk.deg_df(adata).sort_values("logfoldchanges",ascending=False)
```