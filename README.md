# Single_Cell_RNA_Seq
##Project overview:

The intestine is the organ that plays an important role in nutrient digestion, absorption, defense from microorganisms, and hormone secretion. Here, I took the single-cell RNA-seq data of the epithelial cell from the human colon (Wang et al., 2018) and search for the following question:
•	Did each cluster exhibit highly expressed marker genes that matched with the provided cell-type annotation for each cluster?
•	Within the Paneth cell of the colon population, what specific sets of genes were overexpressed and under-expressed between males and females?

import scanpy as sc

import pandas as pd

import numpy as np

import decoupler

import sc_toolbox

import matplotlib.pyplot as plt

import seaborn as sns

import seaborn.objects as so

sc.logging.print_header()

sc.settings.set_figure_params(facecolor="white")

scanpy==1.9.5 anndata==0.10.2 umap==0.5.4 numpy==1.23.1 scipy==1.11.3 pandas==2.1.1 scikit-learn==1.3.2 statsmodels==0.14.0 igraph==0.10.8 louvain==0.8.0 pynndescent==0.5.10

adata = sc.read_h5ad(
    "C:/Users/mdhdu/AppData/Local/Programs/Python310/Scripts/cmd/data/epiithelial_colon.h5ad"
)

adata

AnnData object with n_obs × n_vars = 4329 × 17015

    obs: 'CellType', 'n_counts', 'log1p_n_counts', 'n_genes', 'log1p_n_genes', 'percent_mito', 'percent_ribo', 'percent_hb', 
    
    'percent_top50', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 
    
    'disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id',
    
    'sex_ontology_term_id', 'tissue_ontology_term_id', 'donor_id', 'suspension_type', 'cell_type', 'assay', 'disease',
    
    'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage'
    
    var: 'mito', 'ribo', 'hb', 'n_counts', 'n_cells', 'n_genes', 'highly_variable', 'means', 'dispersions', 
    
    'dispersions_norm', 'gene_symbols', 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype'

    uns: 'default_embedding', 'leiden', 'neighbors_hm', 'pca', 'schema_version', 'title'
    
    obsm: 'X_umap'
    
    varm: 'PCs'

Preprocessing

Before any analysis, we need to perform preprocessing to remove outliers and duplicates within the data

Remove cells with fewer than 200 genes and genes found in fewer than 3 cells

sc.pp.filter_cells(adata, min_genes=200)

sc.pp.filter_genes(adata, min_cells=3)

Visualize the top 20 genes with the highest expression across all cells. feature_id denotes the ID of each gene

sc.pl.highest_expr_genes(adata, n_top=20, gene_symbols="feature_name")

![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/10ee3ad5-b669-48a5-a82e-165021d21d5e)

This code snippet detects any genes with a MT- prefix (mitochondria) and plots statistics of these genes.

adata.var["mt"] = adata.var_names.str.startswith(
    "MT-"
)  # annotate the group of mitochondrial genes as 'mt'

sc.pp.calculate_qc_metrics(

    adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)
sc.pl.violin(

    adata,
    
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    
    jitter=0.4,
    
    multi_panel=True,
)
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/c06f47d6-d88a-4e95-9fe8-379e485822a4)
pct_counts_mt is strictly 0, meaning that our data does not have any fragment in the mitochondria.

We now visualize, for each cell, the correlation between the number of genes with non-zero expression (n_genes_by_counts)

vs the total number of fragments within each cell (total_counts).

sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/61927d38-8f79-41ed-b749-5319a7c70a4a)
We want to filter out some cells with too many total counts and number of expressed genes. Most tutorials suggest hard-coding a threshold, but the actual threshold is really dependent on the total number of fragments (read depth) of the experiment. Therefore, here we extract cells with n_genes_by_counts and total_counts from 2% quartile to 98% quartile

# Remove cells with too many total counts

upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.98)

lower_lim = np.quantile(adata.obs.n_genes_by_counts.values, 0.02)

print(f"{lower_lim} to {upper_lim}")

adata = adata[

    (adata.obs.n_genes_by_counts < upper_lim)
    
    & (adata.obs.n_genes_by_counts > lower_lim)
]

1224.56 to 4859.319999999999

Before we move on to normalizing expression, let's save the raw counts to a different layer for a downstream analysis.

adata.layers["counts"] = adata.X.copy()

C:\Users\mdhdu\AppData\Local\Temp\ipykernel_11188\1517723426.py:1: ImplicitModificationWarning: Setting element `.layers['counts']` of view, initializing view as actual.
  adata.layers["counts"] = adata.X.copy()

Normalization

There are many different normalization techniques, depending on the goal of downstream analyses. For example, the Single-cell Best Practices notebook presents shifted algorithm (divide cell counts by 10000 followed by log-transformed) (also present in Scanpy's tutorial), the Scran's pooling-based size factor estimation method, and the Analytic Pearson residuals method (which is also recommended in this online tutorial).

Since our goal is to identify differentially expressed genes for specific groups, shifted log transformation is sufficient.

# Normalize every cell to 10,000 UMI and log-transform

scales_counts = sc.pp.normalize_total(adata, target_sum=1e4, inplace=False)

adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

Feature selection

Even though we have filtered out genes that appear in fewer than 3 cells or cells with fewer than 200 genes, many of the remaining genes are uninformative, such as having zero expression across a range number of cells.

To select the most relevant genes for downstream analysis, we followed the recommended method proposed by Scanpy's tutorial: Identify highly variable genes with min_mean=0.0125, max_mean=3, and min_disp=0.5:

# Identify highly-variable genes
sc.pp.highly_variable_genes(

    adata, min_mean=0.0125, max_mean=3, min_disp=0.5, layer="log1p_norm"
)
sc.pl.highly_variable_genes(adata)
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/7064f991-fc4b-4835-baa3-fdf3c22527d2)
We can do a brief counting to see how many genes have been marked as highly variable (True means highly variable, False means otherwise). We can see we will retain 4523 genes for downstream analysis.

print(adata.var.highly_variable.value_counts())

adata = adata[:, adata.var.highly_variable]

highly_variable

False    14340

True      2675

Name: count, dtype: int64

We save the raw data for some downstream analysis:

adata.X = adata.layers["log1p_norm"]

adata.raw = adata

We then regress out effects of total counts per cell and the percentage of mitochondrial genes expressed and scale the data to unit variance. This code snippet may take some time to run

Dimensionality reduction and clustering

Dimensionality reduction with PCA

We would like to compress the data matrix into a smaller subspace while still retaining the most meaningful information. PCA is a widely accepted algorithm for this purpose. However, the number of principal components to retain, without keeping too much technical noise, is greatly dependent on the dataset itself (https://satijalab.org/seurat/articles/pbmc3k_tutorial).

(https://github.com/Mustardburger/BRN-assignments/blob/main/Single-cell_RNA-seq/Scanpy_Alzheimer's_reformatted.ipynb)

Without delving too much in this technical aspect, let's use the default value recommended by Scanpy: 50 principal components.

adata.X = adata.layers["log1p_norm"]

sc.tl.pca(adata, svd_solver="arpack")

adata.obsm["log1p_norm_pca"] = adata.obsm["X_pca"]

sc.pl.pca_variance_ratio(adata, log=True)
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/fe0ac09c-070d-4a9d-a289-14440424c47a)
Clustering with neighborhood graph

Our next step is to cluster the cells using a neighborhood graph on a PCA-transformed expression data. Normally, the clustering, followed by visualization by either tSNE or UMAP and identification of marker genes, is used for cell type assignment for each cluster. However, this dataset comes with readily available cell type labels. Our task is simply to confirm that the clustering separates the cell types well. We also want to see whether other attributes (such as development stage or Braak stage) are scattered or aggregated in distinct clusters.

Two main parameters for clustering are the following: n_neighbors which determines the number of neighboring data points, and n_pcs which determines the number of PCs from PCA to use. The choice of n_neighbors and n_pcs can lead to significant change in downstream analysis. Let's follow Scanpy's tutorial's recommended parameters:

sc.pp.neighbors(

    adata,
    
    n_neighbors=10,
    
    n_pcs=30,
    
    use_rep="log1p_norm_pca",
    
    key_added="log1p_norm_neighbors",
)
Visualization with UMAP

After clustering, we can run UMAP to visualize different clusters.

sc.tl.umap(adata, neighbors_key="log1p_norm_neighbors")

sc.pl.umap(adata)

C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:391: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap', 'norm' will be ignored
  cax = scatter(
  ![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/5bd0872f-edfe-47ec-9038-7b318008294a)
  Despite the difference, we can confirm that the clustering effectively captures different annotated cell types:

  sc.pl.umap(adata, color=["cell_type"])

  C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:1207: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
  if not is_categorical_dtype(values):
C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:1216: FutureWarning: The default value of 'ignore' for the `na_action` parameter in pandas.Categorical.map is deprecated and will be changed to 'None' in a future version. Please set na_action to the desired value to avoid seeing this warning
  color_vector = pd.Categorical(values.map(color_map))
C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:391: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
  ![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/4d45b23c-ab80-465f-8cdb-2bce5d0a278d)
  Besides, we can also look at how other attributes of the dataset appear on the cluster. For example, let's highlight sex:

  sc.pl.umap(adata, color=["sex"])

  C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:1207: FutureWarning: is_categorical_dtype is deprecated and will be removed in a future version. Use isinstance(dtype, CategoricalDtype) instead
  if not is_categorical_dtype(values):
C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:1216: FutureWarning: The default value of 'ignore' for the `na_action` parameter in pandas.Categorical.map is deprecated and will be changed to 'None' in a future version. Please set na_action to the desired value to avoid seeing this warning
  color_vector = pd.Categorical(values.map(color_map))
C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\scanpy\plotting\_tools\scatterplots.py:391: UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
  cax = scatter(
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/ac412fb1-1532-4cf9-941e-5ad00a7e15d9)
Find highly expressed genes in each cluster

Analysis across different cell types

After we have preprocessed the data, we would like to run differential gene analysis to identify what genes are most highly expressed in one particular cluster against the rest. The differential analysis can be run using our previous clustering methods: cell_type, and sex.

The results shown here use Scanpy's built-in sc.tl.rank_gene_groups function, which performs statistical testing between the cluster of interest and all other clusters combined. With single-cell analysis, gene expression can be very asymmetrical, so we use Wilcoxon ranked sum test for comparison.

Here, we run sc.tl.rank_gene_groups across cell_type clusters to determine the most differentially expressed genes between one cell_type cluster and the rest.

# This is used to find differentially expressed genes found in each cluster

# compared with the rest

sc.tl.rank_genes_groups(

    adata,
    
    "cell_type",
    
    method="wilcoxon",
    
    key_added="wilcoxon_cluster_cell_type",
)
sc.pl.rank_genes_groups(

    adata,
    
    n_genes=10,
    
    sharey=False,
    
    gene_symbols="feature_name",
    
    key="wilcoxon_cluster_cell_type",
)

C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\numpy\core\fromnumeric.py:84: FutureWarning: The behavior of DataFrame.sum with axis=None is deprecated, in a future version this will reduce over both axes and return a scalar. To retain the old behavior, pass axis=0 (or do not pass axis)
  return reduction(axis=axis, out=out, **passkwargs)
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/842275c1-9f8d-497c-9599-e8bfaf8b6e6c)
results = adata.uns["wilcoxon_cluster_cell_type"]

# Convert rank_gene_groups object to a dataframe

out = np.array([[0, 0, 0, 0, 0]])

for group in results["names"].dtype.names:

    out = np.vstack(
    
        (
            out,
            
            np.vstack(
            
                (
                    results["names"][group],
                    
                    results["scores"][group],
                    
                    results["pvals_adj"][group],
                    
                    results["logfoldchanges"][group],
                    
                    np.array([group] * len(results["names"][group])).astype("object"),
                    
                )
                
            ).T,
        )
    )

markers_df = pd.DataFrame(

    out[1:], columns=["Gene", "scores", "pval_adj", "lfc", "cluster"]
    
)

# Extract only differentially expressed genes

# p-value smaller than 0.05 and a log fold change larger than 1

markers_df = markers_df[(markers_df["pval_adj"] < 0.05) & (abs(markers_df["lfc"]) > 1)]

# Rename genes from its ID to its name

feature_reference_to_gene_name = (

    adata.var.reset_index()
    
    .loc[:, ["feature_reference", "feature_name"]]
    
    .rename(columns={"feature_reference": "Gene", "feature_name": "gene_name"})
)
markers_df = markers_df.merge(feature_reference_to_gene_name, on="Gene")

# Extract the top genes based on absolute log fold change

markers_df["neg_abs_scores"] = -markers_df["scores"].abs()

markers_df = markers_df.sort_values(by=["pval_adj", "neg_abs_scores"], ascending=True)

# markers_df.head()

# Take about 25 top genes

top_genes = markers_df.head(25)["gene_name"].tolist()

With a heatmap and some clustering of the rows, we show that specific cell clusters and their associated cell types highly expressed genes that were not highly expressed in any other cluster, demonstrating the genes' significance to that cluster/cell type.

sc.pl.rank_genes_groups_heatmap(
    adata,
    n_genes=5,
    key="wilcoxon_cluster_cell_type",
    groupby="cell_type",
    show_gene_labels=True,
    gene_symbols="feature_name",
)

WARNING: dendrogram data not found (using key=dendrogram_cell_type). Running `sc.tl.dendrogram` with default parameters. For fine tuning it is recommended to run `sc.tl.dendrogram` independently.
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/86a7d8bb-a127-4633-bf93-491786299ff7)
Using a dot plot, we can visualize whether the marker genes are being expressed in more than one cell type or not.

marker_genes = [
    "BTG2",
    "PHGR1",
    "TFF3",
    "CA7",
    "ADH1C",
    "ASCL2",
    "HMGB2",
]

sc.pl.dotplot(adata, marker_genes, groupby="cell_type", gene_symbols="feature_name")
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/09476bb4-98c7-4d5c-a433-c012b8edf6ba)
On the other hand, using the genes identified in rank_genes_group, we can see a totally different set of genes with high mean expression in each cell type.

sc.pl.rank_genes_groups_dotplot(adata, n_genes=5, key="wilcoxon", gene_symbols="feature_name", groupby="cell_type")
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/ba460ef3-0257-4a6d-b3d5-6ee0a631258d)
Differential expression (DE) across different conditions by pyDEseq2 on paneth cell of colon

After all the preprocessing, we are now interested in running differential expression (DE) analysis in oligodendrocyte cell populations between three Braak stages. For this analysis, we convert the dataset into pseudobulk for bulk differential expression analysis tools such as DEseq2. In particular, we extract cells belonging to a particular cell type, then take the sum of raw counts of genes across all cells within that cluster. We can now identify differentially expressed genes between two conditions.

We extract all paneth cell of colon in our dataset

adata.X = adata.layers["counts"]

adata_colon = adata[adata.obs["cell_type"] == "paneth cell of colon"]

adata_colon

View of AnnData object with n_obs × n_vars = 337 × 2675
    obs: 'CellType', 'n_counts', 'log1p_n_counts', 'n_genes', 'log1p_n_genes', 'percent_mito', 'percent_ribo', 'percent_hb', 'percent_top50', 'assay_ontology_term_id', 'cell_type_ontology_term_id', 'development_stage_ontology_term_id', 'disease_ontology_term_id', 'self_reported_ethnicity_ontology_term_id', 'is_primary_data', 'organism_ontology_term_id', 'sex_ontology_term_id', 'tissue_ontology_term_id', 'donor_id', 'suspension_type', 'cell_type', 'assay', 'disease', 'organism', 'sex', 'tissue', 'self_reported_ethnicity', 'development_stage', 'n_genes_by_counts', 'total_counts', 'total_counts_mt', 'pct_counts_mt'
    var: 'mito', 'ribo', 'hb', 'n_counts', 'n_cells', 'n_genes', 'highly_variable', 'means', 'dispersions', 'dispersions_norm', 'gene_symbols', 'feature_is_filtered', 'feature_name', 'feature_reference', 'feature_biotype', 'mt', 'n_cells_by_counts', 'mean_counts', 'pct_dropout_by_counts', 'total_counts'
    uns: 'default_embedding', 'leiden', 'neighbors_hm', 'pca', 'schema_version', 'title', 'hvg', 'log1p_norm_neighbors', 'umap', 'cell_type_colors', 'development_stage_colors', 'sex_colors', 'wilcoxon_cluster_cell_type', 'dendrogram_cell_type'
    obsm: 'X_umap', 'X_pca', 'log1p_norm_pca'
    varm: 'PCs'
    layers: 'counts', 'log1p_norm'
    obsp: 'log1p_norm_neighbors_distances', 'log1p_norm_neighbors_connectivities'

To increase statistical power, we want to include pseudoreplicates for each condition. The following code snippet creates pseudobulk for each Braak stage by summing counts of each gene across cells. Besides, for each Braak stage, it also creates 3 pseudoreplicates by randomly selecting a subset of cells in that group.

Code is adopted from (https://github.com/mousepixels/sanbomics_scripts/blob/main/pseudobulk_pyDeseq2.ipynb)
(https://github.com/Mustardburger/BRN-assignments/blob/main/Single-cell_RNA-seq/Scanpy_Alzheimer's_reformatted.ipynb)

import random

pseudobulk_list = []

num_pseudorep = 3

for sample in adata_colon.obs["sex"].unique():

    samp_cell_subset = adata_colon[adata_colon.obs["sex"] == sample]

    # Use raw count data
    
    # print(samp_cell_subset.layers['counts'])
    
    samp_cell_subset.X = samp_cell_subset.layers["counts"]

    # Create random indices to split the dataset
    
    indices = list(samp_cell_subset.obs_names)
    
    random.shuffle(indices)
    
    indices = np.array_split(np.array(indices), num_pseudorep)

    for i, pseudo_rep in enumerate(indices):
    
        # Create a pseudoreplicate object by summing gene counts across cells
        
        rep_adata = sc.AnnData(
        
            X=samp_cell_subset[indices[i]].X.sum(axis=0),
            
            var=samp_cell_subset[indices[i]].var[[]],
        )

        rep_adata.obs_names = [sample + "_" + str(i)]
        
        rep_adata.obs["condition"] = samp_cell_subset.obs["sex"].iloc[0]
        
        rep_adata.obs["replicate"] = i

        pseudobulk_list.append(rep_adata)

        Then create a new Scanpy object with obs.condition denoting different Braak stage and obs.replicate indicating its 
        replicate

        pseudobulk_adata = sc.concat(pseudobulk_list)
        
        pseudobulk_adata.obs.head()
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/322363cb-b1c2-4834-b54d-e6cc590ae2a3)
Now ready to run deseq2. The following code snippet runs deseq2 with our newly created pseudobulk object pseudobulk_adata. Explanations are included with comments in each code block.

from pydeseq2.dds import DeseqDataSet

from pydeseq2.ds import DeseqStats

# Extract counts into a dataframe

counts = pd.DataFrame(pseudobulk_adata.X, columns=pseudobulk_adata.var_names)

counts.head()
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/54cd1b2f-d4ed-4b32-b44e-38e3b85cd430)
# Create DEseqDataSet object

dds = DeseqDataSet(

    counts=counts, metadata=pseudobulk_adata.obs, design_factors="condition"
)

# Filter out genes with no expression in any cell

sc.pp.filter_genes(dds, min_cells=1)

# Run DESeq2

dds.deseq2()

Fitting size factors...

... done in 0.00 seconds.

Fitting dispersions...

... done in 0.89 seconds.

Fitting dispersion trend curve...

... done in 0.82 seconds.

C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\anndata\_core\views.py:144: RuntimeWarning: invalid value encountered in log

  results = super().__array_ufunc__(
  
Fitting MAP dispersions...

... done in 2.31 seconds.

C:\Users\mdhdu\AppData\Roaming\Python\Python310\site-packages\pydeseq2\dds.py:695: RuntimeWarning: invalid value encountered in log

  self.varm["_outlier_genes"] = np.log(self.varm["genewise_dispersions"]) > np.log(
Fitting LFCs...

... done in 0.95 seconds.

Refitting 0 outliers.

visualize whether the conditions are clustered well using PCA

sc.tl.pca(dds)

sc.pl.pca(dds, color="condition", size=200, annotate_var_explained=True)
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/c8a779e6-3014-43e5-bbda-51986354a2c3)
Thereafter  calculate the p-value and log2 fold change of differential gene expressions

# sex m vs f
stat_res_mvsf = DeseqStats(dds, n_cpus=8, contrast=("condition", "male", "female"))

stat_res_mvsf.summary()

de_mvsf = stat_res_mvsf.results_df

de_mvsf = de_mvsf.sort_values("stat", ascending=False).dropna()

Running Wald tests...

Log2 fold change & Wald test p-value: condition male vs female

... done in 9.24 seconds.
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/34172280-80f8-49f2-a25c-1e920a32780a)
Visualization

Now visualize the set of genes that are differentially expressed across these conditions

genes1 = sc.get.rank_genes_groups_df(adata_colon, group='male', key='wilcoxon')['names'][:20]

genes2 = sc.get.rank_genes_groups_df(adata_colon, group='female', key='wilcoxon')['names'][:20]

genes = genes1.tolist() +  genes2.tolist() 

sc.pl.dotplot(adata_colon,genes, groupby='sex')
![1](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/df31af1e-8672-4c4f-9bdb-7ae6c417406a)

















