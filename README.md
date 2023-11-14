# Single_Cell_RNA_Seq
##Project overview:

The human body are highly organized system composed of many cells with incredible various types, states, and interactions to maintain cellular homeostasis. As almost all cells in humans and other organisms have the same set of genetic material, it remains to be further revisited to investigate the heterogeneity of the cellular environment. Whereas transcriptome (the entire collection of RNA sequences in a cell ) provides  the unique activity of only a subset of genes. Transcriptome allows researchers to profile the gene expression activity to probe cell identity, state, function, and response. Single-cell RNA-seq (scRNA-seq) offers a glimpse into analyzing the transcriptome at the single‐cell level. Here, mouse brain data were used from the Spateo tutorial to see the cell type, identification of marker genes, and their expression. Initial data analysis was done by Scanpy.

##Data sources:

For this project, I will be analyzing the mouse brain data from the Spateo tutorial. At first, data was loaded from Spateo and started analysis using Scanpy.

##Tools:

Scanpy(https://scanpy.readthedocs.io/en/stable/tutorials.html)

Spateo (https://spateo-release.readthedocs.io/en/latest/)

##Data analysis:

###Load and read the data:

                  import spateo as st
                  
                  fname_bin60 = "mousebrain_bin60.h5ad"
                  
                  adata = st.sample_data.mousebrain(fname_bin60)
                  
                  adata
                  
After reading in the data at Spateo, I'll perform preprocessing, basic filtering, cell quality control, and normalization to remove genes with uninformative and low quality of cells. Subsequently, I will do the clustering to visualize the cell population and find the marker.

###Preprocessing:

import numpy as np

import pandas as pd

import scanpy as sc

sc.pl.highest_expr_genes(adata, n_top=20, )

(It will show the gene with highest fraction of count in each cells)

###Basic filtering:

sc.pp.filter_cells(adata, min_genes=200)

sc.pp.filter_genes(adata, min_cells=3)

(It will show true cells that are high quality)

###Cell quality control:

sc.pl.violin(adata,)['n_genes_by_counts', 'total_conts', 'pct_conts_mt'],
             jitter=0.4, multi_panel=True)

(It will check very low gene count, high gene count, and high mitochondrial gene percentage. As high mitochondrial gene expression is indicative of cell stress, we need to consider for quality of data analysis).

###Normalization:

sc.pp.normalize_total(adata, target_sum=1e4)

sc.pp.log1p(adata)

(As each cell in scRNA-seq has a different number of reads, normalization helps to compare the expression between cells accurately. Here, I normalize the matrix to 10,000 reads per cell)

###UMAP and Leiden Clustering:

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

sc.tl.umap(adata)

sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata, color=['leiden'], legend_fontsize=8, save='_leiden')

(It will show to visualize the cell population)

###Marker gene:

sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

top_markers = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(5)

print(top_markers)

(It will allow us to perform a statistical test to find genes enriched in each cell population)

###Results:

####Visualization of cell populations

![Misti4](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/4266a9a1-70c7-4fa8-a09b-0aee77e71d43)

####Finding marker gene

![Misti2](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/5b5fed80-9f61-4166-96fd-567ba9a9d486)

####Visualization of gene expression in different cell types

![Misti3](https://github.com/Delowarrumana/Single_Cell_RNA_Seq/assets/146145134/4467bae0-aa3c-407c-9f75-e3312fb4771b)



