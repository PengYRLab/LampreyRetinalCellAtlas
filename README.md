# LampreyRetinalCellAtlas

## A framework for cross-species integration of scRNA-Seq of evolutionarily distant species

<img src="/Illustrative_Figures/Figure_cross-species_integration_framework_v0.6.png" width="400" height="600">
Figure. The illustration of the computational framework for cross-species integration of scRNA-Seq of evolutionarily distant species. (A) Orthologous gene assignment using OrthoFinder. Proteomes from different species are downloaded and preprocessed to retrieve the longest peptide for each gene. OrthoFinder is then used to resolve the gene tree and the orthogroup (a set orthologous genes). (B) Selection of the most variable orthologous genes from the orthogroup. First, the gene expression is log-normalized against the sequencing depth. Second, the standard deviation is calculated and used to rank the orthologous genes from the same orthogroup. The gene with the highest standard deviation is selected to assign the gene pair that is used for cross-species integration. (C) The CCA and rPCA approaches in Seurat and the criteria for choosing between CCA and rPCA. Seurat utilize mutual nearest anchors to remove cell type differences between species represented by different colors. However, CCA has tendency leading to over-integration when large proportions of cell types are not shared between datasets. We assess whether cells are over-integrated or not when using CCA by checking if cells from the same cell types fall apart in the integrated space. If so, we opt for rPCA for integration and fine-tuned parameters to enhance integration strength. (D) We also infer the cell-type correspondence using a machine learning technique. The XGBoost is used to train the model using the labeled cell types from one species. Then the model is applied to predict the cell types in another species.

## Directories 
### Functions
This directory includes functions written from scratch. 

### Tutorial
This directory includes subdirectoris with scripts used for the analyses. The scripts perform the analysis tasks.

**Cell type annotation using SingleR**

**Dimension reduction and clustering analysis using Seurat**

**Cross-species integration**

**Transcriptome mapping use XGBoost**

**Protein activity analysis using Rosa**

### Data
This directory includes data used for the analyses.

## How to cite
Please cite https://www.biorxiv.org/content/10.1101/2023.12.10.571000v1.abstract

## Code Author
Junqiang Wang

Email: junqiangwang333@gmail.com

@Peng Lab 
website: [](https://yirongpeng.com)https://yirongpeng.com
