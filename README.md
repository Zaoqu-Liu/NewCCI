### NewCCI: Identify Novel Cell-Cell Interactions

NewCCI is an R package designed to facilitate the identification of novel cell-cell interactions (CCI) using single-cell RNA sequencing (scRNA-seq) data. By allowing users to input hypothesized ligand candidates, NewCCI enhances traditional CCI analysis by discovering previously uncharacterized ligand-receptor pairs. The package integrates differential expression analysis with the CellChat database, providing a comprehensive approach to studying how cells communicate within a biological system.


### Installation

To install and set up all the required packages, and run the example analysis, use the following code:

```r
# Install required packages
install.packages(c("tidyr", "dplyr", "furrr", "future", "tibble"))
devtools::install_github("jinworks/CellChat")
remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
devtools::install_github("shenorrLabTRDF/cellAlign")
devtools::install_github("immunogenomics/presto")
devtools::install_github('Zaoqu-Liu/TimeCCI')
devtools::install_github('Zaoqu-Liu/NewCCI')
```

### FindNewCCI Tutorial

The core function, FindNewCCI, is used to explore potential cell-cell interactions between source and target cell types based on ligand-receptor interactions. NewCCI is ideal for researchers looking to extend their exploration of intercellular communication networks, particularly for latent or less-studied ligand candidates.

```r
# Clear the workspace and load the NewCCI package
rm(list = ls())
library(NewCCI)

# Run the FindNewCCI function
res <- FindNewCCI(
  seu.obj = TimeCCI::seu.example,
  latent.ligands = c('MMP2','MMP9','MMP11'),
  source.celltype = 'mCAF',
  target.celltype = 'EpiT',
  celltype.colname = 'Cell.Type'
)

# Visualize the cell-cell interactions
CellChat::netVisual_bubble(res, sources.use = 'mCAF', targets.use = 'EpiT')
```
