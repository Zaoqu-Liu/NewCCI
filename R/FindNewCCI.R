#' @title Identifying New Cell-Cell Interactions Using Latent Ligands
#' @description This function performs a detailed analysis of cell-cell interactions (CCI) in single-cell RNA sequencing data using latent ligand-receptor pairs. The function integrates cell type information and differentially expressed genes to identify potential cell-cell interactions. It uses the CellChat database and the wilcoxauc method from Presto for differential expression analysis.
#' @author Zaoqu Liu; Email: liuzaoqu@163.com
#'
#' @param seu.obj A Seurat object containing the single-cell RNA-seq data.
#' @param latent.ligands A vector of latent ligand gene symbols that are hypothesized to be involved in cell-cell communication.
#' @param source.celltype The source cell type from which ligand signals are hypothesized to originate.
#' @param target.celltype The target cell type that receives signals via receptors.
#' @param celltype.colname The name of the column in `seu.obj@meta.data` that contains the cell type annotations. Default is 'Cell.Type'.
#' @param logFC.cutoff Numeric value specifying the log fold change threshold for filtering differentially expressed genes. Default is 0.2.
#' @param pct.dff.cutoff Numeric value for the percent difference cutoff between in-group and out-group cells. Default is 0.1.
#' @param padj.cutoff Numeric value specifying the adjusted p-value threshold for filtering significant interactions. Default is 0.05.
#' @param CCI.database The CellChat database to use for CCI analysis. Default is `TimeCCI::LRdf`.
#' @param species The species for the analysis. Default is "human".
#'
#' @return A CellChat object containing identified cell-cell interactions between the specified source and target cell types based on the given ligand-receptor pairs.
#'
#' @details
#' This function operates as follows:
#' 1. Differential expression analysis is performed on the target cell type using the Wilcoxon rank-sum test via the `presto::wilcoxauc()` function.
#' 2. Genes are filtered based on log fold change, percentage difference, and adjusted p-value cutoffs.
#' 3. Receptor genes are identified from the latent ligand-receptor pairs and filtered from the differentially expressed genes.
#' 4. A custom ligand-receptor interaction database is created using the provided latent ligands and the receptors identified from differentially expressed genes.
#' 5. The function creates a `CellChat` object from the Seurat object’s expression data and metadata.
#' 6. The custom ligand-receptor database is applied, and the CellChat pipeline is executed to identify overexpressed genes and interactions, compute communication probabilities, and aggregate interaction networks.
#'
#' The function returns a CellChat object that can be further used for visualizing and interpreting cell-cell communication.
#' @import presto
#' @import dplyr
#' @import CellChat
#' @export
FindNewCCI <- function(
    seu.obj,
    latent.ligands,
    source.celltype,
    target.celltype,
    celltype.colname = "Cell.Type",
    logFC.cutoff = 0.2,
    pct.dff.cutoff = 0.1,
    padj.cutoff = 0.05,
    CCI.database = TimeCCI::LRdf,
    species = "human") {
  # Perform differential expression analysis using the Wilcoxon rank-sum test
  dea <- presto::wilcoxauc(seu.obj)

  # Calculate the percentage difference and filter based on cutoffs
  dea2 <- dea %>%
    dplyr::mutate(pct.diff = pct_in - pct_out) %>%
    dplyr::filter(
      group == target.celltype,
      pct.diff > pct.dff.cutoff,
      logFC > logFC.cutoff,
      padj < padj.cutoff
    )

  # Identify receptor genes from the LR database
  receptors <- CCI.database$receptor.symbol
  dea2_receptor <- dea2[dea2$feature %in% receptors, ]

  # Create ligand-receptor pairs
  lr <- expand.grid(latent.ligands, dea2_receptor$feature)
  colnames(lr) <- c("ligand", "receptor")
  lr <- lr %>%
    dplyr::mutate(
      ligand = as.character(ligand),
      receptor = as.character(receptor),
      interaction_name = paste(ligand, receptor, sep = "_"),
      interaction_name2 = paste(ligand, receptor, sep = "-"),
      pathway_name = "Self"
    )

  # Update the CellChat database with the new ligand-receptor pairs
  require(CellChat)
  db.new <- CellChat::updateCellChatDB(db = lr, species_target = species)

  # Extract RNA counts and metadata
  x <- seu.obj@assays$RNA@counts
  data <- list(data = x, meta = seu.obj@meta.data)

  # Create a CellChat object
  cellchat <- CellChat::createCellChat(object = data$data, meta = data$meta, group.by = celltype.colname)

  # Assign the custom ligand-receptor database to the CellChat object
  cellchat@DB <- db.new

  # Subset data, identify overexpressed genes and interactions, and compute communication probabilities
  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
  cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE, trim = 0.1, nboot = 10, population.size = TRUE)
  cellchat <- CellChat::aggregateNet(cellchat)

  # Return the CellChat object for further analysis
  return(cellchat)
}
