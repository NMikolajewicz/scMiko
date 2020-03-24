
#' UMAP stratified by cluster ID
#'
#' UMAP plot with colors indicating cluster membership
#'
#' @param so Seurat Object
#' @param group.by Character. Metadata feature to group data by. Default is 'seurat_clusters'.
#' @param x.label Character. X axis label.
#' @param y.label Character. Y axis label.
#' @param plot.name Character. Plot title.
#' @param include.labels Logical specifying wheter to plot group IDs on UMAP.
#' @name cluster.UMAP
#' @return ggplot handle
#'
cluster.UMAP <- function(so, group.by = "seurat_clusters", x.label = "UMAP 1", y.label = "UMAP 2", plot.name = "UMAP", include.labels = T, ...){

  plt.handle <- DimPlot(so, reduction = "umap", group.by = group.by, label = include.labels, ...)  +
    ggtitle(label = plot.name) +
    xlab(x.label) + ylab(y.label)

  return(plt.handle)

}




#' Cell-level gene expression projected on UMAP
#'
#' UMAP plot with cell-level gene expression for queried gene.
#'
#' @param so Seurat Object
#' @param query.gene Character. Gene name to plot.
#' @param x.label Character. X axis label.
#' @param y.label Character. Y axis label.
#' @param plot.name Character. Plot title.
#' @name scExpression.UMAP
#' @return ggplot handle
#'
scExpression.UMAP <- function(so, query.gene, x.label = "UMAP 1", y.label = "", plot.name = NULL){

  if (is.null(plot.name)) plot.name <- query.gene

  plt.handle <- FeaturePlot(object = so, features = query.gene, cols =rev(brewer.pal(11,"RdYlBu")),
                            reduction = "umap",pt.size = TRUE, sort.cell = TRUE) +
    xlab(x.label) +
    ylab(y.label) +
    ggtitle(plot.name)

  return(plt.handle)

}


#' Plot variable genes in Seurat Object
#'
#' Using Seurat's VariableFeaturePlot, variable genes are identified and plotted. Top N variable genes are labelled.
#'
#' @param so Seurat Object
#' @param gNames Named gene list; entries are Symbols, names are Ensemble.
#' @param set_name Character specfiying name of dataset. Optional.
#' @param top.n.genes Numeric. Top n genes to label on plot.
#' @name variableGenes.Plot
#' @return ggplot handle
#'
variableGenes.Plot <- function(so, gNames, set_name = NULL, top.n.genes = 10){

  # top 10 most and least variable genes
  top10 <- head(VariableFeatures(so), assay = "gene_name", top.n.genes)

  # Plot variable features
  plt.handle <- VariableFeaturePlot(so)

  # plot title
  if (is.null(set_name)){
    plt.title <- "Variable Genes"
  } else {
    plt.title <- paste(set_name, ": Variable Genes")
  }

  # label genes
  plt.handle <- LabelPoints(plot = plt.handle,
                            points = top10,
                            labels = as.vector(gNames[top10]),
                            repel = TRUE) +
    ggtitle(label = plt.title)

  return(plt.handle)
}
