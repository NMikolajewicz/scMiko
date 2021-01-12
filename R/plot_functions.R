
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
#' @param reduction Character specifying which dimensional reduction to use (e.g., umap, pca). Default is 'umap'.
#' @param ... additional arguments passed to Seurat::DimPlot().
#' @name cluster.UMAP
#' @return ggplot handle
#'
cluster.UMAP <- function(so, group.by = "seurat_clusters", x.label = "UMAP 1", y.label = "UMAP 2", plot.name = "UMAP", include.labels = T, reduction = "umap", ...){

  if (group.by == "seurat_clusters") so@meta.data[["seurat_clusters"]] <- orderedFactor(so@meta.data[["seurat_clusters"]])

  plt.handle <- DimPlot(so, group.by = group.by, label = include.labels,reduction = reduction, ...)  +
    ggtitle(label = plot.name) +
    xlab(x.label) + ylab(y.label)

  return(plt.handle)

}




#' Cell-level gene expression projected on UMAP
#'
#' UMAP plot with cell-level gene expression for queried gene. Uses Seurat::FeaturePlot().
#'
#' @param object Seurat Object
#' @param query.gene Character. Gene name to plot.
#' @param x.label Character. X axis label.
#' @param y.label Character. Y axis label.
#' @param adjust.pt.size Adjust point size for plotting. Logical.
#' @param order.cells Plot cells in order of expression. Logical.
#' @param plot.name Character. Plot title.
#' @param ... additional parameters passed to Seurat::FeaturePlot().
#' @name scExpression.UMAP
#' @return ggplot handle
#'
scExpression.UMAP <- function(object, query.gene, x.label = "UMAP 1", y.label = "", adjust.pt.size = T, order.cells = T, plot.name = NULL, ...){

  if (is.null(plot.name)) plot.name <- query.gene

  plt.handle <- FeaturePlot(object = object, features = query.gene, cols =rev(brewer.pal(11,"RdYlBu")),
                            reduction = "umap",pt.size = adjust.pt.size, sort.cell = order.cells, ...) +
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


#' QC violin plots
#'
#' Number of genes/cell, UMI/cell and mitochondiral content/cell are visualized with violin plots. Two ggplot handles are generated. First contains QC metrics pooled across all cells, while second stratifies dataset by barcode labels.
#'
#' @param so Seurat Object
#' @param plt.log.flag Logical specifying whether data are plotted on log scale. Default is True.
#' @name QC.violinPlot
#' @return list of ggplot handles
#'
QC.violinPlot <- function(so, plt.log.flag = T){

  # Overall and cell-type specific QC parameters
  if (plt.log.flag){
    # logarithmic scale
    plt1 <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = TRUE)
    plt2 <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "subset_group", ncol = 3, log = TRUE)
  } else{
    # raw scale
    plt1 <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log = FALSE)
    plt2 <- VlnPlot(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "subset_group", ncol = 3, log = FALSE)
  }

  output <- list(plt1, plt2)
  return(output)
}

#' QC scatter plots
#'
#' Pairwise relationships between 1) UMI count/cell and mitchondrial content, 2) gene count/cell and mitochondrial content, and 3) UMI count/cell and gene count/cell. Cellular density is computed and overlayed as color on scatter plot.
#'
#' @param so Seurat Object
#' @param legend.width numeric specifying legend width (in cm)
#' @param ... additional parameters passed to scMiko::getDensity().
#' @name QC.scatterPlot
#' @return ggplot handle
#'
QC.scatterPlot <- function(so, legend.width = 1, ...){

  df.meta <- so@meta.data

  # correlations
  rho1p <- signif(cor(x = df.meta$nCount_RNA, y =  df.meta$percent.mt, method = "pearson"), 2)
  rho1s <- signif(cor(x = df.meta$nCount_RNA, y =  df.meta$percent.mt, method = "spearman"), 2)
  rho2p <- signif(cor(x = df.meta$nFeature_RNA, y =  df.meta$percent.mt, method = "pearson"), 2)
  rho2s <- signif(cor(x = df.meta$nFeature_RNA, y =  df.meta$percent.mt, method = "spearman"), 2)
  rho3p <- signif(cor(x = df.meta$nCount_RNA, y =  df.meta$nFeature_RNA, method = "pearson"), 2)
  rho3s <- signif(cor(x = df.meta$nCount_RNA, y =  df.meta$nFeature_RNA, method = "spearman"), 2)

  # cell densities
  df.meta$density1 <- tryCatch({getDensity(df.meta$nCount_RNA, df.meta$percent.mt, n = 100, ...)},
                               error = function(e) {return(1)})
  df.meta$density2 <- tryCatch({getDensity(df.meta$nFeature_RNA, df.meta$percent.mt, n = 100, ...)},
                               error = function(e) {return(1)})
  df.meta$density3 <- tryCatch({getDensity(df.meta$nCount_RNA, df.meta$nFeature_RNA, n = 100, ...)},
                               error = function(e) {return(1)})

  # construct plots
  plt.handle1 <- df.meta %>% ggplot(aes(x = nCount_RNA, y = percent.mt, color = density1)) + geom_point() + theme_miko(legend = T) +
    xlab("UMI/cell") + ylab("Mitochondrial Content (%)") +
    labs(title = paste0("r = ", rho1p, "; rho = ", rho1s)) + scale_color_viridis("Density") +
    theme(legend.position="bottom", legend.key.width=unit(legend.width,"cm"))

  plt.handle2 <- df.meta %>% ggplot(aes(x = nFeature_RNA, y = percent.mt, color = density1)) + geom_point() + theme_miko(legend = T) +
    xlab("Genes/cell") + ylab("Mitochondrial Content (%)")  +
    labs(title = paste0("r = ", rho2p, "; rho = ", rho2s)) + scale_color_viridis("Density") +
    theme(legend.position="bottom", legend.key.width=unit(legend.width,"cm"))

  plt.handle3 <- df.meta %>% ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = density1)) + geom_point() + theme_miko(legend = T) +
    xlab("UMI/cell") + ylab("Genes/cell")  +
    labs(title = paste0("r = ", rho3p, "; rho = ", rho3s)) + scale_color_viridis("Density") +
    theme(legend.position="bottom", legend.key.width=unit(legend.width,"cm"))

  plt.QC_scatter <- cowplot::plot_grid(plt.handle1, plt.handle2, plt.handle3, ncol = 3)

  return(plt.QC_scatter)

}



#' Generate heatmap with plotly interface
#'
#' Generate heatmap with plotly interface. Takes output from scMiko::getHeat() function as input.
#'
#' @param mat matrix. Must be same matrix provided to scMiko::getHeat().
#' @param heat.object output from scMiko::getHeat()
#' @param scale.limit Numeric specifying heatmap limits. If unspecified, defined as max(c(abs(max(mat)), abs(min(mat)))).
#' @param legend.label Character.
#' @param x.label x axis label. Character
#' @param y.label y axis label. Character.
#' @param title plot title. Character.
#' @name generateHeat
#' @return plotly handle
#'
generateHeat <- function(mat, heat.object, scale.limit = NULL, legend.label = "", x.label = "", y.label = "", title = ""){

  if (is.null(scale.limit)) scale.limit <- max(c(abs(max(mat)), abs(min(mat))))

  mat.reconstruct <- mat[heat.object[["rowInd"]], heat.object[["colInd"]]]

  # create data frame from reconstructed matrix
  gene.order <- rownames(mat.reconstruct)
  df.recon <- as.data.frame(mat.reconstruct)
  df.recon$genes <- factor(rownames(df.recon), levels = gene.order)
  df.recon.melt <- melt(df.recon)
  df.recon.melt$variable <- factor(df.recon.melt$variable, levels = colnames(mat.reconstruct))

  # define limits
  df.recon.melt$value[df.recon.melt$value > scale.limit] <- scale.limit
  df.recon.melt$value[df.recon.melt$value < -scale.limit] <- -scale.limit

  # heatmap (ggplot2-based)
  # legend name
  fill.label <- legend.label


  p_heat <- ggplot(df.recon.melt, aes(variable, genes)) +
    geom_tile(aes(fill = (value))) +
    scale_fill_gradientn(colours = rev(brewer.pal(9, "RdBu")), limits = c(-scale.limit, scale.limit)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab(x.label) +
    ggtitle(title) + labs(fill = fill.label)

  # helper function for creating dendograms
  ggdend <- function(df) {
    ggplot() +
      geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
      labs(x = "", y = "") +
      ggdendro::theme_dendro() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
  }

  # get dendromgram data from heat object
  dx.cor <- ggdendro::dendro_data(heat.object[["colDendrogram"]])
  dy.cor <- ggdendro::dendro_data(heat.object[["rowDendrogram"]])

  # x/y dendograms
  px <- ggdend(dx.cor$segments)
  py <- ggdend(dy.cor$segments) + coord_flip()

  # hide axis ticks and grid lines
  eaxis <- list(
    showticklabels = FALSE,
    showgrid = FALSE,
    zeroline = FALSE
  )

  p_empty <- plot_ly() %>%
    layout(margin = list(l = 50),
           xaxis = eaxis,
           yaxis = eaxis)

  output <- list(p_heat = p_heat,
                 px = px,
                 py = py,
                 p_empty =p_empty)


  plotly.heat <- subplot(output$px, output$p_empty, output$p_heat, output$py,
                nrows = 2, margin = 0.01,
                widths = c(0.85, 0.15), heights = c(0.15, 0.85)) %>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "myplot",
        width = 600,
        height = 700
      )
    )

  return(plotly.heat)

}


#' Save ggplot as image file
#'
#' Export ggplot to image file
#'
#' @param plt.handle ggplot handle
#' @param save.as save file name. A character.
#' @param format Format of saved image
#' \itemize{
#' \item "png" - Default
#' \item "pdf"
#' \item "jpeg"
#' \item "ps"
#' \item "bmp"
#' \item "wmf"
#' }
#' @param ... pass additional arguments to graphics device functions.
#' @name savePlot
#' @return
#'
savePlot <- function(plt.handle, save.as, format = "png", ...){

  formats <- c("pdf", "png", "jpeg", "ps", "bmp", "wmf")

  if (!(format %in% formats)) stop("Format incorrectly specified")

  filename <- paste(save.as, ".", format, sep = "")

  if (format == "pdf"){
    pdf(filename, ...)
  } else if (format == "png"){
    png(filename, ...)
  } else if (format == "jpeg"){
    jpeg(filename, ...)
  } else if (format == "ps"){
    postscript(filename, ...)
  } else if (format == "bmp"){
    bmp(filename, ...)
  } else if (format == "wmf"){
    win.metafile(filename, ...)
  }
  print(plt.handle)
  print.suppres <- capture.output(dev.off())

}

#' Plot relationship showing percentage of cells expressing atleast percentage of genes
#'
#' Curve illustrating percentage of cells expressing atleast percentage of genes
#'
#' @param so Seurat object
#' @param which.genes Specify geneset. If null, all gene are included (unless only.variable is set to True)
#' @param only.variable Only used variable genes
#' @param which.data specify which data slot to use (scale or data; data is default)
#' @name geneRepCurve
#' @return ggplot handle
#'
geneRepCurve <- function(so, which.genes = NULL, only.variable = F, which.data = "data"){


  # binarize expression
  exp.mat <- getExpressionMatrix(so, only.variable = only.variable, which.data = which.data)

  # get subset of genes if specified
  if (!is.null(which.genes)){
    exp.mat <- exp.mat[rownames(exp.mat) %in% which.genes, ]
  }

  exp.mat.bin <- exp.mat
  exp.mat.bin[exp.mat.bin > 0] <-  1

  per.vec <- seq(1:100) / 100
  exp.mat.colsum <- apply(exp.mat.bin,2, function(x) sum(x)/length(x))

  per.exp <- c()
  for (i in 1:length(per.vec)){
    per.exp[i] <- sum(exp.mat.colsum > per.vec[i]) / length(exp.mat.colsum)
  }

  df.per <- data.frame(x = per.vec, y = per.exp)

  plt.per <- df.per %>%
    ggplot(aes(x=x, y = y)) +
    geom_point() +
    theme_classic() +
    xlab("% Genes Expressed Per Cell (minimum)") +
    ylab("% of Cell") +
    ggtitle("Genes Detected")

  return(plt.per)

}



#' Plot network properties for different soft thresholds
#'
#' Plot network properties for different soft thresholds
#'
#' @param sft list of soft threshold picks (output from scMiko::getSoftThreshold)
#' @param r2.threshold Numeric specifying threshold at which to draw red horizontal curve. Default is 0.85.
#' @name networkProperties.Plot
#' @return plot
#' @examples
#'
#' sft <- getSoftThreshold(s.mat, networkType = "unsigned")
#' networkProperties.Plot(sft)
#'
#'
networkProperties.Plot <- function(sft, r2.threshold = 0.85){
  # Plot the results
  sizeGrWindow(9, 5)
  par(mfrow = c(1,2));
  cex1 = 0.9;
  powers = c(c(1:10), seq(from = 12, to=20, by=2));

  # Scale-free topology fit index as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");

  # Red line corresponds to using an R^2 cut-off
  abline(h=r2.threshold,col="red")

  # Mean connectivity as a function of the soft-thresholding power
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
}



#' Violin plot of single cell gene expression
#'
#' Draws violin plot of single cell data, stratified by specified group. Hierarchically-clustered dendrogram is overlaid to illustrate relationships between groups.
#'
#' @param so Seurat object
#' @param which.gene query gene. Must be available in Seurat object.
#' @param e.mat expression matrix (gene x cell). If not provided, computed from seurat object (more computationally heavy).
#' @param f.mat expression fraction matrix (gene x cell). If not provided, computed from seurat object (more computationally heavy).
#' @param which.group Group to cluster gene expression by. Must be column header in so meta data. Default is "seurat_clusters".
#' @param which.data Specifies which seurat data used to compute e.mat and f.mat. Must be either "data" (default) or "scale". Used only if e.mat and f.mat are not provided.
#' @param which.assay Specifies which seurat assay to used to compute e.mat. Used only if e.mat is not provided. Default is DefaultAssay(so)
#' @param x.label x axis title (a character). Default is which.group.
#' @param x.label.angle rotation angle for x axis title (a numeric). Default is NULL.
#' @name expression.Plot
#' @return ggplot object
#' @examples
#'
#' # get expression plot of GPNMB stratified by cluster group.
#' plt.sgExp.clust <- expression.Plot.dev(so.query, which.gene = "GPNMB", which.group = "seurat_clusters", which.data = "data")
#' print(plt.sgExp.clust)
#'
#'
expression.Plot <- function(so, which.gene, e.mat = NULL, f.mat = NULL,
                            which.group = "seurat_clusters", which.data = "data", which.assay = DefaultAssay(so),
                            x.label = NULL, x.label.angle = NULL){

  # get expression data
  if (is.null(e.mat)){
    em <- getExpressionMatrix(
      so,
      only.variable = F,
      which.assay = which.assay,
      which.data = which.data,
      use.additional.genes = NA
    )
  } else {
    em <- e.mat
  }

  # get fraction of expressing cells
  if (is.null(f.mat)){
    em.frac <- avgGroupExpression(
      so,
      which.data = "data",
      which.center = "fraction",
      which.group = which.group
    )
  } else {
    em.frac <- f.mat
  }

  if ("genes" %in% colnames(em)){
    rownames(em) <- em$genes
    em <- em %>% dplyr::select(-c("genes"))
  }
  if ("genes" %in% colnames(em.frac)){
    rownames(em.frac) <- em.frac$genes
    em.frac <- em.frac %>% dplyr::select(-c("genes"))
  }

  # get meta data
  so.meta <- so@meta.data

  # ensure grouping variable exists
  if (!(which.group) %in% colnames(so.meta)) stop("Grouping variable does not exist.")

  # merge datasets
  if (sum(rownames(em) %in% (which.gene)) != 1) stop("More than one gene matched query. Cannot continue.")
  df.meta <- data.frame(id = rownames(so.meta), group = so.meta[,which.group])
  em.mark <- as.data.frame(em[(rownames(em)) %in% (which.gene), ])
  em.df <- data.frame(id = rownames(em.mark), query = em.mark[,1])
  em.merge <- merge(em.df, df.meta, by = "id")

  # get summary statistics (for clustering)
  suppressMessages({em.sum <- em.merge %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(x.mean = (mean((query), na.rm = T)),
                     x.sd = sd(query, na.rm = T))})
  em.frac.mark <- as.data.frame(em.frac[(rownames(em.frac)) %in% (which.gene), ])
  if (rownames(em.frac.mark) %in% which.gene) em.frac.mark <-as.data.frame(t(em.frac.mark))
  colnames(em.frac.mark) <- which.gene
  em.frac.mark.df <- data.frame(group = rownames(em.frac.mark), frac = em.frac.mark[,1])

  if (which.group == "seurat_clusters") em.frac.mark.df$group <- gsub("c", "", em.frac.mark.df$group)
  em.sum <- merge(em.sum, em.frac.mark.df, by = "group")

  # hierarchial clustering
  row.names.df <- em.sum$group
  em.sum <- em.sum %>% dplyr::select(-c("group"))
  clust.var <- as.matrix(em.sum)
  rownames(clust.var) <- row.names.df
  d <- dist(as.matrix(clust.var))   # find distance matrix

  clust.success <- F
  try({
    hc <- hclust(d)                # apply hirarchical clustering
    clust.success <- T
  }, silent = T)


  if (is.null(x.label)) x.label <- which.group

  if (clust.success){
    # helper function for creating dendograms
    ggdend.v2 <- function(df) {
      ggplot() +
        geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
        ggdendro::theme_dendro()
    }

    # get dendromgram data from heat object
    d.query.clust <- ggdendro::dendro_data(hc)

    # ggplot dendograms
    p.query.clust <- ggdend.v2(d.query.clust$segments) +
      ggtitle(paste0(which.gene, " Expression")) +
      theme(axis.title.x=element_blank(),
            axis.text=element_blank(),
            axis.line=element_blank(),
            axis.ticks=element_blank(),
            legend.position = "none")

    # reorder data according to clusters
    em.merge$group <- factor(as.character(em.merge$group), levels = (d.query.clust[["labels"]][["label"]]))

    max.query <- max(em.merge$query, na.rm = T)
    em.merge$query.norm <- rescaleValues(values = em.merge$query, new.min = 0, new.max = 1)

    suppressMessages({em.merge.sum <- em.merge %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(ef = mean(query > 0),
                ev = mean(query.norm, na.rm = T))})


    plt.em <- ggplot() +
      geom_bar(data = em.merge.sum, aes(x = group, y = ef, fill = group), stat = "identity", alpha = 0.5) +
      coord_cartesian(ylim = c(0, 1)) +
      geom_violin(data = em.merge, aes(x = group, y = query.norm, fill = group)) +
      geom_point(data = em.merge.sum, aes(x = group, y = (ev), fill = group)) +
      theme_miko() +
      xlab(x.label)  +
      scale_y_continuous(sec.axis = sec_axis(~., name = "Normalized Expression (violin)"), name = "Expressing Fraction (bar)") +
      theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
            legend.position = "none")

    if (!is.null(x.label.angle) && is.numeric(x.label.angle)){
      plt.em <- plt.em + theme(axis.text.x = element_text(angle = x.label.angle, hjust = 1))
    }

    # combine plots
    plt.sgExp <-  cowplot::plot_grid(p.query.clust, plt.em, ncol = 1, align = "v", rel_heights = c(1,3))

  } else {

    max.query <- max(em.merge$query, na.rm = T)
    em.merge.sum <- em.merge %>%
      dplyr::group_by(group) %>%
      dplyr::summarize(ef = mean(query > 0),
                ev = mean(query, na.rm = T))
    plt.sgExp <- ggplot() +
      geom_bar(data = em.merge.sum, aes(x = group, y = ef, fill = group), stat = "identity", alpha = 0.5) +
      coord_cartesian(ylim = c(0, 1)) +
      geom_violin(data = em.merge, aes(x = group, y = query/max.query, fill = group)) +
      geom_point(data = em.merge.sum, aes(x = group, y = (ev)/max.query, fill = group)) +
      theme_miko() +
      xlab(x.label)  +
      scale_y_continuous(sec.axis = sec_axis(~., name = "Expression (violin)"), name = "Expressing Fraction (bar)") +
      theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
            legend.position = "none")

    if (!is.null(x.label.angle) && is.numeric(x.label.angle)){
      plt.em <- plt.em + theme(axis.text.x = element_text(angle = x.label.angle, hjust = 1))
    }
  }

  return(plt.sgExp)

}



#' Split violin plot using ggplot2
#'
#' Draws split violin plot using ggplot2. Credit to jan-glx (https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2)
#'
#' @param mapping ggplot mapping
#' @param data data.frame
#' @param stat ggplot identity
#' @param position ggplot position
#' @param ... additional parameters to pass to layer()
#' @param draw_quantiles draw quantiles. Default is NULL
#' @param trim Trim data. Default is True.
#' @param scale Default is area.
#' @param na.rm Default is na.rm
#' @param show.legend Default is NA
#' @param inherent.aes Default is TRUE
#' @name geom_split_violin
#' @return ggplot object
#' @examples
#'
#' # split violin plot
#' ggplot(my_data, aes(x, y, fill = m)) + geom_split_violin()
#'
#'
geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ...,
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE,
                              show.legend = NA, inherit.aes = TRUE) {


  GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin,
                             draw_group = function(self, data, ..., draw_quantiles = NULL) {
                               data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                               grp <- data[1, "group"]
                               newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                               newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                               newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

                               if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                                 stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                           1))
                                 quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                                 aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                                 aesthetics$alpha <- rep(1, nrow(quantiles))
                                 both <- cbind(quantiles, aesthetics)
                                 quantile_grob <- GeomPath$draw_panel(both, ...)
                                 ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                               }
                               else {
                                 ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                               }
                             })



  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin,
        position = position, show.legend = show.legend, inherit.aes = inherit.aes,
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}



#' Upset plot
#'
#' Generate upset plot, comparing intersection between different (gene)sets. Uses ComplexHeatmap Package.
#'
#' @param gene.sets named list of genesets, where names specify name of gene set, and entries are character vectors specifying genes belongs to the respective set.
#' @param row.title Row title.
#' @param column.title Column title.
#' @name upset.Plot
#' @return plot handle
#' @examples
#'
#' Generate upset plot for list of genesets
#' plt.upset <- upset.Plot(gene.sets, row.title = "Genesets")
#'
#' # print plot
#' print(plt.upset)
#'
upset.Plot <- function(gene.sets, row.title = "", column.title = ""){

  if (!require(ComplexHeatmap)) stop("ComplexHeatmap package not found")
  m = ComplexHeatmap::make_comb_mat(gene.sets)

  plt.upset <- ComplexHeatmap::UpSet(m, top_annotation = HeatmapAnnotation(
    degree = as.character(comb_degree(m)),
    "Intersection\nsize" = anno_barplot(comb_size(m),
                                        border = FALSE,
                                        gp = gpar(fill = "black"),
                                        height = unit(2, "cm")
    ),
    annotation_name_side = "left",
    annotation_name_rot = 0),
    row_title = row.title,
    column_title = column.title)

  return(plt.upset)
}


#' scMiko Theme
#'
#' Custom themes applied to ggplot2-based plots.
#'
#' @param style Choise of 'theme_classic' or 'theme_bw' called using 'classic' or 'bw', respectively. Default is 'bw'
#' @param legend Logical to include legend. Default is F.
#' @param grid Logical to include grid. Default is F.
#' @param bold.title Logical to bold title. Default is T.
#' @param center.title Logical to center justify title. Default is F.
#' @param x.axis.rotation rotation angle (e.g., 45)
#' @param fill.palette palette from ggthemes (e.g., "ptol"). Default is NA.
#' @param color.palette palette from ggthemes (e.g., "ptol"). Default is NA.
#' @name theme_miko
#' @return ggplot2 theme object
#' @examples
#'
#'
theme_miko <- function(style = "bw", legend = F, grid = F, bold.title = T, center.title = F, x.axis.rotation = 0, fill.palette = NA, color.palette = NA){

  if (style == "bw"){
    tm <- theme_bw()
  } else if (style == "classic"){
    tm <- theme_classic()
  }

  if (!legend) tm <- tm + theme(legend.position = "none")

  if (!grid) tm <- tm + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

  if (bold.title) tm <- tm + theme( plot.title = element_text(face = "bold"))

  if (center.title){
    tm <- tm +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5))
  }

  if (x.axis.rotation != 0){
    try({tm <- tm + theme(axis.text.x = element_text(angle = 45, hjust = 1))}, silent = T)
  }

  if (!is.na(fill.palette)){
    require(ggthemes)
    try({tm <- tm + do.call(paste0("scale_fill_", "ptol"), args = list())}, silent =  T)
  }

  if (!is.na(color.palette)){
    require(ggthemes)
    try({tm <- tm + do.call(paste0("scale_colour_", "ptol"), args = list())}, silent =  T)
  }

  return(tm)

}


#' UMAP projection of pseudotimes
#'
#' Visualize pseudotimes using UMAP coordinates.
#'
#' @param x x coordinates (e.g., UMAP 1)
#' @param y y coordinates (e.g., UMAP 2)
#' @param pseudotime pseudotimes. Numerical vector, same length as x and y.
#' @param pt.size Point size. Numeric [0,inf]
#' @param pt.alpha Point alpha. Numeric [0,1]
#' @param x.lab x axis label. Default is 'UMAP 1'.
#' @param y.lab y axis label. Default is 'UMAP 2'.
#' @name pseudotime.UMAP
#' @return ggplot2 handle
#' @examples
#'
#'
pseudotime.UMAP <- function(x, y, pseudotime, pt.size = 1, pt.alpha = 1, x.lab = "UMAP 1", y.lab = "UMAP 2"){

  df <- data.frame(x,y,pseudotime)

  plt.pt <- df %>%
    ggplot(aes(x, y, color = pseudotime)) +
    geom_point(size = pt.size, alpha = pt.alpha) +
    theme_classic() +
    xlab(x.lab) +
    ylab(y.lab) +
    viridis::scale_color_viridis()

  return(plt.pt)

}


#' Apply discrete color palette to ggplot object
#'
#' Apply discrete color palette to ggplot object
#'
#' @param gg ggplot handle
#' @param fc character specifying whether to apply color palette to 'fill', 'color', or 'both'
#' @param n.groups numeric specifying number of groups
#' @name discretePalette
#' @return ggplot2 handle
#' @examples
#'
#' gg.plot <- discretePalette(gg.plot, fc = "color", n.groups = 4)
#'
discretePalette <- function(gg, fc, n.groups){

  if (n.groups <= 10){

    if (fc == "color"){
      gg <- gg + ggsci::scale_color_npg()
    } else if (fc == "fill"){
      gg <- gg + ggsci::scale_fill_npg()
    } else if (fc == "both"){
      gg <- gg + ggsci::scale_color_npg() + ggsci::scale_fill_npg()
    }
  }
  return(gg)
}


#' Automatically determine optimal point size for geom_point()
#'
#' Automatically determine optimal point size for geom_point()
#'
#' @param n.points Number of data points
#' @param scale.factor scaling constant used to determine optimal point size.
#' @name autoPointSize
#' @return optimal point size
#' @examples
#'
#'  gg.plot <-  df.umap %>% dplyr::arrange(get(module.names[i])) %>%
#'       ggplot(aes(x = x, y = y, color = get(module.names[i]))) +
#'       geom_point(size = autoPointSize(nrow(df.umap)))
#'
autoPointSize <- function(n.points, scale.factor = 1583){
  min(scale.factor/n.points, 1)
}
