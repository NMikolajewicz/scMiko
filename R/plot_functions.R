
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
#' @param show.axis show graph axis (omit theme_void from existing ggplot). Default is T.
#' @param ... additional arguments passed to Seurat::DimPlot().
#' @name cluster.UMAP
#' @return ggplot handle
#'
cluster.UMAP <- function(so, group.by = "seurat_clusters", x.label = "UMAP 1", y.label = "UMAP 2", plot.name = "UMAP", include.labels = T, reduction = "umap",pt.size = autoPointSize(ncol(so)), show.axis = T, ...){

  if (group.by %in% colnames(so@meta.data)){
    if (group.by == "seurat_clusters") so@meta.data[["seurat_clusters"]] <- orderedFactor(so@meta.data[["seurat_clusters"]])
  } else {
    group.by = NULL
  }

  plt.handle <- DimPlot(so, group.by = group.by, label = include.labels,reduction = reduction, pt.size = pt.size, ...)  +
    ggtitle(label = plot.name) +
    xlab(x.label) + ylab(y.label) +
    theme_miko(legend = T)

  if (!show.axis){
    plt.handle <- plt.handle +  theme_void()
  }

  return(plt.handle)

}




#' Cell-level gene expression projected on UMAP
#'
#' UMAP plot with cell-level gene expression for queried gene. Uses Seurat::FeaturePlot().
#'
#' @param object Seurat Object
#' @param query.gene Character specifying gene name to plot.
#' @param x.label Character specifying X axis label.
#' @param y.label Character specifying Y axis label.
#' @param adjust.pt.size Adjust point size for plotting. Logical.
#' @param order.cells Plot cells in order of expression. Logical.
#' @param plot.name Character specifying plot title.
#' @param reduction Character specifying reduction name. Default is "umap".
#' @param ... additional parameters passed to Seurat::FeaturePlot().
#' @name scExpression.UMAP
#' @return ggplot handle
#'
scExpression.UMAP <- function(object, query.gene, x.label = "UMAP 1", y.label = "", adjust.pt.size = autoPointSize(ncol(object)), order.cells = T, plot.name = NULL, reduction = "umap",  ...){

  if (is.null(plot.name)) plot.name <- query.gene
  # cols =rev(brewer.pal(11,"RdYlBu")),
  plt.handle <- FeaturePlot(object = object, features = query.gene,
                            reduction = reduction,pt.size = adjust.pt.size, sort.cell = order.cells, ...) +
    xlab(x.label) +
    ylab(y.label) +
    ggtitle(plot.name)

  return(plt.handle)

}


#' Plot variable genes in Seurat Object
#'
#' Using Seurat's VariableFeaturePlot, variable genes are identified and plotted. Top N variable genes are labeled.
#'
#' @param so Seurat Object
#' @param gNames Named gene list; entries are Symbols, names are Ensemble.
#' @param set_name Character specfiying name of dataset. Optional.
#' @param top.n.genes Numeric. Top n genes to label on plot.
#' @param pt.colors colors to specify non-variable/variable status.
#' @param ... additional arguments passed to Seurat::VariableFeaturePlot(...)
#' @name variableGenes.Plot
#' @return ggplot handle
#'
variableGenes.Plot <- function(so, gNames, set_name = NULL, top.n.genes = 10, pt.colors = c("black", "tomato"), ...){

  top10 <- head(VariableFeatures(so), assay = "gene_name", top.n.genes)

  # Plot variable features
  plt.handle <- VariableFeaturePlot(so, ...)

  # plot title
  if (is.null(set_name) || set_name == ""){
    plt.title <- "Variable Genes"
  } else {
    plt.title <- paste(set_name, ": Variable Genes")
  }

  var_labs <- c()
  try({
    gene.rep <- checkGeneRep(reference.genes = gNames, query.genes = top10 )
    if (gene.rep == "ensembl"){
      var_labs <- as.vector(gNames[top10])
    } else if (gene.rep == "symbol"){
      var_labs <- as.vector(top10)
    }
  }, silent = T)

  df.dat <- plt.handle[["data"]]
  gene.rep <- checkGeneRep(reference.genes = gNames, query.genes = rownames(df.dat) )
  if (gene.rep == "ensembl"){
    df.dat$gene <- as.vector(gNames[rownames(df.dat)])
  } else if (gene.rep == "symbol"){
    df.dat$gene  <-  rownames(df.dat)
  }
  df.dat.top <- df.dat %>% dplyr::filter(gene %in% top10)

  gene.rep <- checkGeneRep(reference.genes = gNames, query.genes =  VariableFeatures(so))
  if (gene.rep == "ensembl"){
    var.gene <- as.vector(gNames[VariableFeatures(so)])
  } else if (gene.rep == "symbol"){
    var.gene  <- VariableFeatures(so)
  }

  df.dat <- df.dat[complete.cases(df.dat), ]
  df.dat$is.var <- paste0("Non-Variable: ", nrow(df.dat) - length(var.gene))
  df.dat$is.var[df.dat$gene %in% var.gene] <- paste0("Variable: ", length(var.gene))

  plt.handle <- df.dat %>%
    ggplot(aes(x = gmean, y = residual_variance, color = is.var)) +
    geom_point(size = 1) +
    ggrepel::geom_text_repel(data= df.dat.top, aes(x = gmean, y = residual_variance, label = gene), inherit.aes = F, max.overlaps = Inf) +
    theme_miko(legend = T) +
    scale_color_manual(values = pt.colors) +
    scale_x_log10() +
    labs(x = "Geometric Mean of Expression", y = "Residual Variance", color = NULL)


  return(plt.handle)
}



#' QC violin plots
#'
#' Number of genes/cell, UMI/cell and mitochondiral content/cell are visualized with violin plots. Two ggplot handles are generated. First contains QC metrics pooled across all cells, while second stratifies dataset by grouping variable.
#'
#' @param so Seurat Object
#' @param features meta data features to plot. Default is c("nFeature_RNA", "nCount_RNA", "percent.mt").
#' @param group.by meta data field to group plots by.
#' @param plt.log.flag Logical specifying whether data are plotted on log scale. Default is True.
#' @param ... additional arguments passed to Seurat::VlnPlot(...)
#' @name QC.violinPlot
#' @return list of ggplot handles
#'
QC.violinPlot <- function(so, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = NULL, plt.log.flag = T, ...){

  stopifnot("'so' must be seurat object"  = "Seurat" %in% class(so))
  if (!is.null(group.by)){
    if (!(group.by %in% colnames(so@meta.data))) group.by <- NULL
  }

  features2plot <- features
  features2plot <- features2plot[features2plot %in% colnames(so@meta.data)]

  if (plt.log.flag){
    # logarithmic scale
    plt1 <- VlnPlot(so, features = features2plot, ncol = 3, log = TRUE, ...)
    plt2 <- VlnPlot(so, features = features2plot, group.by = group.by, ncol = 3, log = TRUE, ...)
  } else{
    # raw scale
    plt1 <- VlnPlot(so, features = features2plot, ncol = 3, log = FALSE, ...)
    plt2 <- VlnPlot(so, features = features2plot, group.by = group.by, ncol = 3, log = FALSE, ...)
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

  stopifnot("'so' must be seurat object"  = "Seurat" %in% class(so))
  stopifnot("'nFeature_RNA', 'nCount_RNA', 'percent.mt' must all be provided in seurat object meta data"  = all(c("nFeature_RNA", "nCount_RNA", "percent.mt") %in%  colnames(so@meta.data)))
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

  plt.handle2 <- df.meta %>% ggplot(aes(x = nFeature_RNA, y = percent.mt, color = density2)) + geom_point() + theme_miko(legend = T) +
    xlab("Genes/cell") + ylab("Mitochondrial Content (%)")  +
    labs(title = paste0("r = ", rho2p, "; rho = ", rho2s)) + scale_color_viridis("Density") +
    theme(legend.position="bottom", legend.key.width=unit(legend.width,"cm"))

  plt.handle3 <- df.meta %>% ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = density3)) + geom_point() + theme_miko(legend = T) +
    xlab("UMI/cell") + ylab("Genes/cell")  +
    labs(title = paste0("r = ", rho3p, "; rho = ", rho3s)) + scale_color_viridis("Density") +
    theme(legend.position="bottom", legend.key.width=unit(legend.width,"cm"))

  plt.QC_scatter <- cowplot::plot_grid(plt.handle1, plt.handle2, plt.handle3, ncol = 3)

  return(plt.QC_scatter)

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
#' @param show.violin Show expression as violin plot. Default is T.
#' @param show.full.axis Show full normalized axis [0,1], rather than adapting to data. Default is T.
#' @param font.size Size of font.
#' @param verbose print progress. Default is T.
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
                            x.label = NULL, x.label.angle = NULL, show.violin = T, show.full.axis = T, font.size = NULL,  verbose = T){

  if (!("Seurat" %in% class(so))) stop("'so' is not a Seurat object")
  which.gene <- which.gene[which.gene %in% rownames(so)]
  if (length(which.gene) == 0) stop("No features provided in 'which.gene' are available")

  miko_message("Preparing Expression Data...", verbose = verbose)


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
      which.features = which.gene,
      which.data = "data",
      which.center = "fraction",
      which.group = which.group,
      verbose = F
    )
  } else {
    em.frac <- f.mat
  }

  if ("genes" %in% colnames(em)){
    rownames(em) <- make.unique(em$genes)
    em <- em %>% dplyr::select(-c("genes"))
  }
  if ("genes" %in% colnames(em.frac)){
    rownames(em.frac) <- make.unique(em.frac$genes)
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
  miko_message("Clustering groups...", verbose = verbose)
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

  miko_message("Generating plot...", verbose = verbose)
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
      # coord_cartesian(ylim = c(0, 1)) +
      # geom_violin(data = em.merge, aes(x = group, y = query.norm, fill = group)) +
      geom_point(data = em.merge.sum, aes(x = group, y = (ev), fill = group)) +
      theme_miko() +
      xlab(x.label)  +
      scale_y_continuous(sec.axis = sec_axis(~., name = "Normalized Expression (dot)"), name = "Expressing Fraction (bar)")

      # theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
      #       legend.position = "none")

    if (show.violin){
      plt.em <- plt.em +
        geom_violin(data = em.merge, aes(x = group, y = query.norm, fill = group)) +
        geom_point(data = em.merge.sum, aes(x = group, y = (ev), fill = group))  +
        scale_y_continuous(sec.axis = sec_axis(~., name = "Normalized Expression (violin)"), name = "Expressing Fraction (bar)")
    }

    if (show.full.axis){
      plt.em <- plt.em + coord_cartesian(ylim = c(0, 1)) +
        theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
              legend.position = "none")
    }


    if (!is.null(x.label.angle) && is.numeric(x.label.angle)){
      plt.em <- plt.em + theme(axis.text.x = element_text(angle = x.label.angle, hjust = 1))
    }

    if (!is.null(font.size)){
      plt.em <- plt.em + theme(text = element_text(size=font.size))
      p.query.clust <- p.query.clust + theme(text = element_text(size=font.size))
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
      # geom_violin(data = em.merge, aes(x = group, y = query/max.query, fill = group)) +
      geom_point(data = em.merge.sum, aes(x = group, y = (ev)/max.query, fill = group)) +
      theme_miko() +
      xlab(x.label)  +
      scale_y_continuous(sec.axis = sec_axis(~., name = "Normalized Expression (dot)"), name = "Expressing Fraction (bar)")
      # theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
      #       legend.position = "none")

    if (show.violin){
      plt.sgExp <- plt.sgExp +
        geom_violin(data = em.merge, aes(x = group, y = query/max.query, fill = group)) +
        geom_point(data = em.merge.sum, aes(x = group, y = (ev)/max.query, fill = group))    +
        scale_y_continuous(sec.axis = sec_axis(~., name = "Expression (violin)"), name = "Expressing Fraction (bar)")
    }

    if (show.full.axis){
      plt.sgExp <- plt.sgExp +
        coord_cartesian(ylim = c(0, 1)) +
        theme(plot.margin = unit(c(0, 1, 0, 1), "cm"),
                               legend.position = "none")
    }

    if (!is.null(x.label.angle) && is.numeric(x.label.angle)){
      plt.sgExp <- plt.sgExp + theme(axis.text.x = element_text(angle = x.label.angle, hjust = 1))
    }

    if (!is.null(font.size)){
      plt.sgExp <- plt.sgExp + theme(text = element_text(size=font.size))
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
#' @param color.luminescence color luminescence, passed to `scale_color_hue()`.
#' @param fill.luminescence fill luminescence, passed to `scale_fill_hue()`.
#' @name theme_miko
#' @return ggplot2 theme object
#' @examples
#'
#'
theme_miko <- function(style = "bw", legend = F, grid = F, bold.title = T, center.title = F, x.axis.rotation = 0, fill.palette = NA, color.palette = NA, color.luminescence = NA, fill.luminescence = NA){

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
    try({tm <- tm + theme(axis.text.x = element_text(angle = x.axis.rotation, hjust = 1))}, silent = T)
  }



  # if (!is.na(legend.color.size)){
  #   tm <- tm + guides(color = guide_legend(override.aes = list(size=legend.color.size)))
  # }
  #
  # if (!is.na(legend.fill.size)){
  #   tm <- tm + guides(fill = guide_legend(override.aes = list(size=legend.fill.size)))
  # }
  #
  # if (do.minimal){
  #   tm <- tm +  theme(
  #     panel.border = element_blank(),
  #     axis.text = element_blank(),
  #     panel.grid = element_blank(),
  #     axis.ticks = element_blank()
  #   ) + xlab("") + ylab("")
  # }


  tm.list <- list()
  list.ind <- 1
  if (!is.na(fill.palette)){
    if (fill.palette == "flatly"){
      flatly.colors <- c('#18BC9C','#2C3E50','#F39C12','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_fill_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else if (fill.palette == "flatly2"){
      flatly.colors <- c('#2C3E50','#18BC9C','#F39C12','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_fill_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else if (fill.palette == "flatly3"){
      flatly.colors <- c('#2C3E50','#F39C12','#18BC9C','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_fill_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else if (fill.palette == "flatly4"){
      flatly.colors <- c('#F39C12','#2C3E50','#18BC9C','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_fill_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else {
      require(ggthemes)
      tm.list[[list.ind]] <- do.call(paste0("scale_fill_", fill.palette), args = list())
      list.ind <- list.ind + 1
    }
  }

  if (!is.na(color.palette)){
    if (color.palette == "flatly"){
      flatly.colors <- c('#18BC9C','#2C3E50','#F39C12','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_color_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else if (color.palette == "flatly2"){
      flatly.colors <- c('#2C3E50','#18BC9C','#F39C12','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_color_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else if (color.palette == "flatly3"){
      flatly.colors <- c('#2C3E50','#F39C12','#18BC9C','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_color_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else if (color.palette == "flatly4"){
      flatly.colors <- c('#F39C12','#2C3E50','#18BC9C','#E74C3C','#3498DB','#18BC9C','#2C3E50','#F39C12')
      tm.list[[list.ind]] <- do.call(paste0("scale_color_manual"), args = list(values = flatly.colors))
      list.ind <- list.ind + 1
    } else {
    require(ggthemes)
    tm.list[[list.ind]] <- do.call(paste0("scale_colour_", color.palette), args = list())
    list.ind <- list.ind + 1
    }
  }


  if (!is.na(color.luminescence)){
    tm.list[[list.ind]] <- do.call("scale_color_hue", args = list(l=color.luminescence))
    list.ind <- list.ind + 1

    # tm <- tm + scale_color_hue(l=color.luminescence)
  }

  if (!is.na(fill.luminescence)){
    tm.list[[list.ind]] <- do.call("scale_fill_hue", args = list(l=fill.luminescence))
    list.ind <- list.ind + 1

    # tm <- tm + scale_fill_hue(l=fill.luminescence)
  }

  tm.list[[list.ind]] <- tm

  return(tm.list)

}




#' Automatically determine optimal point size for geom_point()
#'
#' Automatically determine optimal point size for geom_point()
#'
#' @param n.points Number of data points
#' @param scale.factor scaling constant used to determine optimal point size.
#' @param max.size max point size
#' @param min.size min point size
#' @name autoPointSize
#' @return optimal point size
#' @examples
#'
#'  gg.plot <-  df.umap %>% dplyr::arrange(get(module.names[i])) %>%
#'       ggplot(aes(x = x, y = y, color = get(module.names[i]))) +
#'       geom_point(size = autoPointSize(nrow(df.umap)))
#'
autoPointSize <- function(n.points, scale.factor = 10000, max.size = 2, min.size = 0.01){

  # df.scaling <- data.frame(x = 10^(seq(0, 6, by = 0.1)))
  # df.scaling$y = scale.factor/df.scaling$x
  # df.scaling$y[df.scaling$y > max.size] <- max.size
  # df.scaling$y[df.scaling$y < min.size] <- min.size
  #
  # df.scaling %>%
  #   ggplot(aes(x = log10(x), y = y)) +
  #   geom_point()

  p.size <- scale.factor/n.points
  if (p.size > max.size) p.size <- max.size
  if (p.size < min.size) p.size <- min.size
  return(p.size)
}


#' Visualize feature activity/expression gradient overlaid on UMAP
#'
#' Visualize feature activity/expression gradient overlaid on UMAP
#'
#' @param object Seurat object.
#' @param feature scaling constant used to determine optimal point size.
#' @param umap.key name of UMAP reduction. Default is "umap". Requires that UMAP embedding was calculated using Seurat::RunUMAP() prior to running scMiko::featureGradient().
#' @param min.quantile.cutoff upper quantile threshold at which values are winsorized. Default is 0.975
#' @param max.quantile.cutoff lower quantile threshold at which values are winsorized. Default is 0.025
#' @param arrow.color gradient arrow color. Default is "black".
#' @param arrow.size gradient arrow size Default is 1.
#' @name featureGradient
#' @return ggplot handle
#' @examples
#'
#'  gg.plot <-  featureGradient(object = so, feature = "S.Score")
#'
featureGradient <- function(object, feature,  umap.key = "umap", min.quantile.cutoff = 0.025, max.quantile.cutoff = 0.975, arrow.color = "black", arrow.size = 1){

  df.umap <- getUMAP(object, umap.key = umap.key)[["df.umap"]]

  if (sum(rownames(object) %in% feature) == 1){
    df.umap$z <- object@assays[[DefaultAssay(object)]]@data[rownames(object) %in% feature, ]
  } else if (sum(colnames(object@meta.data) %in% feature) == 1){
    df.umap$z <-object@meta.data[ ,feature]
  } else {
    stop(feature, " was not found")
  }

  # df.umap <- df.umap[complete.cases(df.umap), ]
  df.umap <- df.umap[complete.cases(df.umap[ ,c("x", "y", "z")]), ]

  lm.fit <- lm(z ~ x + y, data = df.umap)
  df.umap$z.fit <- lm.fit[["fitted.values"]]
  min.x <- min(df.umap$x)
  max.x <- max(df.umap$x)
  min.y <- min(df.umap$y)
  max.y <- max(df.umap$y)
  x.val <- c(min.x, max.x)
  xy.slope <- ((lm.fit[["coefficients"]][["x"]]*1))/(-lm.fit[["coefficients"]][["y"]]) / 1
  x.val <- seq(min.x, max.x, abs((max.x-min.x)/1000))
  df.grad <- data.frame(x = x.val, y = -(1/xy.slope) *x.val)
  df.grad <- df.grad[ (df.grad$y < max.y) & (df.grad$y > min.y ),]
  df.grad <- df.grad[ c(which.min(df.grad$x), which.max(df.grad$x)), ]

  df.grad$z <- (lm.fit[["coefficients"]][["x"]]*df.grad$x) +  (lm.fit[["coefficients"]][["y"]]*df.grad$y) + (lm.fit[["coefficients"]][["(Intercept)"]])
  df.grad <- df.grad %>% dplyr::arrange(z)
  q99 <- quantile(df.umap$z, max.quantile.cutoff)
  q01 <- quantile(df.umap$z, min.quantile.cutoff)
  df.umap$z[df.umap$z>q99] <- q99
  df.umap$z[df.umap$z<q01] <- q01

  lm.sum <- summary(lm.fit)
  r2.val <- paste0("R2 = ", signif(lm.sum[["r.squared"]], 3))

  df.umap %>% ggplot(aes(x =x , y = y, color=z)) +
    geom_point(size = autoPointSize(nrow(df.umap))) +
    geom_segment(x = df.grad$x[1],y = df.grad$y[1], xend = df.grad$x[2],
                 yend = df.grad$y[2], inherit.aes = F, arrow = arrow(), color = arrow.color, size = arrow.size) +
    scale_color_distiller(palette = "RdYlBu") +
    theme_miko(legend = T) +
    labs(x = "UMAP 1", y = "UMAP 2", title = paste0(feature, " gradient"), subtitle = r2.val, color = "")


}

#' Visualize gene expression on UMAP
#'
#' Visualize gene expression on UMAP. Wrapper for Nebulosa::plot_density() and scMiko::scExpression.UMAP() functions.
#'
#' @param object Seurat object.
#' @param feature feature name.
#' @param plot.subtitle Plot title.
#' @param do.neb Logical to use Nebulosa::plot_density instead of scMiko::scExpression.UMAP. Default is F.
#' @param title.size Size of plot title. Default is 10.
#' @param slot which slot to pull from? Default is "data".
#' @param assay assay to retrieve data from.
#' @param reduction reduction name. Default is "umap".
#' @param size point size for plots.
#' @param scale.color color of scale. Default is "tomato".
#' @param ... additional parameters passed to scExpression.UMAP (do.neb = F) or plot_density (do.neb = T)
#' @name exprUMAP
#' @return ggplot handle
#' @examples
#'
#'  gg.plot <-  exprUMAP(object = so, feature = "Prrx1")
#'
exprUMAP <- function(object, feature, plot.subtitle = NULL, do.neb = F, title.size = 10, slot = "data", assay = DefaultAssay(object), reduction = "umap", size = autoPointSize(ncol(object)),scale.color = "tomato",  ...){

  require(RColorBrewer)

  DefaultAssay(object) <- assay
  if (!do.neb){
    scExpression.UMAP(object,feature, adjust.pt.size =size, slot = slot, reduction = reduction,...) +
      theme_void() + theme(legend.position = "none") +
      scale_color_gradient(low = "grey95", high = scale.color) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      labs(subtitle = plot.subtitle) +
      theme(plot.subtitle=element_text(size=5, color="black")) +
      theme(plot.title =element_text(size=title.size, face="italic", color="black"))
  } else {
    Nebulosa::plot_density(object,feature, size = size, slot = slot,reduction = reduction, ... ) +
      theme_void() + theme(legend.position = "none") +
      scale_color_gradient(low = "grey95", high = scale.color) +
      theme(plot.title = element_text(hjust = 0.5)) +
      theme(plot.subtitle = element_text(hjust = 0.5)) +
      labs(subtitle = plot.subtitle) +
      theme(plot.subtitle=element_text(size=5, color="black")) +
      theme(plot.title =element_text(size=title.size, face="italic", color="black"))
  }


}



#' Draw volcano plot to visualize differential expression.
#'
#' Draw volcano plot to visualize differential expression. Uses data.frame output from getDEG(..., return.list = F).
#'
#' @param df.deg Differential expression data. Dataframe output from presto::wilcoxauc() or scMiko::getDEG(..., return.list = F).
#' @param group group to visaulize data for. Must be entry in 'group' column of df.deg.
#' @param show.n Top n features to label. Default is 10.
#' @param show.sig.only Show significant genes only. Default is T.
#' @param rank.by statistic to rank features by when selecting top feature to label. Default is "auc" (options include "auc", "logFC", "pval").
#' @param features Specific features to label. If specified, show.n is ignored.
#' @param sig.threshold significance (p-value) threshold. Default is 0.05.
#' @param lfc.threshold log fold change (lfc) threhsold. Default is 0.
#' @param label.size Label size.
#' @param pt.size.range vector with two values, the first representing the smallest point size, and the second representing the largest point size. Default is c(1,4).
#' @param correct.p.value use adjusted p-value (FDR) for thresholding significant hits. Default is T.
#' @param plot.title title of volcano plot. Default is "Volcano Plot".
#' @param cols named listed with 3 entries labeled "low", "mid", "high" specifying colors to create a color gradient. Default is list(low = scales::muted("blue"), mid = "white", high = scales::muted("red")).
#' @seealso getDEG
#' @name miko_volcano
#' @return ggplot handle
#' @examples
#'
#'  df.dat <- getDEG(so.query, return.list = F)
#'  plt.volcano <-  miko_volcano(df.deg = df.dat)
#'
miko_volcano <- function(df.deg, group = NULL, show.n = 10, show.sig.only = T, rank.by = c("auc", "logFC", "pval"),
                         features = NULL, sig.threshold = 0.05, lfc.threshold = 0, label.size = NA, pt.size.range = c(1,4),
                         correct.p.value = T, plot.title = "Volcano Plot", cols = list(low = scales::muted("blue"), mid = "white", high = scales::muted("red"))){

  # assertions
  if (!("avgExpr" %in% colnames(df.deg))) {
    df.deg$avgExpr <- 1
    size.label <- NULL
  } else {
    size.label <- "avgExpr"
  }
  stopifnot(all(c("group", "avgExpr", "logFC", "pval",  "feature") %in% colnames(df.deg)))

  if (ulength(df.deg$group) > 1){

    stopifnot(!is.null(group))
    stopifnot(sum(grepl(group, df.deg$group)) > 0)
    # filter
    df.deg$deg.group <- df.deg$group
    df.deg <- df.deg %>% dplyr::select(-c("group"))
    df.deg <- df.deg %>% dplyr::filter(grepl(pattern = group, x = deg.group))
  }

  # get genes to show
  if (!(rank.by %in% colnames(df.deg))){
    rank.by <- c("auc", "pval", "logFC", "avgExpr")[which(c("auc", "pval", "logFC", "avgExpr") %in% colnames(df.deg))[1]]
  }
  if (rank.by %in% "pval"){
    df.deg$z <- (-log10(unlist(df.deg[ ,rank.by]))) * sign(df.deg$logFC)
  } else if (rank.by %in% "logFC"){
    df.deg$z <- df.deg[ ,"logFC"]
  } else if (rank.by %in% "auc"){
    df.deg$z <- df.deg[ ,"auc"]
  } else if (rank.by %in% "avgExpr"){
    df.deg$z <- unlist(df.deg[ ,"avgExpr"])  * sign(df.deg$logFC)
  } else if (rank.by %in% "pct_in"){
    df.deg$z <- unlist(df.deg[ ,"pct_in"])  * sign(df.deg$logFC)
  }

  df.deg$padj <- p.adjust(df.deg$pval, method = "BH")

  deg.top <- NULL
  if (show.sig.only){
    if (correct.p.value){
      df.deg.sig <- df.deg %>% dplyr::filter(padj < sig.threshold, abs(logFC) > lfc.threshold)
    } else {
      df.deg.sig <- df.deg %>% dplyr::filter(pval < sig.threshold, abs(logFC) > lfc.threshold)
    }
  } else {
    df.deg.sig <- df.deg
    lfc.threshold <- 0
  }
  if (is.null(features)){
    deg.top <- bind_rows(df.deg.sig %>% dplyr::top_n(round(show.n/2), z),
                         df.deg.sig %>% dplyr::top_n(round(show.n/2), -z))
  } else {
    deg.top <- df.deg.sig %>% dplyr::filter(feature  %in% features)
  }

  if (is.null(nrow(deg.top)) | nrow(deg.top) == 0){
    deg.top <- bind_rows(df.deg.sig %>% dplyr::top_n(round(10/2), logFC),
                         df.deg.sig %>% dplyr::top_n(round(10/2), logFC))
  }


  if (correct.p.value){
    y.axis.label <- "-log10(FDR)"
    plabel <- "FDR"

  } else {
    y.axis.label <- "-log10(p)"
    deg.top$padj <- deg.top$pval
    df.deg$padj <- df.deg$pval
    plabel <- "p"
  }


  # color by...
  color.by <- "auc"
  if (!(color.by %in% colnames(df.deg))){
    color.by <- c("auc", "pval", "logFC")[which(c("auc", "pval", "logFC") %in% colnames(df.deg))[1]]
  }
  if (color.by %in% "pval"){
    df.deg$w <- (-log10(unlist(df.deg[ ,color.by]))) * sign(df.deg$logFC)
    color.label <- "signed\n-log(p)"
  } else if (color.by %in% "logFC"){
    df.deg$w <- df.deg[ ,"logFC"]
    color.label <- "logFC"
  } else if (color.by %in% "auc"){
    df.deg$w <- abs(df.deg[ ,"auc"] - 0.5)
    color.label <- "|AUC-0.5|"
  }

  # generate plot
  if (nrow(deg.top)> 0)  deg.top$feature <- paste0("italic('", deg.top$feature, "')")

  nup <- sum(df.deg$padj <= sig.threshold & df.deg$logFC > lfc.threshold)
  ndown <- sum(df.deg$padj <= sig.threshold & df.deg$logFC < -lfc.threshold)

  if (lfc.threshold > 0){
    sig.label  <- paste0(nup, " upregulated, ", ndown , " downregulated (", plabel, "<", sig.threshold, ", |logFC|>", lfc.threshold, ")")
  } else {
    sig.label  <- paste0(nup, " upregulated, ", ndown , " downregulated (", plabel, "<", sig.threshold, ")")
  }


  plt.vol <-  ggplot() +
    geom_point(data = df.deg %>% dplyr::filter(padj <= sig.threshold, abs(logFC) >= lfc.threshold),
               aes(x = logFC, y = -log10(padj), color = w, size = abs(avgExpr))) +
    geom_point(data = df.deg %>% dplyr::filter(padj > sig.threshold | abs(logFC) < lfc.threshold),
               aes(x = logFC, y = -log10(padj), size = abs(avgExpr)), color = "grey") +
    scale_size(range = pt.size.range) +
    ggrepel::geom_text_repel(data = deg.top, aes(x = logFC, y = -log10(padj), label = feature),
                             parse = T, size = label.size, min.segment.length = 0, max.overlaps = Inf) +
    geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed") +
    geom_vline(xintercept = 0) +
    scale_color_gradient2(low = cols$low, mid = cols$mid, high = cols$high) +
    theme_miko(legend = T) +
    labs(title = plot.title,
         subtitle = sig.label,
         color = color.label, size = size.label, y = y.axis.label)


  if (lfc.threshold > 0){
    plt.vol <- plt.vol + geom_vline(xintercept = c(-lfc.threshold, lfc.threshold), linetype = "dashed")
  }


  plt.vol

}


#' Gradient color scale
#'
#' Wrapper for scale_color_gradient2, with preferred default parameters
#'
#' @param low color for low end of gradient
#' @param mid color for middle point of gradient
#' @param high  color for high end of gradient
#' @param ... additional arguments passed to scale_color_gradient2(...)
#' @name scale_color_miko
#' @examples
#'
#'  df.dat <- getDEG(so.query, return.list = F)
#'  plt.volcano <-  miko_volcano(df.deg = df.dat) + scale_color_miko()
#'
scale_color_miko <- function(low = scales::muted("blue"), high = scales::muted("red"), mid = "white", ...){
  scale_color_gradient2(low = low, high = high, mid=  mid)
}

#' Gradient fill scale
#'
#' Wrapper for scale_fill_gradient2, with preferred default parameters
#'
#' @param low color for low end of gradient
#' @param mid color for middle point of gradient
#' @param high  color for high end of gradient
#' @param ... additional arguments passed to scale_color_gradient2(...)
#' @name scale_fill_miko
#'
scale_fill_miko <- function(low = scales::muted("blue"), high = scales::muted("red"), mid = "white", ...){
  scale_fill_gradient2(low = low, high = high, mid=  mid)
}

#' Highlight cells on UMAP plot
#'
#' Highlight cells on UMAP plot
#'
#' @param object Seurat object.
#' @param group grouping variable in object's meta data.
#' @param reduction reduction used to visualize cells. Default is "umap".
#' @param highlight.color color of highlighted cells. Default is "tomato".
#' @name highlightUMAP
#' @return list of ggplot handles
#' @examples
#'
#'  plot.highlight <-  highlightUMAP(object = so, group = "seurat_clusters")
#'
highlightUMAP <- function(object, group = "seurat_clusters", reduction = "umap", highlight.color = "tomato"){


  df.red.group <- data.frame(object@reductions[[reduction]]@cell.embeddings)
  colnames(df.red.group) <- c("x", "y")

  stopifnot("group is not specified in object metadata" = group %in% colnames(object@meta.data))
  df.red.group$cluster <- object@meta.data[[group]]
  plt.cluster.umap <- list()

  cluster.lab <- group
  cluster.membership <- object@meta.data[[cluster.lab]]
  #
  # # get unique clusters
  if (group == "seurat_clusters"){
    u.clusters <- unique(as.numeric(as.character((cluster.membership))))
    u.clusters <- u.clusters[order(u.clusters)]
  } else {
    u.clusters <- unique((as.character((cluster.membership))))
  }


  for (i in 1:length(u.clusters)){

    cluster.name <- u.clusters[i]

    df.red.group$do.color <- "grey"
    df.red.group$do.color[df.red.group$cluster %in% cluster.name] <- highlight.color

    ncell <- sum(df.red.group$cluster %in% cluster.name)

    df.red.group <- df.red.group %>% dplyr::arrange(do.color)

    plt.cluster.umap[[paste0("group_", cluster.name)]] <- df.red.group %>%
      ggplot(aes(x = x, y = y)) +
      geom_point(color = df.red.group$do.color, size = autoPointSize(nrow(df.red.group))) +
      theme_miko() +
      theme_void() +
      xlab("UMAP 1") + ylab("UMAP 2") +
      labs(title = "UMAP", subtitle = paste0("Group: ", cluster.name, " (", ncell, " cells)"))


  }


  return(plt.cluster.umap)

}
