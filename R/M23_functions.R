
#' Stratify gene expression into high and low expression groups
#'
#' Creates new metadata field stratifying gene expression (for queried gene) into high and low expresion groups. If return.plt is true, histogram illustrating threshold is generated; Otherwise seurat object with new metadata entry is returned.
#'
#' @param so Seurat Object
#' @param query.gene Character specifying gene of interest.
#' @param quartile.threshold Quantile used to threshold expression data. E.g., 0.25 means that <25th percentile and >75th percentile are flagged as low and high expression groups, respecitvely.
#' @param return.plt Logical indicating whether histogram with thresholds is returned. Otherwise Seurat object is returned.
#' @name m23.binarizeExpression
#' @return Seurat object (unless return.plt is T, in which case ggplot handle)
#'
m23.binarizeExpression <- function(so, query.gene, quantile.threshold = 0.25, return.plt = F){


  which.match <- which(rownames(so@assays[[DefaultAssay(so)]]) %in% query.gene)

  df.exp <- data.frame(data = so@assays[[DefaultAssay(so)]]@data[which.match, ],
                       scale = so@assays[[DefaultAssay(so)]]@scale.data[which.match, ])

  df.exp$cellID <- seq(1, ncol(so@assays[[DefaultAssay(so)]]@data))

  df.exp.filtered <- df.exp[df.exp$data > median(df.exp$data), ]

  df.exp.filtered$rank <- percent_rank(df.exp.filtered$data)

  df.exp.filtered$expression.grp <- "other"
  df.exp.filtered$expression.grp[df.exp.filtered$rank > (1-quantile.threshold)] <- "high"
  df.exp.filtered$expression.grp[df.exp.filtered$rank < (quantile.threshold)] <- "low"

  expression.group <- paste(query.gene, ".Expression", sep  = "")
  so@meta.data[[expression.group]] <- "other"
  so@meta.data[[expression.group]][df.exp.filtered$cellID[df.exp.filtered$expression.grp == "high"]] <- "high"
  so@meta.data[[expression.group]][df.exp.filtered$cellID[df.exp.filtered$expression.grp == "low"]] <- "low"


  if (return.plt){
    plt.hist.filtered <-   df.exp.filtered %>%
      ggplot(aes(x = data)) +
      geom_histogram(aes(fill =  df.exp.filtered$expression.grp)) +
      xlab("Expression") + ggtitle(paste(query.gene, "Expression")) +
      scale_fill_manual(values = c( "red", 'blue', 'grey') )  +
      theme_classic()  + labs(fill = "Expression")

    return(plt.hist.filtered)
  } else {
    return(so)
  }


}
