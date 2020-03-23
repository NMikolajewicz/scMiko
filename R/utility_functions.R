#' Uppercase first letter and lowercase rest
#'
#' Function used to naively convert gene list to mouse representation; i.e., Convert Hs to Mm gene symbols.
#'
#'
#' @param x Character vector
#' @name firstup
#' @return Character vector
#'
firstup <- function(x) {
  x <- tolower(x)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


#' Check gene representation
#'
#' Check whether genes are represented as SYMBOL or ENSEMBL. This is accomplished by comparing query to reference, and using a majority vote.
#'
#' @param reference.genes Named vector of genes; names are ENSEMBL, entries are SYMBOL.
#' @param query.genes vector of query genes to check representation
#' @name checkGeneRep
#' @return Character specifying ensembl or symbol
#'
checkGeneRep <- function(reference.genes, query.genes){

  ensembl_rep <- sum(as.vector(names(reference.genes)) %in% query.genes)
  symbol_rep <- sum(toupper(as.vector(reference.genes)) %in% toupper(query.genes))

  if (ensembl_rep > symbol_rep){
    gene.rep <- "ensembl"
  } else if (symbol_rep > ensembl_rep){
    gene.rep <- "symbol"
  } else {
    gene.rep <- NA
  }

  return(gene.rep)

}


#' Convert gene ensemble to symbol in seurat object
#'
#' Process Seurat object so that genes in Ensemble format are converted to Symbol representation.
#'
#' @param reference.genes Named vector of genes; names are ENSEMBL, entries are SYMBOL.
#' @param query.genes vector of query genes to check representation
#' @name ens2sym.so
#' @return Character specifying ensembl or symbol
#'
ens2sym.so <- function(so, gNames.list, convert.RNA = TRUE){

    # var features
    so_temp <- so@assays[["SCT"]]@var.features
    so@assays[["SCT"]]@var.features <- as.vector((gNames.list[so_temp]))

    # scale data
    so_temp <- so@assays[["SCT"]]@scale.data
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@assays[["SCT"]]@scale.data <- so_temp

    # data
    so_temp <- so@assays[["SCT"]]@data
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@assays[["SCT"]]@data <- so_temp

    # metadata
    so_m <- so@assays[["SCT"]]@meta.features
    rownames(so_m) <-  make.names(as.vector((gNames.list[rownames(so_m)])), unique = T)
    so@assays[["SCT"]]@meta.features <- so_m

    # dimnames
    so_temp <- so@assays[["SCT"]]@counts@Dimnames[[1]]
    so@assays[["SCT"]]@counts@Dimnames[[1]] <- as.vector((gNames.list[so_temp]))

    # pca feature loading
    so_temp <-  so@reductions[["pca"]]@feature.loadings
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@reductions[["pca"]]@feature.loadings <- so_temp

    # RNA ASSAY

    if (convert.RNA == TRUE){
      # var features
      so_temp <- so@assays[["RNA"]]@var.features
      so@assays[["RNA"]]@var.features <- as.vector((gNames.list[so_temp]))

      # data
      so_temp <- so@assays[["RNA"]]@data
      row.names(so_temp) <-   so@assays[["RNA"]]@meta.features[["gene_name"]]
      so@assays[["RNA"]]@data <- so_temp

      # counts
      so_temp <- so@assays[["RNA"]]@counts
      row.names(so_temp) <-   so@assays[["RNA"]]@meta.features[["gene_name"]]
      so@assays[["RNA"]]@counts <- so_temp

      # Dimnames
      so@assays[["RNA"]]@counts@Dimnames[[1]] <- so@assays[["RNA"]]@meta.features[["gene_name"]]
    }

    if (DefaultAssay(so) == "integrated"){

      # var features
      so_ens <- so@assays[["integrated"]]@var.features
      so@assays[["integrated"]]@var.features <-as.vector((gNames.list[so_ens]))

      # scale data
      so_sd <- so@assays[["integrated"]]@scale.data
      row.names(so_sd) <-  as.vector((gNames.list[row.names(so_sd)]))
      so@assays[["integrated"]]@scale.data <- so_sd

      # data
      so_d <- so@assays[["integrated"]]@data
      row.names(so_d) <-  as.vector((gNames.list[row.names(so_d)]))
      so@assays[["integrated"]]@data <- so_d
    }


    return(so)
}


#' Convert gene symbol representation to Hs or Mm
#'
#' Convert gene symbol to approriate species (as specified), while ensuring that gene is present in reference gene list.
#'
#' @param query.gene Character. Query gene symbol
#' @param gene.list Character(s). Vector of all genes under consideration
#' @param which.species Character. Species
#' @name speciesConvert
#' @return Character specifying ensembl or symbol
#'
speciesConvert <- function(query.gene, gene.list, which.species){

  if (which.species == "Hs"){
    if (toupper(query.gene) %in% gene.list){
      cur.marker <- toupper(query.gene)
    } else if (firstup(query.gene) %in% gene.list){
      cur.marker <- firstup(query.gene)
    } else {
      cur.marker <- NA
    }
  } else {
    if (firstup(query.gene) %in% gene.list){
      cur.marker <- firstup(query.gene)
    } else if (toupper(query.gene) %in% gene.list){
      cur.marker <- toupper(query.gene)
    } else {
      cur.marker <- NA
    }

  }

  return(cur.marker)
}


#' Convert gene symbol to entrez
#'
#' Gene symbol is converted to entrez id using org.Hs.eg.db or org.Mm.eg.db annotation databases.
#'
#' @param my.symbols Character. Vector of Gene symbols.
#' @param my.species Character. Species.
#' @name sym2entrez
#' @return data.frame mapping gene Symbols to Entrez
#'
sym2entrez <- function(so, gNames.list, convert.RNA = TRUE){

  my.symbols <- as.vector(my.symbols)
  if (my.species == "Hs"){
    db <- org.Hs.eg.db
  } else if (my.species == "Mm"){
    db <- org.Mm.eg.db
  }

  my.entrez <- AnnotationDbi::select(db,
                                     keys = my.symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")

  return(my.entrez)
}


#' Create heatmap object
#'
#' Uses gplots::heatmap.2 function to generate heatmap object. Run in chunk specifying include = F to suppress plot and only output object.
#'
#' @param mat Matrix. Input matrix for heatmap.
#' @param hmcol Heatmap colors
#' @param scale.limit Numeric. Color scale limit.
#' @name sym2entrez
#' @return data.frame mapping gene Symbols to Entrez
#'
getHeat <- function(mat, hmcol, scale.limit){

  #threw this into a function to suppress plot that is automatically generated
  heat.object <- NULL
  try({
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      breaks =seq(-scale.limit, scale.limit, by = ((2*scale.limit/100))),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = "Gene Exp Matrix",
                                      xlab = "Cluster ID")

  }, silent = T)

  if (is.null(heat.object)){
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = "Gene Exp Matrix",
                                      xlab = "Cluster ID")
  }

  return(heat.object)
}
