#' Uppercase first letter and lowercase rest
#'
#' Function used to naively convert gene list to mouse representation; i.e., Convert Hs to Mm gene symbols.
#'
#'
#' @param x Character vector
#' @name firstup
#' @seealso \code{\link{toupper}}
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
#' @seealso \code{\link{ens2sym.so}}
#' @author Nicholas Mikolajewicz
#' @return Character specifying ensembl or symbol
#' @examples
#'
#' \dontrun{
#' gene.rep <-  checkGeneRep(gNames.list, as.vector(rownames(so.query@assays[[current.assay]]@scale.data)))
#'
#' if (gene.rep == "ensembl"){
#'   so.query <- ens2sym.so(so = so.query, gNames.list = gNames.list, convert.RNA = TRUE)
#'   gene.rep <-  checkGeneRep(gNames.list, as.vector(rownames(so.query@assays[[current.assay]]@scale.data)))
#' }
#' }
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
#' @author Nicholas Mikolajewicz
#' @return Seurat object with specified gene representation.
#'
ens2sym.so <- function(so, gNames.list, convert.RNA = TRUE){

  # var features
  so_temp <- so@assays[["SCT"]]@var.features
  gene.rep <- checkGeneRep (gNames.list, so_temp)
  if (gene.rep == "ensembl") so@assays[["SCT"]]@var.features <- as.vector((gNames.list[so_temp]))

  # scale data
  so_temp <- so@assays[["SCT"]]@scale.data
  gene.rep <- checkGeneRep (gNames.list, row.names(so_temp))
  if (gene.rep == "ensembl") {
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@assays[["SCT"]]@scale.data <- so_temp
  }

  # data
  so_temp <- so@assays[["SCT"]]@data
  gene.rep <- checkGeneRep (gNames.list, row.names(so_temp))
  if (gene.rep == "ensembl") {
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@assays[["SCT"]]@data <- so_temp
  }

  # metadata
  so_m <- so@assays[["SCT"]]@meta.features
  so_m$ENSEMBLE <- rownames(so_m)
  gene.rep <- checkGeneRep (gNames.list,  so_m$ENSEMBLE)
  if (gene.rep == "ensembl") {
    so_m$SYMBOL <- as.vector((gNames.list[so_m$ENSEMBLE]))
    rownames(so_m) <-  make.names(so_m$SYMBOL, unique = T)
    so@assays[["SCT"]]@meta.features <- so_m
  }

  # dimnames
  so_temp <- so@assays[["SCT"]]@counts@Dimnames[[1]]
  gene.rep <- checkGeneRep (gNames.list, so_temp)
  if (gene.rep == "ensembl") {
    so@assays[["SCT"]]@counts@Dimnames[[1]] <- as.vector((gNames.list[so_temp]))
  }

  # pca feature loading
  if ("pca" %in% names(so@reductions)){
    so_temp <-  so@reductions[["pca"]]@feature.loadings
    gene.rep <- checkGeneRep (gNames.list, row.names(so_temp))
    if (gene.rep == "ensembl") {
      row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
      so@reductions[["pca"]]@feature.loadings <- so_temp
    }
  }

  # vst
  gene.rep <- checkGeneRep (gNames.list,  names(so@assays[["SCT"]]@misc[["vst.out"]][["genes_log_gmean_step1"]]))
  if (gene.rep == "ensembl") {
    names(so@assays[["SCT"]]@misc[["vst.out"]][["genes_log_gmean_step1"]]) <- as.vector((gNames.list[names(so@assays[["SCT"]]@misc[["vst.out"]][["genes_log_gmean_step1"]])]))
  }

  if (!is.null(so@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]])){
    gene.rep <- checkGeneRep (gNames.list,   so@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]]@Dimnames[[1]])
    if (gene.rep == "ensembl") {
      so@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]]@Dimnames[[1]] <- as.vector((gNames.list[ so@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]]@Dimnames[[1]]]))
    }
  }


  gene.rep <- checkGeneRep (gNames.list,  rownames(so@assays[["SCT"]]@misc[["vst.out"]][["gene_attr"]]))
  try({
    if (gene.rep == "ensembl") {
      rownames(so@assays[["SCT"]]@misc[["vst.out"]][["gene_attr"]]) <- as.vector((gNames.list[rownames(so@assays[["SCT"]]@misc[["vst.out"]][["gene_attr"]])]))
    }
  }, silent = T)

  try({
    so.temp <- rownames(so@assays[["SCT"]]@misc[["vst.out"]][["model_pars"]])
    gene.rep <- checkGeneRep (gNames.list, so.temp)
    if (gene.rep == "ensembl") {
      rownames(so@assays[["SCT"]]@misc[["vst.out"]][["model_pars"]]) <- as.vector((gNames.list[so.temp]))
    }
  }, silent = T)

  try({
    so.temp <- rownames(so@assays[["SCT"]]@misc[["vst.out"]][["model_pars_fit"]])
    gene.rep <- checkGeneRep (gNames.list,  so.temp)
    if (gene.rep == "ensembl") {
      rownames(so@assays[["SCT"]]@misc[["vst.out"]][["model_pars_fit"]]) <- as.vector((gNames.list[so.temp]))
    }
  }, silent = T)

  # RNA ASSAY

  if (convert.RNA == TRUE){
    # var features
    so_temp <- so@assays[["RNA"]]@var.features
    if (length(so_temp) > 0){
      gene.rep <- checkGeneRep (gNames.list, so_temp)
      if (gene.rep == "ensembl") {
        so@assays[["RNA"]]@var.features <- as.vector((gNames.list[so_temp]))
      }
    }


    # data
    so_temp <- so@assays[["RNA"]]@data
    gene.rep <- checkGeneRep (gNames.list, row.names(so_temp))
    if (gene.rep == "ensembl") {
      row.names(so_temp) <-   as.vector((gNames.list[ row.names(so_temp)]))
      so@assays[["RNA"]]@data <- so_temp
    }

    # counts
    so_temp <- so@assays[["RNA"]]@counts
    gene.rep <- checkGeneRep (gNames.list, row.names(so_temp))
    if (gene.rep == "ensembl") {
      row.names(so_temp) <-   as.vector((gNames.list[ row.names(so_temp)]))
      so@assays[["RNA"]]@counts <- so_temp
    }

  }

  if ("integrated" %in% names(so@assays)){

    # var features
    so_ens <- so@assays[["integrated"]]@var.features
    gene.rep <- checkGeneRep (gNames.list, so_ens)
    if (gene.rep == "ensembl") {
      so@assays[["integrated"]]@var.features <-as.vector((gNames.list[so_ens]))
    }

    # scale data
    so_sd <- so@assays[["integrated"]]@scale.data
    gene.rep <- checkGeneRep (gNames.list,  row.names(so_sd))
    if (gene.rep == "ensembl") {
      row.names(so_sd) <-  as.vector((gNames.list[row.names(so_sd)]))
      so@assays[["integrated"]]@scale.data <- so_sd
    }

    # data
    so_d <- so@assays[["integrated"]]@data
    gene.rep <- checkGeneRep (gNames.list,  row.names(so_d))
    if (gene.rep == "ensembl") {
      row.names(so_d) <-  as.vector((gNames.list[row.names(so_d)]))
      so@assays[["integrated"]]@data <- so_d
    }

  }

  # ensure dim names are correctly specified
  so <- updateDimNames(so)


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
#' @author Nicholas Mikolajewicz
#' @return Gene Symbol (character)
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
#' @author Nicholas Mikolajewicz
#' @return data.frame mapping gene Symbols to Entrez
#'
sym2entrez <- function(my.symbols, my.species){

  my.symbols <- as.vector(my.symbols)
  if (my.species == "Hs"){
    db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (my.species == "Mm"){
    db <- org.Mm.eg.db::org.Mm.eg.db
  }

  my.entrez <- AnnotationDbi::select(db,
                                     keys = my.symbols,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "SYMBOL")

  return(my.entrez)
}

#' Convert gene symbol to ensembl
#'
#' Gene symbol is converted to ensembl id using org.Hs.eg.db or org.Mm.eg.db annotation databases.
#'
#' @param my.symbols Character. Vector of Gene symbols.
#' @param my.species Character. Species.
#' @name sym2ens
#' @author Nicholas Mikolajewicz
#' @return data.frame mapping gene Symbols to Ensembl
#'
sym2ens <- function(my.symbols, my.species){

  my.symbols <- as.vector(my.symbols)
  if (my.species == "Hs"){
    db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (my.species == "Mm"){
    db <- org.Mm.eg.db::org.Mm.eg.db
  }

  my.ensembl <- AnnotationDbi::select(db,
                       keys = my.symbols,
                       columns = c("ENSEMBL", "SYMBOL"),
                       keytype = "SYMBOL",
                       multiVals = first)

  return(my.ensembl)
}


#' Create heatmap object
#'
#' Uses gplots::heatmap.2 function to generate heatmap object. Run in chunk specifying include = F to suppress plot and only output object.
#'
#' @param mat Matrix. Input matrix for heatmap.
#' @param hmcol Heatmap colors
#' @param scale.limit Numeric. Color scale limit.
#' @param main Character. Name of plot.
#' @param xlab Character. X axis label.
#' @param ... additional parameters passed to gplots::heatmap.2(...)
#' @name getHeat
#' @author Nicholas Mikolajewicz
#' @return Heatmap object
#'
getHeat <- function(mat, hmcol = NULL, scale.limit = NULL, main = NULL, xlab = NULL, ...){

  # get color palette
  if (is.null(hmcol)) hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100)

  # get scale limits
  if (is.null(scale.limit)) scale.limit <- max(abs(mat))

  # get plot title
  if (is.null(main)) main <- "Gene Exp Matrix"

  # get x label
  if (is.null(xlab)) xlab <- "Cluster ID"

  #threw this into a function to suppress plot that is automatically generated
  heat.object <- NULL
  try({
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      breaks =seq(-scale.limit, scale.limit, by = ((2*scale.limit/100))),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = main,
                                      xlab = xlab,
                                      ...)

  }, silent = T)

  if (is.null(heat.object)){
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = "Gene Exp Matrix",
                                      xlab = "Cluster ID",
                                      ...)
  }

  return(heat.object)
}


#' Get list of available files
#'
#' Get list of files in specified directory. Default directory is "Preprocessed Datasets/".
#'
#' @param directory Character. Directory name.
#' @name getAvailableFiles
#' @author Nicholas Mikolajewicz
#' @return list of files
#'
getAvailableFiles <- function(directory = "Preprocessed Datasets/"){
return(list.files(directory))
}


#' Normalize and scale data within appropriate assay
#'
#' Ensures data are properly normalized and scaled. If intergrated dataset provided, Seurat's NormalizeData, FindVariableFeatures, ScaleData workflow is applied and default assay is set to 'RNA'. Otherwise, if  non-integrated dataset is provided, SCT transform is assumed and default assay is set to 'SCT'. If non-intergrated data have not been process with SCT workflow, see m1.scNormScale() function.
#'
#'
#' @param so Seurat Object
#' @name prepExpression
#' @author Nicholas Mikolajewicz
#' @return Seurat Object
#'
prepExpression <- function(so){

  if (DefaultAssay(so) == "integrated"){
    DefaultAssay(so) <- "RNA"
    so <-NormalizeData(so, verbose = FALSE)
    so <- ScaleData(so, verbose = FALSE)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 10000)
  } else {
    if ("SCT" %in% names(so@assays)) DefaultAssay(so) <- "SCT"
  }

  return(so)

}




#' Check if gene is avaialble in Seurat Object
#'
#' Return logical flag indicated whether query gene is present in seurat expression matrix.
#'
#' @param so Seurat Object
#' @param query.gene Character. Gene of interest.
#' @param reference.genes Named vector of genes; names are ENSEMBL, entries are SYMBOL.
#' @name isGeneAvailable
#' @author Nicholas Mikolajewicz
#' @return Logical
#'
isGeneAvailable <- function(so, query.gene, reference.genes){
  all.genes <- rownames(so@assays[[DefaultAssay(so)]])
  gene.rep <-  checkGeneRep(reference.genes, all.genes)
  stopifnot(gene.rep == "symbol")
  geneAvailable <- query.gene %in% all.genes

  return(geneAvailable)
}



#' Reload scMiko package
#'
#' Function that detaches and attaches scMiko package.
#'
#' @name scMikoReload
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # reload scMiko package
#' scMikoReload()
#'
scMikoReload <- function(){

  try({detach("package:scMiko", unload=TRUE)}, silent = T)
  library(scMiko)

}



#' Merge list of seurat objects
#'
#' Merges list of seurat objects without any normalization of batch correction
#'
#' @param so.list List of seurat objects
#' @name mergeSeuratList
#' @author Nicholas Mikolajewicz
#' @return Seurat Object
#'
mergeSeuratList <- function(so.list){

  # merge seurat objects if multiple present
  if (length(so.list) > 1){
    for (i in 1:length(so.list)){
      if (i == 1){
        so <- so.list[[i]]
      } else {
        so <- merge(so, so.list[[i]])
      }
    }
  } else {
    so <- so.list[[1]]
  }

  return(so)
}




#' Fix barcode labels
#'
#' Rename metadata entry 'CellTypes' to 'Barcode'. This is a fix implemented to correct an error from an earlier analysis pipeline. Analyses post January 2020 do not require this fix.
#'
#'
#' @param so Seurat Object
#' @name fixBarcodeLabel
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#'
fixBarcodeLabel <- function (so){
  # merge CellType and Barcode, if necessary
  meta.data.names <- names(so@meta.data)

  if (("CellType" %in% meta.data.names) & ("Barcode" %in% meta.data.names)){
    if (DefaultAssay(so) == "integrated"){
      barcode <- so@meta.data[["Barcode"]]
      celltype <- so@meta.data[["CellType"]]
      barcode[is.na(barcode)] <- celltype[is.na(barcode)]
    } else {
      barcode <- so@meta.data[["CellType"]]
    }
  } else if (!("CellType" %in% meta.data.names) & ("Barcode" %in% meta.data.names)) {
    barcode <- so@meta.data[["Barcode"]]
  } else if (("CellType" %in% meta.data.names) & !("Barcode" %in% meta.data.names)) {
    barcode <- so@meta.data[["CellType"]]

  } else {stop("Problem with CellType/Barcode metadata detected. Troubleshooting required")}

  so@meta.data[["Barcode"]] <- barcode

  return(so)
}

#' Set cluster resolution
#'
#' Set 'Seurat_Clusters' metadata entry to specified cluster resolution [0, inf]. See Seurat::FindClusters() for details.
#'
#'
#' @param so Seurat Object
#' @param cluster.resolution Numeric [0, inf] specifying cluster resolution. Values [0.1,1] typically perform well.
#' @name setResolution
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#'
setResolution <- function (so, cluster.resolution){

  so <- FindClusters(object = so, resolution = cluster.resolution, verbose = 0, algorithm = 1)

  return(so)
}


#' prep Gene List
#'
#' Ensure is available and represented correctly.
#'
#' @param so Seurat Object
#' @param global.enviroment global.enviroment <- objects()
#' @name prepGeneList
#' @author Nicholas Mikolajewicz
#' @return Named vector of available genes
#' @examples
#'
#' gNames.list <- prepGeneList(so, objects())
#'
prepGeneList <- function (so, global.enviroment){

  # global.enviroment <- object()


  if (("gNames.list_master" %in% global.enviroment)){
    gNames.list <- NULL
    for (i in 1:length(gNames.list_master)){
      gNames.list <- c(gNames.list, gNames.list_master[[i]] )
    }

    gNames.df <-  data.frame(n = gsub("\\..*","",as.vector(names(gNames.list))), g = as.vector(gNames.list))
    gNames.df <- unique(gNames.df)
    gNames.list <- as.vector(gNames.df$g)
    names(gNames.list) <- as.vector(gNames.df$n)
  }  else {

    # check if gene-ensemble pair are present in meta-data
    av.meta <- so@assays[["RNA"]]@meta.features

    if (all(c("SYMBOL", "ENSEMBL") %in% colnames(av.meta))){
      gNames.list <- as.vector(av.meta$SYMBOL)
      names(gNames.list) <- as.vector(av.meta$ENSEMBL)
    }

  }

  # ensure gene list is available
  stopifnot(exists("gNames.list"))

  return(gNames.list)
}


#' Return load path
#'
#' For specified file and directory, return a load path
#'
#'
#' @param file file name
#' @param directory directory
#' @name getLoadPath
#' @return Character. Load path.
#'
getLoadPath <- function (file, directory = NULL){
  if (is.null(directory)) directory <- ""
  load.path <- paste(directory, file, sep = "")
  return(load.path)
}

#' prep Seurat
#'
#' Preprocess using fixBarcodeLabel(), UpdateSeuratObject() and updateDimNames() Functions.
#'
#' @param so Seurat objects
#' @name prepSeurat
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#'
prepSeurat <- function (so){

  if (class(so) != "Seurat") stop("input must be Seurat Object")

  so <- fixBarcodeLabel(so)

  if (!("SCT_nn" %in% names(so@graphs)))  {
    so <- UpdateSeuratObject(so) # required after Seurat 3.1.2 update
  }
  so <- updateDimNames(so) # required after Seurat 3.1.2 update

  return(so)
}


#' Subset Seurat Object
#'
#' Subset Seurat object according to specific metadata field. Only specified metadata entries are retained, while remaining of data is omitted.
#'
#'
#' @param so Seurat object
#' @param subset.df Data.frame specifying which field (subset.df$field) to subset on, and which field entries to retain (subset.df$subgroups).
#' @name subsetSeurat
#' @return Seurat object
#' @examples
#'
#' # define subset parameters
#' subset.df <- data.frame(field = "seurat_clusters", subgroups = c(0,2,3,5,6,7,9,15)) # tumor population from Renca 4000ds T12 sample
#'
#' # subset data
#' so <- subsetSeurat(so, subset.df)
#'
subsetSeurat <- function (so, subset.df){

  # check if subset input is validd
  if (is.na(unique(subset.df$field))){
    subset.flag <- FALSE
  } else if ( unique(subset.df$field) %in% names(so@meta.data)) {
    subset.flag <- TRUE
  } else {
    subset.flag <- FALSE
  }


  # subset data
  if (subset.flag){
    pattern <- paste( as.vector(subset.df$subgroups), collapse="|")
    pattern <- gsub(" ", "", pattern)
    cur.field <- as.vector(unique(subset.df$field))
    match.ind <- grepl(pattern, as.character(so@meta.data[[cur.field]]))
    so <- subset(x = so, cells = which(match.ind))
    so <- UpdateSeuratObject(so)
  }

  return(so)
}


#' Returns list of annotations for given Entrez gene IDs
#'
#' Returns list of Reactome or GO annotations for given Entrez gene IDs
#'
#' @param query.genes Entrez IDs of query genes
#' @param db Database to retrieve annotations from. One of:
#' \itemize{
#' \item "Bader" - Default. List of pathway annotations curated by Bader lab (http://baderlab.org/GeneSets)
#' \item "Reactome"
#' \item "GO"
#' }
#' @param ontology GO ontologies to retrieve if GO db is selected. One of:
#' \itemize{
#' \item "BP" - Default. Biological processes
#' \item "MF" - Molecular functions
#' \item "CC" - Cellular components
#' }
#' @param species "Mm" or "Hs". Default is "Hs".
#' @name getAnnotationPathways
#' @return Named list of vectors with gene sets (Entrez format).
#'
getAnnotationPathways <- function(query.genes, db = c("Bader"), ontology = c("BP"), species = c("Hs")){

  if (db == "GO"){

    which.ontology <- ontology

    if (species == "Hs"){
      go.e2g <- org.Hs.egGO
      go.g2e <- as.list(org.Hs.egGO2EG)
    } else if (species == "Mm"){
      go.e2g <- org.Mm.egGO
      go.g2e <- as.list(org.Mm.egGO2EG)
    }

    # Get the entrez gene identifiers that are mapped to a GO ID
    mapped_genes <- mappedkeys(go.e2g)

    # Convert to a list
    entrez2go.list <- as.list(go.e2g[mapped_genes])

    # get matching GO terms
    which.match <- names(entrez2go.list) %in% query.genes

    # get GO ids that overlap with geneset
    entrez2go.subset <- unique(unlist(lapply(entrez2go.list[which.match], names)))

    # map GO id to term
    annots <-  AnnotationDbi::select(GO.db::GO.db, keys=entrez2go.subset, columns=c("GOID", "TERM", "ONTOLOGY"), keytype="GOID")

    keep.which <- c()
    for (i in 1:length(which.ontology)){
      stopifnot(which.ontology[i] %in% c("BP", "MF", "CC"))
      keep.which <- c(keep.which, annots$ONTOLOGY %in% which.ontology[i])
    }
    annots.ontology <- annots[keep.which, ]


    # filter genesets according to whats left
    go.g2e.filtered <- go.g2e[names(go.g2e) %in% annots.ontology$GOID]
    g2t.map <- annots.ontology$TERM
    names(g2t.map) <- annots.ontology$GOID
    names(go.g2e.filtered) <- g2t.map[names(go.g2e.filtered)]

    pathways <- go.g2e.filtered

  } else if (db == c("Reactome")){
    pathways <- reactomePathways(query.genes)

  } else if (db == c("Bader")){
    pathway.list <- scMiko::baderPathways

    if (species == "Hs"){
      pathways <- pathway.list[["Hs.entrez"]]
    } else if (species == "Mm"){
      pathways <- pathway.list[["Mm.entrez"]]
    }

    include.which <- lapply(pathways, function(x) sum(query.genes %in% unlist(x)) > 0 )
    pathways <- pathways[as.vector(unlist(include.which))]

  }

  return(pathways)
}


#' Map Reactome/GO term to ID
#'
#' Map Reactome/GO term to ID
#'
#' @param term Reactome/GO term
#' @param db Database to use for mapping. One of:
#' \itemize{
#' \item "Reactome" - Default
#' \item "GO"
#' }
#' @param species Character specifying species. If specified, one of:
#' \itemize{
#' \item "Hs"
#' \item "Mm"
#' }
#' @name term2id
#' @return Reactome/GO ID
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # map GO term to GO id
#' id <- term2id(term, db = "GO", species = "Hs")
#'
#' # get geneset for specified GO id (returns character of genes)
#' gene.set <- id2geneset(id, "Hs")
#'
term2id <- function(term, db = "Reactome", species = NULL){

  if (db == "GO"){

    annots <-  AnnotationDbi::select(GO.db::GO.db, keys=term, columns=c("GOID", "TERM", "ONTOLOGY"), keytype="TERM")
    id <- annots$GOID

  } else if (db == "Reactome"){

    if (is.null(species)) stop("Must specify species")

    if (species == "Hs"){
      if (!(grepl("Homo sapiens: ", term))){
        term.Hs <- paste("Homo sapiens: ",term, sep = "")
      } else {
        term.Hs <- term
      }
      annots <- AnnotationDbi::select(reactome.db::reactome.db, keys=term.Hs, columns=c("PATHID","PATHNAME"), keytype="PATHNAME")
    } else if (species == "Mm"){
      if (!(grepl("Mus musculus: ", term))){
        term.Mm <- paste("Mus musculus: ",term, sep = "")
      } else {
        term.Mm <- term
      }
      annots <- AnnotationDbi::select(reactome.db::reactome.db, keys=term.Mm, columns=c("PATHID","PATHNAME"), keytype="PATHNAME")
    }
    id <- annots$PATHID
    return(id)
  }
}


#' Map Reactome/GO ID to term
#'
#' Map Reactome/GO ID to term
#'
#' @param term Reactome/GO ID
#' @param db Database to use for mapping. One of:
#' \itemize{
#' \item "Reactome" - Default
#' \item "GO"
#' }
#' @param species Character specifying species. If specified, one of:
#' \itemize{
#' \item "Hs"
#' \item "Mm"
#' }
#' @name id2term
#' @return Reactome/GO term
#'
id2term <- function(id, db = "Reactome", species = NULL){

  if (db == "GO"){

    annots <-  AnnotationDbi::select(GO.db::GO.db, keys=id, columns=c("GOID", "TERM", "ONTOLOGY"), keytype="GOID")
    term <- annots$TERM

  } else if (db == "Reactome"){

    annots <- AnnotationDbi::select(reactome.db::reactome.db, keys=id, columns=c("PATHID","PATHNAME"), keytype="PATHID")
    term.out <- annots$PATHNAME

    if (!is.null(species)){
      if (species == "Hs"){
        term <- strsplit(term.out, "Homo sapiens: ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
      } else if (species == "Mm"){
        term <- strsplit(term.out, "Mus musculus: ", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][2]
      }
    } else {
      term <- term.out
    }
  }

  return(term)

}


#' Search Reactome/GO databases for terms that match query
#'
#' Searches Reactome/GO databases using query stem (i.e., incomplete term) and returns all matches.
#'
#' @param query Query term. A character.
#' @param db Database to search, if specified. One of:
#' \itemize{
#' \item "Reactome"
#' \item "GO"
#' }
#' @param species Character specifying species, if specified. One of:
#' \itemize{
#' \item "Hs"
#' \item "Mm"
#' }
#' @param ontology Character specifying GO ontology to filter terms by. If unspecified, one or more of:
#' \itemize{
#' \item "BP" - biological processes
#' \item "MF" - molecular functions
#' \item "CC" - cellular compartment
#' }
#' @name searchAnnotations
#' @author Nicholas Mikolajewicz
#' @return dataframe of terms that match query
#'
searchAnnotations <- function(query, db = NULL, species = NULL, ontology = NULL){


  if (is.null(db)) db <- c("Reactome", "GO")

  if ( "Reactome" %in% db){

    term.id.list <- as.list(reactome.db::reactomePATHNAME2ID)

    # filter by query
    term.id.names <- names(term.id.list)
    which.match <- grepl(query, term.id.names)
    term.id.list.match <- term.id.list[which.match]

    # filter by species
    term.id.names <- names(term.id.list.match)

    if (is.null(species)){
      term.id.list.match <- term.id.list.match
    } else if (species == "Mm"){
      which.match <- grepl("Mus musculus: ", term.id.names)
      term.id.list.match <-  term.id.list.match[which.match]
    } else if (species == "Hs"){
      which.match <- grepl("Homo sapiens: ", term.id.names)
      term.id.list.match <-  term.id.list.match[which.match]
    }

    query.match.reactome <- data.frame(id = unlist(term.id.list.match), term = names(term.id.list.match))
    rownames(query.match.reactome) <- c()

  }

  if ("GO" %in% db){

    term.id.list <- as.list(GO.db::GOTERM)

    go.term.list <- lapply(term.id.list, function(x) x@Term)
    go.ontology.list <- lapply(term.id.list, function(x) x@Ontology)

    go.df <- data.frame(id = names(go.term.list), term = unlist(go.term.list), ontology = unlist(go.ontology.list))
    go.df$ontology <- as.character(go.df$ontology)

    go.df.filter <- go.df[grepl(query, go.df$term), ]

    if (!is.null(ontology)){
      for(i in 1:length(ontology)){
        go.df.filter <- go.df.filter[(go.df.filter$ontology %in% ontology[i]), ]
      }
    }

    query.match.go <- go.df.filter
    rownames(query.match.go) <- c()

  }


  if ("KEGG" %in% db){

    term.id.list <- as.list(GO.db::GOTERM)

    go.term.list <- lapply(term.id.list, function(x) x@Term)
    go.ontology.list <- lapply(term.id.list, function(x) x@Ontology)

    go.df <- data.frame(id = names(go.term.list), term = unlist(go.term.list), ontology = unlist(go.ontology.list))
    go.df$ontology <- as.character(go.df$ontology)

    go.df.filter <- go.df[grepl(query, go.df$term), ]

    if (!is.null(ontology)){
      for(i in 1:length(ontology)){
        go.df.filter <- go.df.filter[(go.df.filter$ontology %in% ontology[i]), ]
      }
    }

    query.match.go <- go.df.filter
    rownames(query.match.go) <- c()

  }

  if (exists("query.match.go") & exists("query.match.reactome") ){
    q1 <- data.frame(db = "GO", id = query.match.go$id, term = query.match.go$term, ontology = query.match.go$ontology)
    q2 <- data.frame(db = "Reactome", id = query.match.reactome$id, term = query.match.reactome$term, ontology = NA)
    q.final <- bind_rows(q1, q2)
  } else if (exists("query.match.go") & (!exists("query.match.reactome")) ){
    q.final <- data.frame(db = "GO", id = query.match.go$id, term = query.match.go$term, ontology = query.match.go$ontology)
  } else if ((!exists("query.match.go")) & exists("query.match.reactome") ){
    q.final <- data.frame(db = "Reactome", id = query.match.reactome$id, term = query.match.reactome$term, ontology = NA)
  }

  return(q.final)
}



#' Get Reactome/GO geneset
#'
#' Retrieves Reactome/GO geneset using specified annotation ID
#'
#' @param id Reactome/GO identifier. Reactome has "R-" prefix; GO has "GO" prefix.
#' @param my.species Specify species; ensures that correct gene symbols are retrieved.
#' @name id2geneset
#' @author Nicholas Mikolajewicz
#' @return Character vector of gene symbols belonging to Reactome/GO geneset.
#'
id2geneset <- function(id, my.species){

  if (grepl("GO", id)){
    db = "GO"
  } else if (grepl("R-", id)){
    db = "Reactome"
  } else {
    stop("ID is invalid. Check input.")
  }

  if (db == "GO"){

    if (my.species == "Hs"){
      my.entrez.Hs <- revmap(org.Hs.egGO)[[id]]
      my.symbol <- entrez2sym(my.entrez = my.entrez.Hs, my.species = "Hs")
    } else if (my.species == "Mm"){
      my.entrez.Mm <- revmap(org.Mm.egGO)[[id]]
      my.symbol <- entrez2sym(my.entrez = my.entrez.Mm, my.species = "Mm")
    }

  } else if (db == "Reactome"){

    reactome.list <- as.list(reactomePATHID2EXTID)
    reactome.query <- reactome.list[names(reactome.list) %in% id]
    stopifnot(length(reactome.query) == 1)
    my.entrez <- unlist(reactome.query)

    if (my.species == "Hs"){
      my.symbol <- entrez2sym(my.entrez = my.entrez, my.species = "Hs")
    } else if (my.species == "Mm"){
      my.symbol <- entrez2sym(my.entrez = my.entrez, my.species = "Mm")
    }

  }

  gene.set <- my.symbol$SYMBOL

  return(gene.set)

}



#' Convert entrez id to gene symbol
#'
#' Entrez ID is converted to gene symbol using org.Hs.eg.db or org.Mm.eg.db annotation databases.
#'
#' @param my.symbols Character. Vector of Entrez IDs.
#' @param my.species Character. Species.
#' @name entrez2sym
#' @return data.frame mapping Entrez IDs to gene Symbols.
#'
entrez2sym <- function(my.entrez, my.species){

  my.entrez <- as.vector(my.entrez)
  if (my.species == "Hs"){
    db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (my.species == "Mm"){
    db <- org.Mm.eg.db::org.Mm.eg.db
  }

  my.symbol <- AnnotationDbi::select(db,
                                     keys = my.entrez,
                                     columns = c("ENTREZID", "SYMBOL"),
                                     keytype = "ENTREZID")

  return(my.symbol)
}


#' Get expression matrix from Seurat Object
#'
#' Get expression matrix from Seurat Object
#'
#' @param so Seurat Object
#' @param only.variable Logical indicating whether to include variable features only or not.
#' @param which.assay Seurat assay to get data frome. Default is DefaultAssay(so).
#' @param which.data Specify which data to use (refers to slots in Seurat object assay). One of:
#' \itemize{
#' \item "scale" - Default
#' \item "data"
#' }
#' @param use.additional.genes Character vector of additional genes to include (in addition to varibale, if variable flag is specificed). Default is NA.
#' @name getExpressionMatrix
#' @author Nicholas Mikolajewicz
#' @return gene x cell expression matrix
#'
getExpressionMatrix <- function(so, only.variable = F, which.assay = NULL, which.data = "scale", use.additional.genes = NA){

  # specify assay
  if (is.null(which.assay)) which.assay <- DefaultAssay(so)

  # get complete matrix
  if (which.data == "scale"){
    exp.mat.complete <- so@assays[[which.assay]]@scale.data
  } else if (which.data == "data"){
    exp.mat.complete <- as.matrix(so@assays[[which.assay]]@data)
  }



  if (only.variable){
    var.feat <-  VariableFeatures(so)
    if (length(var.feat) == 0){
      if (which.assay == "SCT"){
        if ("integrated" %in% names(so@assays)){
          var.feat <-  so@assays[["integrated"]]@var.features
        } else {
          var.feat <-  so@assays[["SCT"]]@var.features
        }
      } else {
        stop("Could not find variable features in Seurat Object. ")
      }
    }

    exp.mat <- exp.mat.complete[rownames(exp.mat.complete) %in% var.feat, ]
  } else {
    exp.mat <- exp.mat.complete
  }


  if (!is.na(use.additional.genes)){
    which.missing <- which(!(use.additional.genes %in% rownames(exp.mat)))
    missing.additional.genes <- (use.additional.genes)[which.missing]
    exp.mat.additional <- exp.mat.complete[rownames(exp.mat.complete) %in% missing.additional.genes, ]
    exp.mat <- rbind(exp.mat, exp.mat.additional)
  }

  return(exp.mat)
}



#' Remove "ï.." prefix that is appended to csv header
#'
#' Remove "ï.." prefix that is appended to csv header
#'
#' @param x Character vector
#' @name rmvCSVprefix
#' @return Character vector
#'
rmvCSVprefix <- function(x){
  pattern <- "ï.."
  x <- stringr::str_replace(x, pattern, "")
  return(x)
}

#' Get cells that express query gene
#'
#' Get cells that express query gene above specified threshold
#'
#' @param so Seurat Object
#' @param query Query gene. A character.
#' @param expression.threshold Numeric. Return cells that express query above threshold value. Default is 0.
#' @param which.data Seurat data slot ("data", "scale")
#' @name getExpressingCells
#' @author Nicholas Mikolajewicz
#' @return Character vector of cell barcodes/ids
#'
getExpressingCells <- function(so, query, expression.threshold = 0, which.data = "data"){
  exp.data.all<- getExpressionMatrix(so, which.data = which.data)
  exp.data.query <- exp.data.all[rownames(exp.data.all) %in% query, ]
  expressing.cells <- names(exp.data.query)[exp.data.query > expression.threshold]

  return(expressing.cells)
}




#' Create/update list of gene sets for scMiko package.
#'
#' Takes genesets stored in Excel sheets, and converts them to dataframes stored in lists.
#'
#' @param input.file Excel file (input). Must have ".xlsx" suffix. A character.
#' @param output.file Rdata file (output). Must have ".rda" suffix. A character.
#' @param dir Directory of input and output file (same folder). A character.
#' @param dev.directory.flag Logical indicating whether to use developer specific director. Default is False. If true dir is ignored.
#' @name updateGeneSets
#' @author Nicholas Mikolajewicz
#' @return List of data.frames saved as Rdata file.
#' @examples
#' \dontrun{
#' updateGeneSets("geneSets_MASTER_260420update.xlsx", "geneSets.rda", dev.directory.flag = T)
#' }
#'
updateGeneSets <- function(input.file, output.file, dir = "", dev.directory.flag = F){

  geneSets <- list()

  if (dev.directory.flag){
    input.dir <- "D:/Users/Nick/Dropbox/PDF Projects - JM/Data/scRNA-seq/01_sci-RNA-seq3_Hong_Kevin_Jason/Nick/scRNAseq/Reference Datasets/"
    output.dir <- "D:/Users/Nick/Dropbox/PDF Projects - JM/R Packages/scMiko/data/"
  } else {
    input.dir <- dir
    output.dir <- dir
  }

  sheetNames <- openxlsx::getSheetNames(getLoadPath(input.file, input.dir))

  for (i in 1:length(sheetNames)){
    current.sheet <-  sheetNames[i]
    gene.set <- readxl::read_excel(getLoadPath(input.file, input.dir), sheet = current.sheet)

    geneSets[[current.sheet]] <- gene.set
  }


  save(geneSets, file=getLoadPath(output.file, output.dir))

  return(paste(length(geneSets), " genesets successfully saved to '", getLoadPath(output.file, output.dir), "'", sep= ""))

}




#' Clean and filter gene list
#'
#' Clean gene list (i.e., remove NAs), convert to appropraite species, and keep only those present in seurat object.
#'
#' @param genes Character vector of genes
#' @param so Seurat Object
#' @param which.species Species. Either 'Hs' or 'Mm'.
#' @name cleanFilterGenes
#' @seealso \code{\link{speciesConvert}}
#' @author Nicholas Mikolajewicz
#' @return Character vector of genes
#'
cleanFilterGenes <- function(genes, so, which.species){

  # clean dataset and include only those available in seurat object
  cur.features <- genes
  cur.features <- cur.features[!is.na(cur.features)]
  cur.features <- lapply(cur.features,  speciesConvert, rownames(so@assays[[DefaultAssay(so)]]@data), which.species)
  cur.features <-as.vector(unlist(cur.features))
  cur.features <- cur.features[!is.na(cur.features)]

  return(cur.features)

}


#' Downsample single cell data
#'
#' Downsample number of cells in Seurat object by specified factor
#'
#' @param so Seurat Object
#' @param subsample.factor Numeric [0,1]. Factor to downsample data by.
#' @param subsample.n Numeric [1,ncol(so)]. Number of cells to subsample. If specified, overides subsample.factor.
#' @name downsampleSeurat
#' @author Nicholas Mikolajewicz
#' @return Seurat Object
#'
downsampleSeurat <- function(so, subsample.factor = 1, subsample.n = NULL){


  if (subsample.factor < 1){
    if (subsample.factor<0) stop("subsample.factor must be numeric between 0 and 1")

    so <-  tryCatch({

      if (is.null(subsample.n)){
        n.subset <- round(subsample.factor *ncol(so))
      } else {
        if (subsample.n > ncol(so)) {
          warning(paste("subsample.n exceeded number of cells in seurat object. Data was not subsampled"))
          subsample.n <- ncol(so)
        }
        n.subset <- subsample.n
      }

      cell.ind <- sample(x = seq(1, ncol(so)), size = n.subset, replace = FALSE, prob = NULL)
      so <- subset(so , cells = cell.ind)
    }, error = function(e){
      warning("Failed to downsample seurat object")
      return(so)
    })

  }

  return(so)

}

#' Get vector of unique ordered group names from Seurat Object
#'
#' Get vector of unique ordered group names from Seurat Object
#'
#' @param so Seurat Object
#' @param which.group Character specfying group field in Seurat metadata. Default is "seurat_clusters".
#' @param is.number Logical specifying whether groups IDs are numeric. Default is True.
#' @name getOrderedGroups
#' @return Ordered vector
#'
getOrderedGroups <- function(so, which.group = "seurat_clusters", is.number = T){

  groups <- so@meta.data[[which.group]]

  if (is.number){
    u.groups <- as.numeric(as.character(unique(groups)))
  } else {
    u.groups <- as.character(unique(groups))
  }

  u.groups <- u.groups[order(u.groups)]

  return(u.groups)

}

#' Get summary of group expression in Seurat object
#'
#' Get summary group expression in Seurat object. Can include mean, median, fraction (of expressing cells), sd, or cv. Calls scMiko::avgGroupExpression().
#'
#' @param so Seurat Object
#' @param which.data Character specfying which data slot. Default is "data".
#' @param which.assay Character specigin which assay to use.
#' @param which.center Character indicating which summary measure to use. Must be one of "mean", "median", "fraction", "sum", "sd", or "cv". If unspecified, default is "mean".
#' @param which.group Character specfying group field in Seurat metadata. Default is "seurat_clusters".
#' @param do.parallel Logical specifying whether to perform computations in parallel. Default is F. Uses future.apply package.
#' @name aggGroupExpression
#' @author Nicholas Mikolajewicz
#' @return data.frame (gene rows, group columns)
#'
aggGroupExpression <-  function(so, which.data = "data", which.assay = DefaultAssay(so), which.center = "mean", which.group = "seurat_clusters", do.parallel = F){
  return(avgGroupExpression(so, which.data, which.assay, which.center, which.group, do.parallel))
}

#' Get summary of group expression in Seurat object
#'
#' Get summary group expression in Seurat object. Can include mean, median, fraction (of expressing cells), sd, or cv.
#'
#' @param so Seurat Object
#' @param which.data Character specfying which data slot. Default is "data".
#' @param which.assay Character specigin which assay to use.
#' @param which.center Character indicating which summary measure to use. Must be one of "mean", "median", "fraction", "sum", "sd", or "cv". If unspecified, default is "mean".
#' @param which.group Character specfying group field in Seurat metadata. Default is "seurat_clusters".
#' @param do.parallel Logical specifying whether to perform computations in parallel. Default is F. Uses future.apply package.
#' @name avgGroupExpression
#' @author Nicholas Mikolajewicz
#' @return data.frame (gene rows, group columns)
#'
avgGroupExpression <-  function(so, which.data = "data", which.assay = DefaultAssay(so), which.center = "mean", which.group = "seurat_clusters", do.parallel = F){
  # which.center options: "mean", "fraction", "median", "sum", "sd", "cv"

  # inititate parallel processes
  if (do.parallel){
    library(future.apply)
    plan(multisession) ## Run in parallel on local computer
  }

  # entire matrix
  exp.mat.complete <- getExpressionMatrix(so, which.data = which.data, which.assay = which.assay)

  # group ID vector
  cluster.membership <- so@meta.data[[which.group]]

  # gene list
  gene.list <- rownames(exp.mat.complete)

  # ordered vector of unique groups
  u.clusters <- getOrderedGroups(so, which.group, is.number = F)

  if ((which.center == "fraction") & (which.data != "data")){
    warning("\nData from 'data' slot used to compute expressing fraction")
    which.data <- "data"
    exp.mat.complete <- getExpressionMatrix(so, which.data = which.data, which.assay = which.assay)
    gene.list <- rownames(exp.mat.complete)
  }

  # clear some memory
  rm(so)

  # compute measure of centrality
  avg.mat <- matrix(nrow = length(gene.list), ncol = length(u.clusters))
  for (i in 1:length(u.clusters)){
    warning("\nComputing measures of centrality...")
    if (which.center == "mean"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x, na.rm = T))
      } else {
        avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x, na.rm = T))
      }
    } else if (which.center == "median"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) median(x, na.rm = T))
      } else {
        avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) median(x, na.rm = T))
      }
      } else if (which.center == "sum"){
        if (do.parallel){
          avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sum(x, na.rm = T))
        } else {
          avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sum(x, na.rm = T))
        }
      } else if (which.center == "sd"){
        if (do.parallel){
          avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x, na.rm = T))
        } else {
          avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x, na.rm = T))
        }
      } else if (which.center == "cv"){
        if (do.parallel){
          sd.cur <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x, na.rm = T))
          av.cur <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x, na.rm = T))
        } else {
          sd.cur <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x, na.rm = T))
          av.cur <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x, na.rm = T))
        }
        avg.mat[,i] <- sd.cur / abs(av.cur)
      } else if (which.center == "fraction"){
        if (do.parallel){
          avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sum(x>0)/length(x))
        } else {
          avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sum(x>0)/length(x))
        }
      } else {
        stop("which.center must be specified as 'mean', 'median', 'fraction', 'sd', or 'cv'")
      }
    }

    df.avg <- as.data.frame(avg.mat)

    if (which.group == "seurat_clusters"){
      colnames(df.avg) <- paste("c", u.clusters, sep = "")
    } else {
      colnames(df.avg) <- u.clusters
    }

    df.avg <- bind_cols(data.frame(genes = gene.list), df.avg)

    return(df.avg)


  }




#' Over-representation enrichment of GO terms using weighted fisher method
#'
#' Over-representation enrichment of GO terms using weighted fisher method
#'
#' @param gene.candidates Vector of candidate genes.
#' @param gene.universe Vector of all genes.
#' @param which.species Character specifying species. "Mm" or "Hs".
#' @param which.ontology Character specifying ontology. "BP", "MF", or "CC".
#' @param p.threshold p value threshold. Default is 0.001.
#' @param padj.threshold adjusted p value threhsold (BH). Default is 0.05.
#' @param topGO.object topGO object (Optional). If unspecified, new topGO object is new created. If specified, gene list is used to updated existing object.
#' @name enrichGO.fisher
#' @import topGo
#' @return list of 2 data.frames (unadjusted.results (p<0.001) and adjusted.results (padj < 0.05)) and a topGO object.
#'
enrichGO.fisher <- function(gene.candidates, gene.universe, which.species , which.ontology = "BP", p.threshold = 0.001,
                            padj.threshold = 0.05, topGO.object = NULL){

  if (!(which.species %in% c("Hs", "Mm"))) stop("Species incorrectly specified. Must be either Hs or Mm")

  if (which.species == "Hs"){
    library(org.Hs.eg.db)
    db <- "org.Hs.eg.db"
  } else if (which.species == "Mm"){
    library(org.Mm.eg.db)
    db <- "org.Mm.eg.db"
  }


  allGO2genes <- topGO::annFUN.org(whichOnto=which.ontology, feasibleGenes=NULL, mapping=db, ID="symbol")
  # all.genes <- colnames(datExpr)

  # make named factor showing which genes are of interest
  geneList=factor(as.integer(gene.universe %in% gene.candidates))
  names(geneList)= gene.universe

  if (is.null(topGO.object)){
    # if no GO object provided, create new
    GOdata <- new("topGOdata",
                  ontology=which.ontology,
                  allGenes=geneList,
                  annot=topGO::annFUN.GO2genes,
                  GO2genes=allGO2genes,
                  nodeSize=10)
  } else {
    # if GO object exists, update gene list
    GOdata <- topGO::updateGenes(topGO.object, geneList)
  }

  # define test using the weight01 algorithm (default) with fisher
  res.wfisher <- runTest(GOdata, algorithm='weight01', statistic='fisher')

  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO <- usedGO(GOdata)
  res.table <- GenTable(GOdata, weightFisher=res.wfisher, orderBy='weightFisher', topNodes=length(allGO))

  # ensure p values are numeric
  res.table$weightFisher <- as.numeric(res.table$weightFisher)

  #performing BH correction on our p values
  p.adj <- round(p.adjust(res.table$weightFisher,method="BH"),digits = 4)

  # create the file with all the statistics from GO analysis
  res.table.adj <- cbind(res.table,p.adj)
  res.table.adj=res.table.adj[order(res.table.adj$p.adj),]

  #get list of significant GO before multiple testing correction
  results.table.p <-  res.table.adj[which(res.table.adj$weightFisher<=p.threshold),]

  #get list of significant GO after multiple testing correction
  results.table.bh <- res.table.adj[which(res.table.adj$p.adj<=padj.threshold),]


  return(list(unadjusted.results = results.table.p, adjusted.results = results.table.bh, topGo.object = GOdata))
}




#' Filter seurat object according to specified meta data entries
#'
#' Using a mapping.list, entries from an existing meta data field are relabeled and mapped to new meta data field, and cells which are not specified in this mapping are omitted from the seurat object.
#'
#' @param so Seurat object.
#' @param old.field Existing seurat meta field that will be relabelled.
#' @param new.field New seurat meta field that will be created with relablled entries from old.field.
#' @param mapping.list Named list specifying how to relabel entries in old.field (mapping.list entries) to new.field (mapping.list names) in seurat meta data. Note that entries in mapping list do not have to match old.field entires exactly; entries are used as pattern arguemtn to grepl() function.
#' @name mapSubsetSeurat
#' @author Nicholas Mikolajewicz
#' @return Seurat object in which old field was mapped to new field, according to mapping specified in mapping.list. Any cells that are not mapped are omitted from seurat object.
#' @examples
#'
#' # define mapping lists
#' mapping.list.1 <- list(NSG = "NSG", BALBc = "BALB")
#' mapping.list.2 <- list(Early = "Early", Mid = "Mid", Late = "Late")
#'
#' # Create new field called "Mouse" using relabled entries from "Condition" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query, old.field = "Condition", new.field = "Mouse", mapping.list = mapping.list.1)
#'
#' # Create new field called "Time" using relabled entries from "Group" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query.2, old.field = "Group", new.field = "Time", mapping.list = mapping.list.2)
#'
mapSubsetSeurat <- function(so, old.field, new.field, mapping.list){

  keep.this <- NULL
  for (i in 1:length(mapping.list)){

    cur.group <- names(mapping.list)[order(names(mapping.list))][i]
    cur.pattern <- (mapping.list)[order(names(mapping.list))][[i]]

    if (is.null(keep.this)){
      keep.this <- grepl(cur.pattern, so@meta.data[[old.field]])
    } else {
      keep.this <- (keep.this | grepl(cur.pattern, so@meta.data[[old.field]]))
    }

    so@meta.data[[new.field]][grepl(cur.pattern, so@meta.data[[old.field]])] <- cur.group

  }


  # ensure groups are ordered and that order is maintained throughout analysis.
  u.groups <- as.character(unique(so@meta.data[[new.field]]))
  u.groups <- u.groups[order(u.groups)]
  so@meta.data[[new.field]] <- factor(so@meta.data[[new.field]], levels = u.groups)
  so <-so[, keep.this]

  return(so)
}


#' Perform differential expression analysis
#'
#' Give two existing fields within a Seurat Object, data are stratified by the first, and pairwise comparisons between groups in the second are performed to identify differential gene expression. Uses Seurat's FindMarkers() function.
#'
#' @param so Seurat object.
#' @param ordered.levels Vector of one or two seurat meta fields. If one provided, the second is set to "seurat_clusters".
#' @param which.assay Seurat assay to use. If unspecified, set to DefaultAssay(so).
#' @param ... Additional arguments passed to Seurat's FindMarkers() function.
#' @name multiLevel.FindMarkers
#' @return list of data.frames where each data.fram contains results from an individual pairwise comparison.
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # define mapping lists
#' mapping.list.1 <- list(NSG = "NSG", BALBc = "BALB")
#' mapping.list.2 <- list(Early = "Early", Mid = "Mid", Late = "Late")
#'
#' # Create new field called "Mouse" using relabled entries from "Condition" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query, old.field = "Condition", new.field = "Mouse", mapping.list = mapping.list.1)
#'
#' # Create new field called "Time" using relabled entries from "Group" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query.2, old.field = "Group", new.field = "Time", mapping.list = mapping.list.2)
#'
#' # Perform multilevel differential gene expression.
#' # Data are stratified by time, and then pairwise comparison between "Mouse" groups are performed.
#' # in this example, 3 total pairwise comparisons were made: BALBc vs NSG (Early), BALBc vs NSG (Mid) and BALBc vs NSG (Late)
#' deg.list <- multiLevel.FindMarkers(so.query.2, ordered.levels = c("Time", "Mouse"),
#'                                    logfc.threshold = 0,
#'                                    min.pct = 0,
#'                                    test.use = "wilcox")
#'
multiLevel.FindMarkers <- function(so, ordered.levels, which.assay = NULL, ...){

  if (length(ordered.levels) == 1) ordered.levels <- c(ordered.levels, "seurat_clusters")
  if (is.null(which.assay)) which.assay <- DefaultAssay(so)

  # get unique level factors
  level.1 <- unique(as.character(so@meta.data[[ordered.levels[1]]]))
  level.2 <-  unique(as.character(so@meta.data[[ordered.levels[2]]]))

  # set idents to level.2
  Idents(so) <- ordered.levels[2]

  # get all level 2 combinations
  level.2.combinations <- t(combinat::combn(level.2, 2))

  #N comparisons
  n.comparisons <- ncol(level.2.combinations) * length(level.1)

  #initiate results list
  deg.list <- list()
  for (i in 1:length(level.1)){

    # get current level 1 strata
    level.1.cur <- level.1[i]

    # subset seurat according to level 1
    keep.this <- as.character(so@meta.data[[ordered.levels[1]]]) %in% level.1.cur
    so.cur <-so[, keep.this]

    for (j in 1:nrow(level.2.combinations)){

      # define comparison level 2 comparison group
      group1 <- level.2.combinations[j, 1]
      group2 <- level.2.combinations[j, 2]

      comparison.label <- paste(level.1.cur, "|" , group1, ".vs.", group2, sep = "")

      suppressWarnings({
        deg.list[[comparison.label]] <- FindMarkers(so.cur,
                                                    assay = which.assay,
                                                    ident.1 = group1,
                                                    ident.2 = group2,
                                                    verbose = F,
                                                    ...)

      })

    }

  }

  return(deg.list)

}




#' Gene connectivity within network.
#'
#' Computes gene connectivity from adjaceny/topological overlap matix. Part of WCGNA network analysis workflow.
#'
#' @param w.mat Adjaceny or toplogoical overlap matrix. Symmetrical matrix, with row and column-wise gene entries.
#' @param gene.names Vector of gene names corresponding to genes in coloumns/rows of w.mat.
#' @param flag.top.n Numerical specifying top N genes to flag (for subsequent plotting). Default is 20.
#' @name getConnectivity
#' @return Dataframe with gene connectivity, rank and label flag.
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # run WCGNA
#' output.all <- runWCGNA(datExpr.noz, cor.metric = "rho_p", soft.power = 2, use.TOM = T)
#'
#' # unpack output
#' s.mat <- output.all[["s.mat"]] # similar matrix
#' a.mat <- output.all[["a.mat"]] # adjacency matrix
#' w.mat <- output.all[["w.mat"]] # topological overlap matix
#' d.mat <- output.all[["d.mat"]] # disimilarity matix
#'
#' # get connectivity
#' df.con <- getConnectivity(w.mat, gene.names = colnames(a.mat))
#'
getConnectivity <- function(w.mat, gene.names, flag.top.n = 20){

  # compute connectivity and store in dataframe
  df.con <- data.frame(genes = gene.names, wi = apply(w.mat, 1, function(x) sum(x)))

  # rank
  df.con$rank <- rank(df.con$wi)

  # flag top n genes to label
  df.con$label <- F
  df.con$label[df.con$rank > nrow(df.con)-flag.top.n] <- T
  df.con$label[df.con$rank < flag.top.n] <- T

  return(df.con)

}



#' run WGCNA analysis on scRNAseq expression matrix
#'
#' Run WGCNA analysis on scRNAseq expression matrix, using WGCNA R package.
#'
#' @param e.mat Expression matrix. Row entries are cells, column entries are genes. Colnames and rownames are expected.
#' @param s.mat Similarity matrix (optional). If not provided, will be computed using method specfiied by cor.metric. If provided, s.mat is not recomputed.
#' @param cor.metric Correction measure to use. Default is "rho_p". See "dismay" package for additional options.
#' @param soft.power Soft power used to scale s.mat to a.mat (e.g., a.mat = s.mat ^ soft.power)
#' @param use.TOM Logical flag specifying whether to compute topoligical overlap matrix. If false, w.mat = a.mat.
#' @param network.type Network type. Allowed values are (unique abbreviations of) "unsigned", "signed" (default), "signed hybrid"
#' @param TOM.type TOM type. Allowed values are "unsigned" (default) or "signed"
#' @param rescale.adjacency Logical indicate whether adjacency matrix is rescaled to [0,1]. Default is False.
#' @param ... Additional arguments passessed to TOMsimilarity {WGCNA package}
#' @name runWGCNA
#' @return List containing  similarity matrix (s.mat), adacency matrix (a.mat), topological overlap matrix (w.mat) and disimilarity matrix (d.mat)
#' @examples
#'
#' # Get expression matrix
#' which.data <- "scale"
#'
#' # variable gene only matrix
#' use.var <- T
#' if (use.var){
#'   exp.mat <- getExpressionMatrix(so.query, only.variable = use.var, which.data = which.data, use.additional.genes = NA)
#' } else {
#'   exp.mat <- exp.mat.complete
#' }
#'
#'
#' # transpose expressio matrix (genes are columns)
#' t.exp.mat <- t(exp.mat)
#' datExpr <- as.matrix(t.exp.mat)
#' SubGeneNames=colnames(datExpr)
#'
#' # capture output used to hide undesired print statement
#' print2hide <- capture.output(allowWGCNAThreads())
#'
#' # transform matrix if necessary
#' if (min(datExpr) < 0) {
#'   datExpr.noz <- datExpr + abs(min(datExpr))
#' } else {
#'   datExpr.noz <- datExpr
#' }
#'
#' # run WGCNA
#' output.all <- runWGCNA(datExpr.noz, cor.metric = "rho_p", soft.power = 2, use.TOM = T)
#'
runWGCNA <- function(e.mat, s.mat = NULL, cor.metric = "rho_p", soft.power = 2, use.TOM = T, network.type = "signed", TOM.type = "unsigned", rescale.adjacency = F, ...){

  # similarity matrix - using proportionality metric for scRNAseq data.
  if (is.null(s.mat)){
    warning("Computing similarity matrix...")
    s.mat <-  dismay::dismay(e.mat, metric = cor.metric)
  }

  # adjacency matrix
  warning("\nComputing adjacency matrix...")
  a.mat <-  sim2adj(s.mat, soft.power, network.type)

  # rescale value if needed
  if (rescale.adjacency) a.mat <- recaleValues(a.mat, new.min = 0, new.max = 1)

  # compute topological overlap matix (TOM)
  if (use.TOM){
    warning("\nComputing topological overlap matix...")
    if ((TOM.type) == "signed" & (network.type == "unsigned")) {
      a.mat.tom <- a.mat * sign(s.mat)
    } else {
      a.mat.tom <- a.mat
    }

    # ensure matrix is symmetric (override symmetry checking performed by TOMsimiliary())
    if (isSymmetric(a.mat.tom)){
      a.mat.tom <- (a.mat.tom + t(a.mat.tom)) / 2
    } else {
      stop("Adjacency matrix is not symmetric")
    }

    # ensure adjacency matrix satisfied criterion
    a.mat.tom <- 0.9999*a.mat.tom

    print2hide <-  capture.output(w.mat <- TOMsimilarity(a.mat.tom, TOMType = TOM.type, ...))
  } else {
    w.mat <- a.mat
  }

  # assign row and col names
  rownames(w.mat) <- rownames(a.mat)
  colnames(w.mat) <- colnames(a.mat)

  # dissimilarity measure
  warning("\nComputing dissimilarity matix...")
  d.mat <- 1- w.mat

  output <- list(
    s.mat = s.mat,
    a.mat = a.mat,
    w.mat = w.mat,
    d.mat = d.mat
  )

  return(output)
}



#' Hierarchially-cluster distance matrix
#'
#' Hierarchially-cluster distance matrix using flashClust package.
#'
#' @param d.mat distance matrix
#' @param method hierarchial-clustering method. Default is "average". See flashClust for options.
#' @param ... Additional arguments passed to flashClust {flashClust package}
#' @name dist2hclust
#' @seealso \code{\link{flashClust}}
#' @return hclust object
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # Get expression matrix
#' which.data <- "scale"
#'
#' # variable gene only matrix
#' use.var <- T
#' if (use.var){
#'   exp.mat <- getExpressionMatrix(so.query, only.variable = use.var, which.data = which.data, use.additional.genes = NA)
#' } else {
#'   exp.mat <- exp.mat.complete
#' }
#'
#'
#' # transpose expressio matrix (genes are columns)
#' t.exp.mat <- t(exp.mat)
#' datExpr <- as.matrix(t.exp.mat)
#' SubGeneNames=colnames(datExpr)
#'
#' # capture output used to hide undesired print statement
#' print2hide <- capture.output(allowWGCNAThreads())
#'
#' # transform matrix if necessary
#' if (min(datExpr) < 0) {
#'   datExpr.noz <- datExpr + abs(min(datExpr))
#' } else {
#'   datExpr.noz <- datExpr
#' }
#'
#' # run WGCNA
#' output.all <- runWGCNA(datExpr.noz, cor.metric = "rho_p", soft.power = 2, use.TOM = T)
#'
#' # unpack output
#' s.mat <- output.all[["s.mat"]] # similar matrix
#' a.mat <- output.all[["a.mat"]] # adjacency matrix
#' w.mat <- output.all[["w.mat"]] # topological overlap matix
#' d.mat <- output.all[["d.mat"]] # disimilarity matix
#'
#' geneTree <- dist2hclust(d.mat)
#'
dist2hclust <- function(d.mat, method = "average", ...){

  stopifnot(any(class(d.mat) == "matrix"))
  library(flashClust)

  tree <- flashClust(as.dist(d.mat), method = method, ...)

  return(tree)
}


#' Cut clustered tree at varying heights and overlay with dendrogram to find optimal parameter set.
#'
#' WGCNA::cutreeDynamic is run for varying levels of deepSplit parameter (0:4), and cluster membership is assigned for each parameter set.
#'
#' @param tree h.clust object generated by dist2hclust.
#' @param d.mat distance matrix used to generate tree.
#' @param genes vector of gene names corresponding rows/col of distance matrix (d.mat). If specified, additional "genes" column is provided in output.
#' @param ... additional arguments passed to dynamicTreeCut::cutreeDynamic()
#' @name optimalDS
#' @seealso \code{\link{cutreeDynamic}}
#' @return mColorh, matrix specfying module membership at varying deep split parameter specifications (0:4)
#' @examples
#'
#' geneTree <- dist2hclust(d.mat)
#' # determine number of modules based on refrence dataset
#' print2hide <- capture.output(mColorh <- optimalDS(tree = geneTree, d.mat = d.mat, genes  = rownames(a.mat)))
#'
optimalDS <- function(tree, d.mat, genes = NULL, ...){

  library(dynamicTreeCut)

  mColorh = NULL
  for (ds in 0:4){
    cut.tree <- cutreeDynamic(dendro = tree,distM= d.mat, cutHeight = 0.998, deepSplit=ds, pamRespectsDendro = FALSE, ...)
    mColorh <- cbind(mColorh, labels2colors(cut.tree))
  }

  colnames(mColorh) <- paste0("ds.", seq(0,4))

  if (!is.null(genes) & length(genes) == nrow(mColorh)){
    mColorh <- cbind(mColorh, genes)
    colnames(mColorh)[length(colnames(mColorh))] <- "genes"
  }

  return(mColorh)
}




#' Determine module preservation between reference and query network
#'
#' Determine module preservation between reference and query network using expression matrices from two scRNAseq comparison groups. Uses WGCNA::modulePreservation() to assess how well a module in one sample is preserved in another. 5<Z<10 indicates moderate presevation, while Z>10 indicates high preservation. Grey module contains uncharacterized genes while gold module contains random genes (these are used as controls). Note that future updates will extend the functionality to accomodate >2 networks.
#'
#' @param ref.mat reference data (expression matrix, cols are genes and rows and samples)
#' @param query.mat query data (expression matrix, cols are genes and rows and samples)
#' @param ref.modules reference module membership. Vector of colors (length is equal to number of samples), specfiying sample membership to each module.
#' @param query.modules query module membership. Vector of colors (length is equal to number of samples), specfiying sample membership to each module.
#' @param networkType Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". See WGCNA::adjacency.
#' @param referenceNetworks a vector giving the indices of expression data to be used as reference networks. Reference networks must have their module labels given in multiColor.
#' @param ... Additional arguments passessed to modulePreservation {WGCNA package}
#' @name getModulePreservation
#' @seealso \code{\link{modulePreservation}}
#' @return data.frame of module preservation statistics
#' @import WGCNA
#' @examples
#'
#' mColorh.1 <- optimalDS(dist2hclust(d.1), d.1, pamStage = F,cutHeight = 0.99)
#' modules.1 = mColorh.1[,4]
#'
#' stats <- getModulePreservation(ref.mat = de.1, query.mat = de.2, modules.1, nPermutations=30,maxGoldModuleSize=100,maxModuleSize=400, verbose=3)
#'
getModulePreservation <- function(ref.mat, query.mat, ref.module, query.modules = NULL, networkType = "unsigned", referenceNetworks = 1, ...){

  library(WGCNA)
  multiExpr <- list(A1=list(data=ref.mat),A2=list(data=query.mat))

  if (is.null(query.modules)){
    multiColor <- list(A1 = ref.module)
  } else {
    multiColor <- list(A1 = ref.module, A2 = query.modules)
  }

  mp <- modulePreservation(multiExpr,multiColor,referenceNetworks=referenceNetworks,networkType=networkType, ...)
  stats <- mp$preservation$Z$ref.A1$inColumnsAlsoPresentIn.A2


  return(stats)
}


#' Analysis of scale free topology for soft-threshold
#'
#' Analysis of scale free topology for multiple soft thresholding powers. The aim is to help the user pick an appropriate soft-thresholding power for network construction. This is an adaptation of the WGCNA::pickSoftThreshold function which has been customized for scRNAseq applications.
#'
#' @param s.mat similarity matrix
#' @param RsquaredCut Rsq cutoff. Default is 0.85.
#' @param networkType Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". See WGCNA::adjacency.
#' @param ... Additional arguments passessed to pickSoftThreshold {pickSoftThreshold}
#' @name getSoftThreshold
#' @return list of soft threshold picks
#' @import doParallel
#' @examples
#'
#' # determine optimal soft threshold
#' sft <- getSoftThreshold(s.mat)
#'
#' # Plot the results
#' sizeGrWindow(9, 5)
#' par(mfrow = c(1,2));
#' cex1 = 0.9;
#'
#' # Scale-free topology fit index as a function of the soft-thresholding power
#' plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
#' text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
#'
#' # Red line corresponds to using an R^2 cut-off
#' abline(h=0.80,col="red")
#'
#' # Mean connectivity as a function of the soft-thresholding power
#' plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
#' text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#'
getSoftThreshold <- function (s.mat, dataIsExpr = F, weights = NULL, RsquaredCut = 0.85,
                              powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)),
                              removeFirst = FALSE, nBreaks = 10, blockSize = 1000, corFnc = cor,
                              corOptions = list(use = "p"), networkType = "signed",
                              moreNetworkConcepts = FALSE, gcInterval = NULL, verbose = 0,
                              indent = 0)
{
  data <- s.mat
  powerVector = sort(powerVector)
  networkTypes <- c("unsigned", "signed", "signed hybrid")
  intType = charmatch(networkType, networkTypes)
  if (is.na(intType))
    stop(paste("Unrecognized 'networkType'. Recognized values are",
               paste(networkTypes, collapse = ", ")))
  nGenes = ncol(data)
  if (nGenes < 3) {
    stop("The input data data contain fewer than 3 rows (nodes).",
         "\nThis would result in a trivial correlation network.")
  }
  if (!dataIsExpr) {
    checkSimilarity(data)
    if (any(diag(data) != 1))
      diag(data) = 1
  }
  if (is.null(blockSize)) {
    blockSize = blockSize(nGenes, rectangularBlocks = TRUE,
                          maxMemoryAllocation = 2^30)
    if (verbose > 0)
      printFlush(spaste("pickSoftThreshold: will use block size ",
                        blockSize, "."))
  }
  if (length(gcInterval) == 0)
    gcInterval = 4 * blockSize
  colname1 = c("Power", "SFT.R.sq", "slope", "truncated R.sq",
               "mean(k)", "median(k)", "max(k)")
  if (moreNetworkConcepts) {
    colname1 = c(colname1, "Density", "Centralization",
                 "Heterogeneity")
  }
  datout = data.frame(matrix(666, nrow = length(powerVector),
                             ncol = length(colname1)))
  names(datout) = colname1
  datout[, 1] = powerVector
  spaces = indentSpaces(indent)
  if (verbose > 0) {
    cat(paste(spaces, "pickSoftThreshold: calculating connectivity for given powers..."))
    if (verbose == 1)
      pind = initProgInd()
    else cat("\n")
  }
  corFnc = match.fun(corFnc)
  corFormals = formals(corFnc)
  if ("nThreads" %in% names(corFormals))
    corOptions$nThreads = 1
  datk = matrix(0, nrow = nGenes, ncol = length(powerVector))
  nThreads = WGCNAnThreads()
  nPowers = length(powerVector)
  startG = 1
  lastGC = 0
  corOptions$x = data
  if (!is.null(weights)) {
    if (!dataIsExpr)
      stop("Weights can only be used when 'data' represents expression data ('dataIsExpr' must be TRUE).")
    if (!isTRUE(all.equal(dim(data), dim(weights))))
      stop("When 'weights' are given, dimensions of 'data' and 'weights' must be the same.")
    corOptions$weights.x = weights
  }
  while (startG <= nGenes) {
    endG = min(startG + blockSize - 1, nGenes)
    if (verbose > 1)
      printFlush(paste(spaces, "  ..working on genes",
                       startG, "through", endG, "of", nGenes))
    nBlockGenes = endG - startG + 1
    jobs = allocateJobs(nBlockGenes, nThreads)
    actualThreads = which(sapply(jobs, length) > 0)
    datk[c(startG:endG), ] = foreach(t = actualThreads,
                                     .combine = rbind) %dopar% {
                                       useGenes = c(startG:endG)[jobs[[t]]]
                                       nGenes1 = length(useGenes)
                                       if (dataIsExpr) {
                                         corOptions$y = data[, useGenes]
                                         if (!is.null(weights))
                                           corOptions$weights.y = weights[, useGenes]
                                         corx = do.call(corFnc, corOptions)
                                         if (intType == 1) {
                                           corx = abs(corx)
                                         }
                                         else if (intType == 2) {
                                           corx = (1 + corx)/2
                                         }
                                         else if (intType == 3) {
                                           corx[corx < 0] = 0
                                         }
                                         if (sum(is.na(corx)) != 0)
                                           warning(paste("Some correlations are NA in block",
                                                         startG, ":", endG, "."))
                                       }
                                       else {
                                         corx = data[, useGenes]
                                       }
                                       ind = cbind(useGenes, 1:length(useGenes))
                                       corx[ind] = 1
                                       datk.local = matrix(NA, nGenes1, nPowers)
                                       corxPrev = matrix(1, nrow = nrow(corx), ncol = ncol(corx))
                                       powerVector1 <- c(0, head(powerVector, -1))
                                       powerSteps <- powerVector - powerVector1
                                       uniquePowerSteps <- unique(powerSteps)
                                       corxPowers <- lapply(uniquePowerSteps, function(p) corx^p)
                                       names(corxPowers) <- uniquePowerSteps
                                       for (j in 1:nPowers) {
                                         corxCur <- corxPrev * corxPowers[[as.character(powerSteps[j])]]
                                         datk.local[, j] = colSums(corxCur, na.rm = TRUE) -
                                           1
                                         corxPrev <- corxCur
                                       }
                                       datk.local
                                     }
    startG = endG + 1
    if ((gcInterval > 0) && (startG - lastGC > gcInterval)) {
      gc()
      lastGC = startG
    }
    if (verbose == 1)
      pind = updateProgInd(endG/nGenes, pind)
  }
  if (verbose == 1)
    printFlush("")
  for (i in c(1:length(powerVector))) {
    khelp = datk[, i]
    if (any(khelp < 0))
      # browser()
      SFT1 = scaleFreeFitIndex(k = khelp, nBreaks = nBreaks,
                               removeFirst = removeFirst)
    datout[i, 2] = SFT1$Rsquared.SFT
    datout[i, 3] = SFT1$slope.SFT
    datout[i, 4] = SFT1$truncatedExponentialAdjRsquared
    datout[i, 5] = mean(khelp, na.rm = TRUE)
    datout[i, 6] = median(khelp, na.rm = TRUE)
    datout[i, 7] = max(khelp, na.rm = TRUE)
    if (moreNetworkConcepts) {
      Density = sum(khelp)/(nGenes * (nGenes - 1))
      datout[i, 8] = Density
      Centralization = nGenes * (max(khelp) - mean(khelp))/((nGenes -
                                                               1) * (nGenes - 2))
      datout[i, 9] = Centralization
      Heterogeneity = sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 -
                             1)
      datout[i, 10] = Heterogeneity
    }
  }
  # print(signif(data.frame(datout), 3))
  ind1 = datout[, 2] > RsquaredCut
  indcut = NA
  indcut = if (sum(ind1) > 0)
    min(c(1:length(ind1))[ind1])
  else indcut
  powerEstimate = powerVector[indcut][[1]]
  gc()

  output <- list(powerEstimate = powerEstimate, fitIndices = data.frame(datout))
  return(output)
}



#' Wrapper to run  GO enrichment, using fisher enrichment method
#'
#' Wrapper to run  GO enrichment, using fisher enrichment method. Builds on TopGo functionality.
#'
#' @param gene.list named list of gene sets to query, where name specify name of gene set, and entry is vector of gene symbols.
#' @param gene.universe background genes
#' @param species Species: Hs or Mm
#' @param p.threshold numeric specifying p value threshold. Default is 0.01.
#' @param p.adj.threshold numeric specifing adjusted p value threshold. Defulat is 0.05.
#' @name runEnrichment
#' @seealso \code{\link{enrichGO.fisher}}
#' @return list of data.frames. results.table.p contains unadjusted results, results.table.bh contains BH-adjusted results.
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' enrichment.list <- runEnrichment(module.list, gene.universe = all.genes, species = "Hs")
#' results.p <- enrichment.list$results.table.p
#' results.bh <- enrichment.list$results.table.bh
#'
runEnrichment <- function(gene.list, gene.universe, species, p.threshold = 0.01, p.adj.threshold = 0.05){

  library(topGO)

  #get list of significant GO before multiple testing correction
  results.table.p <- NULL

  #get list of significant GO after multiple testing correction
  results.table.bh <- NULL

  # enrich each module
  topGO.data <- NULL
  for (i in 1:length(gene.list)){

    if (species == "Hs"){

      enrich.list <- enrichGO.fisher(gene.candidates = toupper(gene.list[[i]]),
                                     gene.universe = toupper(gene.universe),
                                     which.species = species,
                                     p.threshold = p.threshold,
                                     padj.threshold = p.adj.threshold,
                                     topGO.object = topGO.data)

    } else if (species == "Mm"){

      enrich.list <- enrichGO.fisher(gene.candidates = gene.list[[i]],
                                     gene.universe = gene.universe,
                                     which.species = species,
                                     p.threshold = p.threshold,
                                     padj.threshold = p.adj.threshold,
                                     topGO.object = topGO.data)

    }

    # get results
    results.table.p.cur <- enrich.list$unadjusted.results
    results.table.bh.cur <- enrich.list$adjusted.results

    # get and reuse topGo object by updating gene list.
    topGO.data <- enrich.list$topGo.object


    if (nrow(results.table.p.cur) > 0) {
      results.table.p.cur$module <- names(gene.list)[i]
      results.table.p <- bind_rows(results.table.p, results.table.p.cur)
    }

    if (nrow(results.table.bh.cur) > 0) {
      results.table.bh.cur$module <-  names(gene.list)[i]
      results.table.bh <- bind_rows(results.table.bh, results.table.bh.cur)
    }
  }

  output <- list(
    results.table.p = results.table.p,
    results.table.bh = results.table.bh
  )

  return(output)
}


#' Get list of module genes
#'
#' Map module membership to gene list and return module gene list.
#'
#' @param module.colors Vector of module membership (usually colors). Must be same length as gene list.
#' @param genes Vector of genes.
#' @param add.prefix Logical specifying whether to add "M#." prefix to each module name.
#' @name getModuleGenes
#' @return Named list, where names are module names and entries are module genes.
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' module.list.1 <- getModuleGenes(modules.1, SubGeneNames, add.prefix = F)
#' module.list.2 <- getModuleGenes(modules.2_new, SubGeneNames, add.prefix = F)
#'
#' # get genes that are common to both network modules
#' module.union.list <- list()
#' common.modules <- intersect(modules.1, modules.2_new)
#' for (i in 1:length(common.modules)){
#'   common.genes <- intersect(module.list.1[[common.modules[i]]], module.list.2[[common.modules[i]]])
#'   module.union.list[[common.modules[i]]] <- common.genes
#' }
#'
#'
getModuleGenes <- function(module.colors, genes, add.prefix = T){

  module_colors= setdiff(unique(module.colors), "grey")
  module.list <- list()
  for (i in 1:length(module_colors) ){
    color <- module_colors[i]

    if (add.prefix){
      module.name <- paste("M", i, ".", color, sep = "")
    } else {
      module.name <- color
    }

    module.list[[module.name]]=genes[which(module.colors==color)]
  }

  return(module.list)
}


#' Convert adjacency or TOM matrix to igraph data.frame
#'
#' Convert TOM matrix (from WGCNA analysis) to igraph data.frame and filter top connections. Uses igraph::graph_from_adjacency_matrix to convert matrix to data.frame.
#'
#' @param w.mat TOM or adjacency matrix.
#' @param top.n Numeric specifying N top connections to return (N < number of connections). If top.n and top.percentile are NULL, all connections retained.
#' @param top.percentile Numeric [0,1] specifying Nth percetile of top connections to return. If top.n and top.percentile are NULL, all connections retained.
#' @param graph.type Type of graph. Default is "undirected"
#' @name wgcna2graphDF
#' @return igraph data.frame
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # get connectivity for specified module
#' module.name <- names(module.list.all)[names(module.list.all) %in% which.modules]
#' module.gene.cur <- module.list.all[[module.name]]
#' w.cur <- w.mat[rownames(w.mat) %in% module.gene.cur, colnames(w.mat) %in% module.gene.cur]
#'
#' # get igraph data.frame for subset of connections
#' w.df.top <- wgcna2graphDF(w.cur, top.n = top.n.interactions)
#' w.df.top$module.membership <- module.name
#'
#'
wgcna2graphDF <- function(w.mat, top.n = NULL, top.percentile = NULL, graph.type = "undirected"){

  # matrix to graph
  w.graph <- igraph::graph_from_adjacency_matrix(w.mat, mode = c( graph.type), diag = F, weighted = T)

  # graph to data.frame
  w.df <- igraph::as_data_frame(w.graph)

  # return top connections
  if (is.null(top.n) & is.null(top.percentile)){
    return(w.df)
  } else if (!is.null(top.n)) {
    w.df.top <- w.df[rank(w.df$weight) > (nrow(w.df) - top.n), ]
    return(w.df.top)
  } else if (!is.null(top.percentile)){
    w.df.top <- w.df[rank(w.df$weight)/nrow(w.df) > (1-top.percentile), ]
    return(w.df.top)
  }

}



#' Get nodes and edges from igraph data.frame for visNetwork
#'
#' Convert igraph data.frame to dataframe of nodes and edges (used as input to visNetwork). Only applicable for undirected networks.
#'
#' @param df.data igraph data.frame (output from scMiko::wgcna2graphDF)
#' @name getNodesEdges
#' @return named list of nodes and edges, where each entry is a data.frame.
#' @examples
#'
#' # get connectivity for specified module
#' module.name <- names(module.list.all)[names(module.list.all) %in% which.modules]
#' module.gene.cur <- module.list.all[[module.name]]
#' w.cur <- w.mat[rownames(w.mat) %in% module.gene.cur, colnames(w.mat) %in% module.gene.cur]
#'
#' # get igraph data.frame for subset of connections
#' w.df.top <- wgcna2graphDF(w.cur, top.n = top.n.interactions)
#' w.df.top$module.membership <- module.name
#'
#' # get nodes and edges
#' node.edge.output <- getNodesEdges(w.df.top)
#' vis.nodes <- node.edge.output$nodes
#' vis.links <- node.edge.output$edges
#'
getNodesEdges <- function(df.data){

  gD.cur <- igraph::simplify(igraph::graph.data.frame(df.data, directed=FALSE))

  # convert to VisNet representation
  visNet.data <- toVisNetworkData(gD.cur)

  # get node and link components
  vis.nodes <- visNet.data$nodes
  vis.links <- visNet.data$edges

  output <- list(nodes = vis.nodes, edges = vis.links)

  return(output)
}


#' Convert values to color gradient
#'
#' Converts ranges of values to corresponding color along color gradient, specified by 3 colors (low color, middle color and high color)
#'
#' @param Values vector of numerical values to convert to colors
#' @param limit numeric specifying limit of color range. If unspecified, limit <- max(c(abs(min(Values)), abs(max(Values))))
#' @param gradient.length Numeric specifying number of bins to split gradient into. Default is 100.
#' @param low.col Color representing low values. Default is "skyblue"
#' @param mid.col Color representing mid values. Default is "grey"
#' @param high.col Color representing high values. Default is "tomato"
#' @name value2col
#' @return vector of colors.
#' @seealso \code{\link{colorRampPalette}}
#' @examples
#'
#' # get edge colors
#' col.lim <- max(c(abs(min(vis.links$value)), abs(max(vis.links$value))))
#' edge.colors <- value2col(vis.links$value, limit = col.lim)
#'
value2col <- function(values, limit = NULL, gradient.length = 100, low.col = "skyblue", mid.col = "grey", high.col = "tomato"){

  # truncate values if necessary
  values[values > limit] <- limit
  values[values < -limit] <- -limit

  jj <- cut(values, breaks = seq(-limit, limit, len = gradient.length),
            include.lowest = TRUE)
  cols <- colorRampPalette(c(low.col, mid.col, high.col))(gradient.length-1)[jj]

  return(cols)
}


#' Convert named list to long data.frame
#'
#' Convert named list to long data.frame. Resulting dataframe will have two columns, the first corresponding to the names within the list, and the second to the corersponding list entries.
#'
#' @param my.list named list
#' @param name.header character specifying header name that will be assigned to name column. If unspecified, defaults to "name".
#' @param value.header character specifying header name that will be assigned to value column. If unspecified, defaults to "value".
#' @name namedList2longDF
#' @return long data.frame
#' @seealso \code{\link{namedList2wideDF}}
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # specify named list
#' my.list <- list(group.1 = c("a", "b", "c"), group.2 = c("d", "e", "f"))
#'
#' # convert named list to long data.frame
#' my.df <- namedList2longDF(my.list)
#'
namedList2longDF <- function(my.list, name.header = NULL, value.header = NULL){

  if (class(my.list) != "list") stop("input is not a list")
  my.list <- my.list[lapply(my.list, length) > 0]

  # unlist and assign values to data.frame
  my.df <- NULL
  for (i in 1:length(my.list)){
    my.df <- bind_rows(my.df, data.frame(name = names(my.list)[i], value = as.vector(unlist(my.list[i]))))
  }

  # get data.frame column names
  if (is.null(value.header)) value.header <- "value"
  if (is.null(name.header)) name.header <- "name"
  colnames(my.df) <- c(name.header, value.header)

  # return long data.frame
  return(my.df)

}


#' Convert wide data.frame to named list
#'
#' Convert wide data.frame to named list. Column entries are used as list entries, and each entry is named using the corresponding column name in the data.frame.
#'
#' @param df.wide wide data.frame
#' @name wideDF2namedList
#' @return named list
#' @seealso \code{\link{namedList2wideDF}}
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # get wide data frame
#' CancerSEA_Hs <-geneSets[["CancerSEA_Hs"]]
#'
#' # convert wide data frame to named list
#' my.list <- namedList2longDF(CancerSEA_Hs)
#'
wideDF2namedList <- function(df.wide){
  try({df.wide <- as.data.frame(df.wide)}, silent = T)
  if (!("data.frame" %in% class(df.wide))) stop("Input must be a data frame")

  n.list <- list()
  for (i in 1:ncol(df.wide)){

    col.name <- colnames(df.wide)[i]
    entries <- df.wide[ ,i]
    entries <- entries[!is.na(entries)]
    entries <- entries[entries != ""]
    n.list[[col.name]] <- entries

  }

  return(n.list)
}


#' Convert named list to wide data.frame
#'
#' Convert named list to wide data.frame. Resulting dataframe will have the same number of columns as there are names in the list.
#'
#' @param my.list named list
#' @name namedList2wideDF
#' @return wide data.frame
#' @seealso \code{\link{namedList2longDF}}
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # specify named list
#' my.list <- list(group.1 = c("a", "b", "c"), group.2 = c("d", "e", "f"))
#'
#' # convert named list to wide data.frame
#' my.df <- namedList2wideDF(my.list)
#'
namedList2wideDF <- function(my.list){

  if (class(my.list) != "list") stop("input is not a list")

  df.long <- scMiko::namedList2longDF(my.list, name.header = "name", value.header = "value")
  df.long$name <- as.character(df.long$name)
  df.long$value <- as.character(df.long$value)
  df.wide <- df.long %>%
    dplyr::group_by(name) %>%
    dplyr::mutate(row = row_number()) %>%
    tidyr::pivot_wider(names_from = name, values_from = value) %>%
    dplyr::select(-row)

  # return wide data.frame
  return(df.wide)
}



#' Compute adjaceny matrix from similary (correlation) matrix
#'
#' Compute adjaceny matrix from similary (correlation) matrix for different network types.
#'
#' @param s.mat Similarity matrix. Symmetric matrix with rows and columns corresponding to genes, and diagonal = 1.
#' @param soft.power Numeric specifying soft power used to calculate adjacency matrix.
#' @param network.type character specifying network type. Must be one of "unsigned", "signed", or "signed hybrid".
#' @name sim2adj
#' @author Nicholas Mikolajewicz
#' @return adjaceny matrix (a.mat)
#' @examples
#'
#' # Convert similary to adjacency
#' a.mat <- sim2adj(s.mat, soft.power = 2, network.type = "signed")
#'
sim2adj <- function(s.mat, soft.power, network.type){

  # adjacency matrix
  softPower <- soft.power
  if ( network.type == "unsigned"){
    a.mat <- abs(s.mat)^softPower
  } else if ( network.type == "signed"){
    a.mat <-  (0.5 * (1+s.mat) )^softPower
  } else if ( network.type == "signed hybrid"){
    a.mat <- s.mat
    a.mat[a.mat <= 0] <- 0
    a.mat <- (a.mat)^softPower
  } else {
    stop("Network type incorrectly specified. Must be one of 'signed', 'unsighed', or 'signed hybrid'")
  }

  return(a.mat)

}



#' Rescale values to specified range.
#'
#' Rescales vector of values to span between specified range [new.min, new.max]
#'
#' @param values Numeric vector of values to rescale.
#' @param new.min Numeric specfiying new minimum. If unspecified, default is 0.
#' @param new.max Numeric specfiying new maximum. If unspecified, default is 1.
#' @name rescaleValues
#' @author Nicholas Mikolajewicz
#' @return Numeric vector of rescaled values
#' @examples
#'
#' # rescale values
#' values <- rescaleValues(values)
#'
rescaleValues <- function(values, new.min = 0, new.max = 1){

  # set lower bound to zero
  old.min <- min(values)
  if (old.min < 0) {
    values <- values + abs(old.min)
  } else if (old.min > 0) {
    values <- values - abs(old.min)
  }
  stopifnot( min(values) == 0)

  # set upper bound to one
  old.max <- (max(values))
  values <- values/old.max
  stopifnot( max(values) == 1)

  new.range <- new.max - new.min
  values <- values * new.range
  values <- values + new.min

  stopifnot(min(values) == new.min)
  stopifnot(max(values) == new.max)

  return(values)

}



#' Bayesian Correlation algorithm
#'
#' Bayesian correlation scheme that assigns low similarity to genes that have low confidence expression estimates. Shown to be more reproducible than Pearson correlations. Source: https://www.biorxiv.org/content/10.1101/714824v1
#'
#' @param X Expression matrix
#' @name BaCo
#' @return Correlation matrix
#' @examples
#'
#' #create a matrix (or load your own)
#' X <- matrix(1:1000, ncol=20)
#'
#' #compute the Bayesian correlation matrix
#' B <- BaCo(X)
#'
BaCo <- function(X){

  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE)
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
}




#' Analysis of scale free topology for soft-threshold. Modified from getSoftThreshold.
#'
#' Analysis of scale free topology for multiple soft thresholding powers. The aim is to help the user pick an appropriate soft-thresholding power for network construction. Inspired by WGCNA::pickSoftThreshold and updated from first version (scMiko::getSoftThreshold)
#'
#' @param s.mat similarity matrix
#' @param power Numeric vector of powers to evaluate. Default is c(seq(0.5,5, by = 0.5), seq(6,10))
#' @param networkType Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". See WGCNA::adjacency.
#' @param nBreaks Number of bins in connectivity histograms. Default is 20.
#' @param removeFirst Logical specifying whether the first bin should be removed from the connectivity histogram. Default is True.
#' @param rescale.adjacency Logical indicating if s.mat should be rescaled to [0,1]
#' @param n.cores Number of cores to use for parallelization. Default is 4.
#' @name getSoftThreshold2
#' @seealso \code{\link{getConnectivity}}
#' @return named list containing power estimates, r2 estimates, distribution plots, optimiation plot and results data.frame.
#' @examples
#'
#' # determine optimal soft threshold
#' sft <- getSoftThreshold2(s.mat, power =c(seq(0.5,5, by = 0.5), seq(6,10)),
#' network.type = "signed", rescale.adjacency = F)
#'
#' # visualize optimization plot
#' print(sft$optimization.plot)
#'
#' # visualize node-linkage density plots
#' cowplot::plot_grid(plotlist = sft$distribution.plot, ncol = 5)

getSoftThreshold2 <- function(s.mat, power =c(seq(0.5,5, by = 0.5), seq(6,10)), network.type = "signed", nBreaks = 20, removeFirst = T, rescale.adjacency = F, n.cores = 4){


  plt.sf.list <- list()
  powers <- power
  r2.sf <- c()

  # start cluster
  cl <- makeCluster(n.cores)
  registerDoParallel(cl)

  sf.list <- list()

  sf.list <- (foreach(i = 1:length(powers), .packages = c("scMiko", "dplyr", "ggplot2", "ggpmisc"))) %dopar% {

    power.cur <- powers[i]
    a.cur <-  sim2adj(s.mat, soft.power = power.cur, network.type = network.type)

    if (rescale.adjacency)  a.cur <- recaleValues(a.cur, new.min = 0, new.max = 1)

    # get connectivity for specified power
    net.connectivity.df <- getConnectivity(a.cur, rownames(a.cur), flag.top.n = 20)

    # nBreaks <- 20
    # removeFirst <- T
    k <- net.connectivity.df$wi
    discretized.k <- cut(k, nBreaks)
    dk <- tapply(k, discretized.k, mean)
    p.dk <- as.vector(tapply(k, discretized.k, length)/length(k))
    breaks1 <- seq(from = min(k), to = max(k), length = nBreaks + 1)
    hist1 <- suppressWarnings(hist(k, breaks = breaks1, equidist = FALSE,
                                   plot = FALSE, right = TRUE))
    dk2 <- hist1$mids
    dk <- ifelse(is.na(dk), dk2, dk)
    dk = ifelse(dk == 0, dk2, dk)
    p.dk <- ifelse(is.na(p.dk), 0, p.dk)
    log.dk <- as.vector(log10(dk))
    if (removeFirst) {
      p.dk = p.dk[-1]
      log.dk = log.dk[-1]
    }
    log.p.dk <- as.numeric(log10(p.dk + 1e-09))
    lm1 <- lm(log.p.dk ~ log.dk)

    r2.sf[i] <- summary(lm1)[["r.squared"]] * sign(lm1[["coefficients"]][["log.dk"]])

    df.sf <- data.frame(x = log.dk, y = log.p.dk)

    # store node linkage distribution plot
    plt.sf.list[[as.character(power.cur)]] <- df.sf %>%
      ggplot(aes(x=x, y=y)) +
      geom_smooth(method = "lm", color = "tomato", fill = "tomato") +
      geom_point(size = 3) +
      xlab("N Links (Log)") +
      ylab("N Nodes (Log)") +
      ggtitle(paste0("Node Linkages\nSoft Power = ", power.cur)) +
      theme_classic() +
      stat_fit_glance(method = "lm",
                      label.y = "bottom",
                      method.args = list(formula = y ~ x),
                      mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',
                                                    stat(r.squared), stat(p.value))),
                      parse = TRUE)


    sf.list[[i]] <- list(
      r2 = r2.sf[i],
      df = df.sf,
      plt = plt.sf.list[[as.character(power.cur)]]
    )

    return(sf.list)

  }

  stopCluster(cl)

  # unpack results
  r2.sf <- c(); plt.sf.list <- list();
  for (i in 1:length(sf.list)){
    power.cur <- powers[i]
    r2.sf[i] <- sf.list[[i]][[i]][["r2"]]
    plt.sf.list[[as.character(power.cur)]] <- sf.list[[i]][[i]][["plt"]]
  }

  # store powers and r2
  df.r2.sf <- data.frame(sf = powers, r2 = r2.sf)

  # optimizatio plot
  plt.opt.sf <- df.r2.sf %>%
    ggplot(aes(x = sf, y = r2)) +
    geom_hline(yintercept = df.r2.sf$r2[which.min(df.r2.sf$r2)], color = "tomato") +
    geom_vline(xintercept = df.r2.sf$sf[which.min(df.r2.sf$r2)], color = "tomato") +
    geom_smooth(method = "loess", color = "black", fill = "grey") +
    geom_point(size = 3) +
    xlab("Soft Power") +
    ylab("R2 (Scale Free Topology)") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(title = "Soft Power Optimization",
         subtitle = paste0("Optimal power = ", df.r2.sf$sf[which.min(df.r2.sf$r2)], ", r2 = ", signif(df.r2.sf$r2[which.min(df.r2.sf$r2)], 3)))

  # store results
  output <- list(
    powerEstimate = df.r2.sf$sf[which.min(df.r2.sf$r2)],
    r2Estimate = df.r2.sf$r2[which.min(df.r2.sf$r2)],
    distribution.plot = plt.sf.list,
    optimization.plot = plt.opt.sf,
    results = df.r2.sf
  )

  return(output)


}



#' Balance matrix dimensions
#'
#' Resamples two matrices to match number of rows in each. Number of rows is matched to minimum or maximum, as specified.
#'
#' @param s.mat similarity matrix
#' @param mat.1 First input matrix.
#' @param mat.2 Second input matrix.
#' @param method Character specifying which matching method to use. Must be one of "match.max" (default) or "match.min". For either option, matrix with min/max number of rows is resamples to match the other matrix.
#' @name balanceMatrixSize
#' @author Nicholas Mikolajewicz
#' @return list of resampled matrices.
#' @examples
#'
#' # ensure sample sizes are balanced across groups
#' output.mat <- balanceMatrixSize(de.orig.1, de.orig.2, method = "match.max")
#' de.1 <- output.mat[["de.1"]] # resampled mat.1
#' de.2 <- output.mat[["de.2"]] # resampled mat.2
#' datExpr.noz <- output.mat[["de.all"]] # concatenated matrices
#'
balanceMatrixSize <- function(mat.1, mat.2, method = "match.max"){

  # ensure sample sizes are balanced across groups
  match.sample.size <- method #option: "match.max", match to max sample size; "match.min", match to min sample size; "none", no matcing

  size.1 <- nrow(mat.1)
  size.2 <- nrow(mat.2)

  if (match.sample.size == "match.max"){

    target.size <- max(c(size.1, size.2))

    fill.1 <- target.size - size.1
    fill.2 <- target.size - size.2

    sample.ind.1 <- sample(seq(1,size.1), fill.1, replace = T)
    sample.ind.2 <- sample(seq(1,size.2), fill.2, replace = T)

    if (fill.1 > 0){
      de.1 <- rbind(mat.1, mat.1[sample.ind.1 , ])
    } else {
      de.1 <- mat.1
    }

    if (fill.2 > 0){
      de.2 <- rbind(mat.2, mat.2[sample.ind.2 , ])
    } else {
      de.2 <- mat.2
    }

    de.all <- rbind(de.1, de.2)


  } else if (match.sample.size == "match.min"){

    target.size <- min(c(size.1, size.2))

    fill.1 <- size.1 - target.size
    fill.2 <- size.2 - target.size

    sample.ind.1 <- sample(seq(1,size.1), target.size, replace = F)
    sample.ind.2 <- sample(seq(1,size.2), target.size, replace = F)

    de.1 <- mat.1[sample.ind.1, ]
    de.2 <- mat.2[sample.ind.2, ]

    de.all <- rbind(de.1, de.2)


  } else if (match.sample.size == "none"){

    de.1 <- mat.1
    de.2 <- mat.2

    de.all <- rbind(de.1, de.2)
  }


  output <- list(
    de.1 = de.1,
    de.2 = de.2,
    de.all = de.all
  )

  return(output)

}



#' TOM matrix rescaling prior to computing consensus topological overlap.
#'
#' Query TOM matrix is rescaled to match nth percentile of reference TOM matrix. Since consensus is defined as component-wise minimum of two TOMS, a bias may result without scaling. This is mitigated by scaling the matrices.
#'
#' @param query.TOM Query topological overlap matrix. This is the matrix that will be scaled. Corresponds to w.mat output from scMiko::runWGCNA() function.
#' @param reference.TOM Reference topological overlap matrix. This is the matrix that will be used as the reference and will not be scaled. Corresponds to w.mat output from scMiko::runWGCNA() function.
#' @param reference.percentile Numerical. Matrix is scaled so that reference.percentile is matched across both matrices. Default is 0.95.
#' @name scaleTOM
#' @return scaled query.TOM
#' @examples
#'
#' # scale TOM matrices
#' w.1 <- scaleTOM(query.TOM = w.1, reference.TOM = w.mat, reference.percentile = 0.95)
#' w.2 <- scaleTOM(query.TOM = w.2, reference.TOM = w.mat, reference.percentile = 0.95)
#'
#' # recompute the distance matrix
#' d.mat <- 1-w.mat
#' d.1 <- 1-w.1
#' d.2 <- 1-w.2
#'
#' # compute consensus topological overlap
#' w.co <- pmin(w.1, w.2) # component-wise parallel minimum of TOMs
#' d.co <- 1-w.co
#'
scaleTOM <- function(query.TOM, reference.TOM, reference.percentile = 0.95){

  # TOM matrices of different datasets may have different statistical properties. Since consensus is defined as teh component-wise minimum of two-TOMs, a bias may results. Simple scaling can mitigate the effects of different statistical properties to some degree. TOM are scales such that 95th percentile equals the 95th percentile of the female TOM.

  # reference: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.pdf

  # Define the reference percentile
  scaleP <- reference.percentile

  # Set RNG seed for reproducibility of sampling
  set.seed(12345)

  # Sample sufficiently large number of TOM entries
  nSamples <- as.integer(1/(1-scaleP) * 1000);

  # Choose the sampled TOM entries
  nGenes <- ncol(reference.TOM)
  scaleSample <- sample(nGenes*(nGenes-1)/2, size = nSamples)

  # Select the sampled TOM entries
  TOMScalingSample.ref <- as.dist(reference.TOM)[scaleSample]
  TOMScalingSample.query <- as.dist(query.TOM)[scaleSample]

  # Calculate the 95th percentile
  scaleQuant.ref <- quantile(TOMScalingSample.ref, probs = scaleP, type = 8)
  scaleQuant.query <- quantile(TOMScalingSample.query, probs = scaleP, type = 8)

  # Get scaling power
  scalePowers <- log(scaleQuant.ref)/log(scaleQuant.query);

  # scale TOM matrices
  scaled.TOM <- query.TOM^scalePowers

  return(scaled.TOM)
}



#' Remove duplicate genes from Seurat Object
#'
#' Remove duplicate genes from Seurat Object
#'
#' @param so Seurat Object
#' @name rmDuplicateGenes
#' @return seurat object
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' so.query <- rmDuplicateGenes(so.query)
#'
rmDuplicateGenes <- function(so){

  # get gene names
  all.genes <- rownames(so)

  # if duplicates exist, subset seurat object
  if (sum(duplicated(all.genes)) > 0){
    # which.dup <- all.genes[duplicated(all.genes)]
    which.unique <- all.genes[!duplicated(all.genes)]
    so <- subset(so, features = which.unique)
  }

  return(so)

}


#' Quantile Normalization of 2 Vectors
#'
#' Performs quantile normalization of 2 vectors (Hicks 2014).
#'
#' @param x Numeric vector.
#' @param y Numeric vector.
#' @param genes Character vector of gene names. Used to label entries in data.frame output.
#' @param flag.top.n Numeric indicating top n genes to flag in data.frame output. Default is 15.
#' @name qNorm
#' @author Nicholas Mikolajewicz
#' @return data.frame of quantile normalized values.
#' @examples
#'
#' # quantile normalization
#' x.s1 <- getConnectivity(s.1, gene.names = colnames(a.1))$wi
#' y.s2 <- getConnectivity(s.2, gene.names = colnames(a.2))$wi
#' df.xy.s <- qNorm(x.s1, y.s2, genes = colnames(a.1))
#'
qNorm <- function(x, y, genes = NULL, flag.top.n = 15){

  # check input
  stopifnot(length(x) == length(y))

  # get original vector orders
  x.orig.order <- seq(1,length(x))
  y.orig.order <- seq(1,length(y))

  # get new vector orders
  x.new.order <-order(x)
  y.new.order <-order(y)

  # sort vector
  x.sort <- x[x.new.order]
  y.sort <- y[y.new.order]

  # get average
  xy.mean <- (x.sort + y.sort)/2

  # restore original order
  x.new <- xy.mean[x.orig.order[x.new.order]]
  y.new <- xy.mean[y.orig.order[y.new.order]]


  if (!is.null(genes)){
    stopifnot(length(genes) == length(x))
    df.xy <- data.frame(genes = genes, x.old = x, y.old = y, x.new = x.new, y.new = y.new)
  } else {
    df.xy <- data.frame(x.old = x, y.old = y, x.new = x.new, y.new = y.new)
  }

  # flag.top.n <- 10
  df.xy$top.old <- F
  df.xy$top.old[rank(df.xy$x.old) >  nrow(df.xy)-flag.top.n] <- T
  df.xy$top.old[rank(df.xy$y.old) >  nrow(df.xy)-flag.top.n] <- T
  df.xy$top.old[rank((df.xy$x.old + df.xy$y.old)/2) >  nrow(df.xy)-flag.top.n] <- T

  df.xy$top.new <- F
  df.xy$top.new[rank(df.xy$x.new) >  nrow(df.xy)-flag.top.n] <- T
  df.xy$top.new[rank(df.xy$y.new) >  nrow(df.xy)-flag.top.n] <- T
  df.xy$top.new[rank((df.xy$x.new + df.xy$y.new)/2) >  nrow(df.xy)-flag.top.n] <- T

  return(df.xy)

}


#' Jaccard Similarity
#'
#' Computes Jaccard Similarity between two sets, x1 and x2.
#'
#' @param x1 set 1
#' @param x2 set 2
#' @param assert.unique Logical flag specifying whether to ensure that x1 and x2 are unique. Default is TRUE.
#' @name getJaccard
#' @return Jaccard similarity score (numeric)
#' @seealso \code{\link{pheatmap}}
#' @examples
#'
#' # compute jaccard similarity matrix for (named) list of genesets.
#' n.sets <- length(gene.sets)
#' j.mat <- matrix(nrow = n.sets, ncol = n.sets)
#'for (i in 1:n.sets){
#'   for (j in 1:n.sets){
#'     i.name <- names(gene.sets)[i]
#'     j.name <- names(gene.sets)[j]
#'     j.mat[i, j] <- getJaccard(gene.sets[[i.name]], gene.sets[[j.name]])
#'   }
#' }
#'
#' # assign row and col names
#' rownames(j.mat) <- names(gene.sets)
#' colnames(j.mat) <- names(gene.sets)
#'
#' # generate heatmap
#' pheatmap::pheatmap(j.mat, show_colnames = F, main = "Jaccard Similarity")
#'
getJaccard <- function(x1, x2, assert.unique = T){

  # ensure each set has unique entries
  if (assert.unique){
    x1 <- unique(x1)
    x2 <- unique(x2)
  }

  # compute intersection
  x.I <- length(intersect(x1, x2))

  # compute similarity
  x.S <- x.I/(length(x1) + length(x2) - x.I)

  # return result
  return(x.S)
}

#' Jaccard Similarity Matrix
#'
#' Computes Jaccard similarity mtrix for list of genesets
#'
#' @param gene.sets named list of genesets, where names specify name of gene set, and entries are character vectors specifying genes belongs to the respective set.
#' @param assert.unique Logical flag specifying whether to remove duplicate entries within individual sets. Default is TRUE.
#' @name jaccardSimilarityMatrix
#' @return Jaccard similarity matrix
#' @seealso \code{\link{pheatmap}}
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # compute jaccard similarity matrix for (named) list of genesets.
#' j.mat <- jaccardSimilarityMatrix(gene.sets)
#'
#' # generate heatmap
#' pheatmap::pheatmap(j.mat, show_colnames = F, main = "Jaccard Similarity")
#'
jaccardSimilarityMatrix <- function(gene.sets, assert.unique = T){
  n.sets <- length(gene.sets)
  j.mat <- matrix(nrow = n.sets, ncol = n.sets)
  for (i in 1:n.sets){
    for (j in 1:n.sets){
      i.name <- names(gene.sets)[i]
      j.name <- names(gene.sets)[j]
      j.mat[i, j] <- scMiko::getJaccard(gene.sets[[i.name]], gene.sets[[j.name]], assert.unique = assert.unique)
    }
  }

  rownames(j.mat) <- names(gene.sets)
  colnames(j.mat) <- names(gene.sets)

  return(j.mat)
}


#' Get cluster centers
#'
#' For given dimensional reduction (e.g., UMAP), return cluster center coordinates.
#'
#' @param df data.frame with columns 'x' and 'y' specifying dimensional reduction coordinates, and column 'cluster' or 'clusters' specifying cluster membership.
#' @param which.mean Specify which central values to use. One of 'mean' or 'median'. Default is 'mean'.
#' @return df.centers
#' @name getClusterCenters
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # compute cluster centers
#' which.center <- "mean"
#' df.red.centers <- get.cluster.centers(df.red, which.center)
#'
getClusterCenters <- function(df, which.center = "mean"){

  if (!("x" %in% colnames(df))) stop("x column is not specified")
  if (!("y" %in% colnames(df))) stop("y column is not specified")
  if (!(("cluster" %in% colnames(df)) | ("clusters" %in% colnames(df)))) stop("cluster column is not specified")

  # helper function for computing central value
  get.center <- function(x, which.center = "mean"){
    if (which.center == "mean" ) return(mean(x))
    if (which.center == "median" ) return(median(x))
  }

  # compute cluster centers
  df.centers <- df %>%
    group_by(cluster) %>%
    summarize(x.center = get.center(x, which.center),
              y.center = get.center(y, which.center))

  return(df.centers)
}


#' Identify pseudotime-dependent genes using Random Forest (RF) Model
#'
#' Identify pseudotime-dependent genes using Random Forest (RF) Model. Given Seurat object and highly variable genes (hvg), expression data are split into training and test set, and random forest model is trained to identify genes that vary with pseudotime. Model performance is evaluated on test set. Random forest model is fit using R 'parsnip' package.
#'
#' @param so Seurat Object
#' @param hvg Genes used to fit model (character vector; must be available in rows of seurat object). It is suggested to keep number of genes low (~200) for optimal performance.
#' @param pseudotimes Numeric vector of pseudotimes. Length must be equal to number of cells in seurat object (ncol(so)).
#' @param lineage.name Name of pseudotime lineage; used to label results.
#' @param slot A character specifying which slot to pull data from; default is 'Data'
#' @param assay A character specifying which assay to use (e.g., 'RNA' or 'SCT'). If unspecified, set to DefaultAssay(so)
#' @param mtry An integer for the number of predictors that will be randomly sampled at each split when creating the tree models.
#' @param trees An integer for the number of trees contained in the ensemble.
#' @param min_n An integer for the minimum number of data points in a node that are required for the node to be split further.
#' @param mode Specfiy type of RF to fit: 'regression' or 'classification'. Regression is default and it is not recommended to change this argument.
#' @param importance Type of importance. Default is 'impurity'.
#' @param num.threads An integer for the number of threads to use when fitting RF model
#' @return List of results
#' @seealso \code{\link{rand_forest}}
#' @name pseudotimeRF
#' @examples
#'
pseudotimeRF <- function(so, hvg, pseudotimes, lineage.name, slot = "data", assay = DefaultAssay(so),
                             mtry = length(hvg)/10, trees = 1000, min_n = 15, mode = "regression", importance = "impurity", num.threads = 3){

  # get data for highly variable genes (hvg)
  cur.data <- GetAssayData(so, slot = slot, assay = assay)
  match.ind <- which(rownames(cur.data) %in% top_hvg)
  dat_use <- as.data.frame(t(GetAssayData(so, slot = slot, assay = assay)[match.ind,]))

  # merge expression data and pseudotime
  dat_use_df <- cbind(pseudotimes, dat_use)
  colnames(dat_use_df)[1] <- "pseudotime"
  dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

  # Define training, testing and validation sets
  colnames(dat_use_df) <- make.names(colnames(dat_use_df) , unique = TRUE, allow_ = TRUE)
  dat_split <- initial_split(dat_use_df)
  dat_train <- training(dat_split)
  dat_val <- testing(dat_split)

  # Train Model
  model <- rand_forest(mtry =mtry, trees = trees, min_n = min_n, mode = mode) %>%
    set_engine("ranger", importance =importance, num.threads = num.threads) %>%
    fit(pseudotime ~ ., data = dat_train)

  # Evaluate Model
  val_results <- dat_val %>%
    mutate(estimate = predict(model, .[,-1]) %>% pull()) %>%
    select(truth = pseudotime, estimate)
  model.metrics <- metrics(data = val_results, truth, estimate)

  # store results
  # model.list[[lineage.name]] <- model
  # results.list[[lineage.name]] <- val_results
  df.metrics <- data.frame(lineage = lineage.name,
                           rmse = signif(model.metrics[[".estimate"]][1],3),
                           rsq = signif(model.metrics[[".estimate"]][2],3),
                           mae = signif(model.metrics[[".estimate"]][3],3))

  output <- list(
    model = model,
    prediction = val_results,
    performance = df.metrics
  )

  return(output)

}


#' Infer initial trajectory through space
#'
#' `inferInitialTrajectory` infers an initial trajectory for  `princurve::principal_curve` by clustering the points and calculating the shortest path through cluster centers. The shortest path takes into account the euclidean distance between cluster centers, and the density between those two points. Based on `infer_initial_trajectory` function from `SCORPIUS` package
#'
#' @param space A numeric matrix or a data frame containing the coordinates of samples.
#' @param k The number of clusters
#' @return the initial trajectory obtained by this method
#' @name inferInitialTrajectory
#' @seealso \code{\link{subsetDimRed}}
#' @examples
#'
#' # specify features (i.e., clusters of interest)
#' which.clusters <- c(0,4,5)
#'
#' # get dimensional reduction for specified clusters
#' dimSubset <- subsetDimRed(so.query, which.features = which.clusters)
#'
#' # get initial trajectory
#' start.traj <- inferInitialTrajectory(as.matrix(dimSubset[["reduction"]]), k = length(unique(dimSubset[["features"]])))
#'
inferInitialTrajectory <- function (space, k) {
  # check_numeric_matrix(space, "space", finite = TRUE)
  # check_numeric_vector(k, "k", whole = TRUE, finite = TRUE,
  # range = c(1, nrow(space) - 1), length = 1)
  fit <- stats::kmeans(space, k)
  centers <- fit$centers
  eucl_dist <- as.matrix(stats::dist(centers))
  i <- j <- NULL
  pts <- crossing(i = seq_len(k), j = seq_len(k), pct = seq(0,
                                                            1, length.out = 21)) %>% filter(i < j)
  pts_space <- (1 - pts$pct) * centers[pts$i, ] + pts$pct *
    centers[pts$j, ]
  pts$dist <- rowMeans(RANN::nn2(space, pts_space, k = 10)$nn.dist)
  dendis <- pts %>% group_by(i, j) %>% summarise(dist = mean(dist)) %>%
    ungroup()
  density_dist <- matrix(0, nrow = k, ncol = k)
  density_dist[cbind(dendis$i, dendis$j)] <- dendis$dist
  density_dist[cbind(dendis$j, dendis$i)] <- dendis$dist
  cluster_distances <- eucl_dist * density_dist
  tsp <- TSP::insert_dummy(TSP::TSP(cluster_distances))
  tour <- as.vector(TSP::solve_TSP(tsp))
  tour2 <- c(tour, tour)
  start <- min(which(tour2 == k + 1))
  stop <- max(which(tour2 == k + 1))
  best_ord <- tour2[(start + 1):(stop - 1)]
  init_traj <- centers[best_ord, , drop = FALSE]
  init_traj
}


#' Get dimensional reduction from Seurat Object for subset of data
#'
#' Get dimensional reduction from Seurat Object for subset of data.
#'
#' @param so Seurat Object
#' @param which.features Character vector of features to subset on.
#' @param groups A character specifying metadata column name (in Seurat object) that contains features of interest. \
#' @param reduction A character specifying which reduction to retrieve from Seurat Object. Default is 'umap'.
#' @return list of results
#' @author Nicholas Mikolajewicz
#' @name subsetDimRed
subsetDimRed <- function(so, which.features, groups = "seurat_clusters", reduction = "umap"){

  if (!(tolower(reduction) %in% tolower(names(so@reductions)) )) stop(paste0(reduction, " does not exist"))
  df.meta <- so@meta.data
  if (!(groups %in% colnames(df.meta) )) stop(paste0(groups, " does not exist"))

  match.ind <- which(df.meta[ ,groups] %in% which.features)
  u.group <- as.character(unique(df.meta[match.ind ,groups]))

  # order groups
  if (is.factor(unique(df.meta[match.ind ,groups]))){
    u.group <-  levels(unique(df.meta[match.ind ,groups]))[levels(unique(df.meta[match.ind ,groups])) %in% u.group]
  }

  df.reduction <- data.frame(so@reductions[[reduction]]@cell.embeddings)

  df.red.cc <- data.frame(df.reduction[ ,c(1,2)])
  colnames(df.red.cc) <- c("x", "y")
  df.red.cc$cluster <- df.meta[ ,groups]
  df.red.cc <- df.red.cc[ match.ind,]

  df.centers <- getClusterCenters(df.red.cc, which.center = "mean")
  df.reduction <- df.reduction[ match.ind,]

  output <- list(
    groups = groups,
    features = df.meta[match.ind ,groups],
    reduction = df.reduction,
    centers = df.centers
  )

  return(output)

}



#' Get lineage trajectories using principal curves and compute pseudotimes
#'
#' For given space (e.g., UMAP, PCA, etc.), lineage trajectories are fit using prinicpal curves, and these are then used to derive pseudotimes.
#'
#' @param space A numeric matrix or a data frame containing the coordinates of samples.
#' @param start either a previously fit principal curve, or else a matrix of points that in row order define a starting curve. If missing or NULL, then the first principal component is used. If the smoother is "periodic_lowess", then a circle is used as the start.
#' @param group.labels Character vector (same length as number of rows in space) specifying group membership. Optional.
#' @param pseudotimes Either logical (TRUE) specifying whether to use internally generated pseudotimes (from principal curve fitting) to color generated umap points, or, numeric vector with pseudotimes. If provided, plotted data are colored by pseudotime. Otherwise colored by group membership.
#' @param thresh convergence threshold on shortest distances to the curve. Default is 0.001.
#' @param maxit maximum number of iterations.
#' @param stretch A stretch factor for the endpoints of the curve, allowing the curve to grow to avoid bunching at the end. Must be a numeric value between 0 and 2.
#' @param smoother choice of smoother. The default is "smooth_spline", and other choices are "lowess" and "periodic_lowess". The latter allows one to fit closed curves. Beware, you may want to use iter = 0 with lowess().
#' @param approx_points Approximate curve after smoothing to reduce computational time. If FALSE, no approximation of the curve occurs. Otherwise, approx_points must be equal to the number of points the curve gets approximated to; preferably about 100.
#' @param trace If TRUE, the iteration information is printed
#' @param plot_iteractions If TRUE the iterations are plotted.
#' @return list of results containing principal curve fits and coordinates, pseutimes, plots
#' @seealso \code{\link{inferInitialTrajectory}}
#' @name lineageTrajectory
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' # get lineage name
#' lineage.name <- names(ss.lineages)[1]
#'
#' # get dimnesional reduction for subset of clusters belonging to lineage
#' dimSubset <- subsetDimRed(so.query, which.features = ss.lineages[[lineage.name]])
#'
#' # get initial trajectory
#' start.traj <- inferInitialTrajectory(as.matrix(dimSubset[["reduction"]]), k = length(unique(dimSubset[["features"]])))
#'
#' # get lineage trajectories
#' LT.results <- lineageTrajectory(space = as.matrix(dimSubset[["reduction"]]),
#'               start = start.traj,
#'               group.labels = dimSubset[["features"]]
#'
#'
#'
lineageTrajectory <- function(space, start = NULL, group.labels = NULL, pseudotimes = NULL, thresh = 0.001, maxit = 10, stretch = 2, smoother = "smooth_spline", approx_points = 100, trace = FALSE, plot_iterations = FALSE){


  # fit prinicpal curves
  fit <- princurve::principal_curve(as.matrix(space), start = start,
                                    thresh = thresh, maxit = maxit, stretch = stretch, smoother = smoother,
                                    approx_points = approx_points, trace = trace, plot_iterations = plot_iterations)
  traj.path <- fit$s[fit$ord, , drop = FALSE]

  # get trajectories
  df.traj <- data.frame(traj.path)

  # get pseudotimes
  ps <- dynutils::scale_minmax(fit$lambda)

  # generate plots
  orig.names <- colnames(space)
  colnames(space)[c(1,2)] <- c("x", "y")
  colnames(df.traj)[c(1,2)] <- c("x", "y")




  if (!is.null(group.labels) & is.null(pseudotimes)){
    df.space <- data.frame(space, group = group.labels)
    plt.space <- df.space %>% ggplot() + geom_point(aes(x,y, color = group))
  } else if (class(pseudotimes) == "logical"){
    if (pseudotimes){
      df.space <- data.frame(space, pt = ps)
      plt.space <- df.space %>% ggplot() + geom_point(aes(x,y, color = pt)) + scale_color_viridis("pseudotime")
    } else {
      df.space <- data.frame(space)
      plt.space <- df.space %>% ggplot() + geom_point(aes(x,y))
    }
  } else if (class(pseudotimes) == "numeric"){
    df.space <- data.frame(space, pt = pseudotimes)
    plt.space <- df.space %>% ggplot() + geom_point(aes(x,y, color = pt)) + scale_color_viridis("pseudotime")
  } else {
    df.space <- data.frame(space)
    plt.space <- df.space %>% ggplot() + geom_point(aes(x,y))
  }

  plt.trajectory <- plt.space +
    geom_path(data = df.traj, aes(x, y), size = 1) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    xlab(orig.names[1]) +
    ylab(orig.names[2])


  output <- list(
    principal.curve.fit = fit,
    principal.curve.coordinates = df.traj,
    pseudotime = ps,
    plt.trajectory = plt.trajectory,
    space = space,
    groups = group.labels,
    df.space = df.space
  )

  return(output)


}




#' Variance explained by each principal component.
#'
#' For given Seurat Object, retrieve principal components and compute proportion of explained variance.
#'
#' @param so
#' @return data.frame summarize proportion of variance explained by each principal component
#' @author Nicholas Mikolajewicz
#' @name propVarPCA
propVarPCA <- function(so){

  # get pca reduction
  pc.std <- so@reductions[["pca"]]@stdev

  # variance explained
  pc.var <- pc.std^2

  # proportion
  pc.prop_var <- pc.var/sum(pc.var)

  # cumulative
  pc.cum_sum <- cumsum(pc.prop_var)

  # store results
  pc.id <- c(1:length(pc.std))
  df.pca <- data.frame(pc.id, pc.prop_var, pc.cum_sum)

  return(df.pca)
}


#' Get branch-specific pseudotimes from diffusion map
#'
#' Get branch-specific pseudotimes from diffusion map. Adopted from destiny package.
#'
#' @param dpt
#' @param branch_id
#' @return numeric vector; pseudotimes
#' @name dpt_for_branch
#' @seealso \code{\link{DPT}}
#' @examples
#'
#' # PCA embedding
#' so.pca <- so.query@reductions[["pca"]]@cell.embeddings
#'
#' # diffusion map
#' dm <- DiffusionMap(so.pca, n_eigs = 50)
#'
#' # diffusion pseudotime
#' dpt <- DPT(dm, tips =tip.cluster)
#'
#' # flatten
#' dpt_flat <- branch_divide(dpt, divide = 1)
#'
#' # get branch specific pseudotime
#' pt_vec <- dpt_for_branch(dpt_flat, branch_id = paths_to)
#'
dpt_for_branch <- function(dpt, branch_id) {
  branch_idx <- dpt@branch[, 1L] == branch_id
  stopifnot(any(branch_idx))
  tip_cells <- which(branch_idx & dpt@tips[, 1L])
  if (length(tip_cells) == 0L) tip_cells <- which(branch_idx)
  dpt[tip_cells[[1L]], ]
}


#' Infer initial trajectory through space, using pseudotime priors
#'
#' Infer initial trajectory through space, using pseudotime priors. Adopted from destiny package, uses smth.gaussian function; The specific function for smoothing using the gaussian window function.
#'
#' @param pt pseudotimes. Numeric vector.
#' @param space A numeric matrix or a data frame containing the coordinates of samples.
#' @param w.width the length of the smoothing window, if an integer, represents number of items, else, if a value between 0 and 1, represents the proportion of the input vector
#' @return the initial trajectory
#' @name inferInitialTrajectory.v2
#' @seealso \code{\link{principal_curve}}
#' @importFrom smoother smth.gaussian
#' @importFrom graphics plot
#' @examples
#'
#'  # lineage-specific pseudotime
#'  pt_vec <- dpt_for_branch(dpt_flat, branch_id = paths_to)
#'
#' # get indices
#' idx <- dpt_flat@branch[, 1] %in% c(root, paths_to)
#'
#' # get initial trajectory path
#' umap.path <- inferInitialTrajectory.v2(pt_vec[idx], so.query@reductions[["umap"]]@cell.embeddings[idx, ], w_width = 0.15)
#'
#'    # get umap coordinates
#' df.umap.all <- data.frame(x = so.query@reductions[["umap"]]@cell.embeddings[ ,1],
#'                           y = so.query@reductions[["umap"]]@cell.embeddings[ ,2])
#' df.umap.all$col <- "grey"
#'
#' df.umap.sub <- df.umap.all[idx, ]
#' df.umap.sub$pt <- pt_vec[idx]
#'
#' # compute prinicpal curves
#' pc.fit <- princurve::principal_curve(as.matrix(df.umap.sub[,c("x", "y")]), start = as.matrix(umap.path[ ,c(1,2)]),
#'                                      thresh = 0.001, maxit = 10, stretch = 2, smoother = "smooth_spline",
#'                                      approx_points = 100, trace = FALSE, plot_iterations = FALSE)
#' traj.path <- pc.fit$s[pc.fit$ord, , drop = FALSE]
#' df.tp <- data.frame(traj.path)
#'
inferInitialTrajectory.v2 <- function(pt, space, w_width = .1) {
  stopifnot(identical(nrow(space), length(pt)))
  as.data.frame(apply(space[order(pt), ], 2, function(col) smoother::smth.gaussian(col, w_width, tails = TRUE)))
}


#' Get cell index corresponding to center of specified root cluster.
#'
#' Get cell index corresponding to center of specified root cluster. For specified cluster, cluster center is computed and index of nearest cell to cluster center is returned.
#'
#' @param x x coordinates (e.g., UMAP 1)
#' @param y y coordinates (e.g., UMAP 2)
#' @param cluster.membership group id's corrorespodning to cluster membership (e.g., seurat_cluster entries). Must be same length as x and y.
#' @param which.cluster specify which cluster to get root cell fot.
#' @return root cell index
#' @name getClusterRoot
#' @author Nicholas Mikolajewicz
#' @examples
#'
getClusterRoot <- function(x, y, cluster.membership, which.cluster){

  # assertions
  stopifnot(length(which.cluster) == 1)
  stopifnot(length(cluster.membership) == length(y))
  stopifnot(length(cluster.membership) == length(x))

  # data.frame
  df.cc <- data.frame(x, y, cluster = cluster.membership)

  # get centers
  df.centers <- getClusterCenters(df.cc, which.center = "mean")

  # filter
  df.coi <- df.centers[df.centers$cluster %in% which.cluster, ]

  # get nearest cell to center
  df.cc$x.dif <- (df.cc[,1] - df.coi$x.center)^2
  df.cc$y.dif <- (df.cc[,2] - df.coi$y.center)^2
  df.cc$xtDist <- sqrt(df.cc$x.dif + df.cc$y.dif)
  root_cell <- which.min( df.cc$xtDist)

  return(root_cell)
}

#' Ensure that all dimNames are correctly specified in Seurat Object
#'
#' Ensure that all dimNames are correctly specified in Seurat Object. If incorrectly specified, subsetting functions will lead to unexpected results.
#'
#' @param so Seurat Object
#' @return Seurat Object
#' @name updateDimNames
#' @author Nicholas Mikolajewicz
#' @examples
#'
#'so <- UpdateDimNames(so)
#'
updateDimNames <- function(so){

  # get all assays
  all.assays <- names(so@assays)

  # ensure dim names are correctly specified
  for (i in 1:length(all.assays)){
    c.dim <- so@assays[[all.assays[i]]]@counts@Dimnames
    d.dim <- so@assays[[all.assays[i]]]@data@Dimnames

    if ((length(c.dim) == 1) & (length(d.dim) == 2)){
      if ((so@assays[[all.assays[i]]]@counts@Dim[1] == so@assays[[all.assays[i]]]@data@Dim[1]) &
          (so@assays[[all.assays[i]]]@counts@Dim[2] == so@assays[[all.assays[i]]]@data@Dim[2])){
        so@assays[[all.assays[i]]]@counts@Dimnames <- so@assays[[all.assays[i]]]@data@Dimnames
      }
    } else if ((length(d.dim) == 1) & (length(c.dim) == 2)){
      if ((so@assays[[all.assays[i]]]@counts@Dim[1] == so@assays[[all.assays[i]]]@data@Dim[1]) &
          (so@assays[[all.assays[i]]]@counts@Dim[2] == so@assays[[all.assays[i]]]@data@Dim[2])){
        so@assays[[all.assays[i]]]@data@Dimnames <- so@assays[[all.assays[i]]]@counts@Dimnames
      }
    }
  }

  return(so)
}




#' Filter seurat object by specified cluster ids
#'
#' Filter seurat object by specified cluster ids
#'
#' @param so Seurat Object
#' @param include Cluster ids to include. Optional, NULL if unspecified.
#' @param omit Cluster ids to omit. Optional, NULL if unspecified.
#' @param which.field Metadata field to filter by. Default is 'seurat_clusters'
#' @return Seurat Object
#' @author Nicholas Mikolajewicz
#' @name clusterFilter
#' @examples
#'
#' # specify filtering parameters
#' filter.parameters <- list(
#'      include = NULL,
#'      omit = c(8,9)
#'      )
#'
#' # filter seurat
#' so.filtered <- clusterFilter(so, include = filter.parameters$include, omit = filter.parameters$omit)
#'
clusterFilter <- function(so, include = NULL, omit = NULL, which.field = "seurat_clusters"){

  # ensure seurat dim names are up to date to ensure proper subsetting
  so <- updateDimNames(so)

  # get seurat meta data
  df.meta <- so@meta.data
  df.meta$cells <- rownames(df.meta)

  # which cells to include (according to inclusion parameters)
  if (!is.null(include)){
    include.which.include <- df.meta$cells[as.character(df.meta[,which.field]) %in% as.character(include)]
  } else {
    include.which.include <- NULL
  }

  # which cells to include (according to omission parameters)
  if (!is.null(omit)){
    include.which.omit <-df.meta$cells[!(as.character(df.meta[,which.field]) %in% as.character(omit))]
  } else {
    include.which.omit <- NULL
  }

  # combine inclusion and omission indices
  if (is.null(include.which.include) & is.null(include.which.omit)){
    include.which.all <- df.meta$cells
  } else {
    include.which.all <- unique(c(include.which.include, include.which.omit))
  }

  # subset and return seurat object
  return(SubsetData(so, cells = WhichCells(so, cells = include.which.all)))
}


# quickly choose an elbow for a PC.
# at variance below 5% per component, choose the largest % drop
# designed for variance percentages, but will also work given a full set of Evalues
#' Quickly estimate the 'elbow' of a scree plot (PCA)
#'
#' This function uses a rough algorithm to estimate a sensible 'elbow' to
#' choose for a PCA scree plot of eigenvalues. The function looks at an initial arbitrarily 'low'
#' level of variance and looks for the first eigenvalue lower than this. If the very first eigenvalue
#' is actually lower than this (i.e, when the PCs are not very explanatory) then this 'low' value is
#' iteratively halved until this is no longer the case. After starting below this arbitrary threshold
#' the drop in variance explained by each pair of consecutive PCs is standardized by dividing over the
#' larger of the pair. The largest percentage drop in the series below 'low' % is selected as the 'elbow'.
#' @param varpc numeric, vector of eigenvalues, or 'percentage of variance' explained datapoints for
#'  each principle component. If only using a partial set of components, should first pass to
#'  estimate.eig.vpcs() to estimate any missing eigenvalues.
#' @param low numeric, between zero and one, the threshold to define that a principle component
#'  does not explain much 'of the variance'.
#' @param max.pc maximum percentage of the variance to capture before the elbow (cumulative sum to PC 'n')
#' @return The number of last principle component to keep, prior to the determined elbow cutoff
#' @export
#' @author Nicholas Cooper
#' @name pcaElbow
#' @examples
#' # correlated data
#' mat <- sim.cor(100,50)
#' result <- princomp(mat)
#' eig <- result$sdev^2
#' elb.a <- quick.elbow(eig)
#' pca.scree.plot(eig,elbow=elb.a,M=mat)
#' elb.b <- quick.elbow(eig,low=.05) # decrease 'low' to select more components
#' pca.scree.plot(eig,elbow=elb.b,M=mat)
#' # random (largely independent) data, usually higher elbow #
#' mat2 <- generate.test.matrix(5,3)
#' result2 <- princomp(mat2)
#' eig2 <- result2$sdev^2
#' elb2 <- quick.elbow(result2$sdev^2)
#' pca.scree.plot(eig2,elbow=elb2,M=mat2)
pcaElbow <- function(varpc,low=.08,max.pc=.9) {
  ee <- varpc/sum(varpc) # ensure sums to 1
  #print(round(log(ee),3))
  while(low>=max(ee)) { low <- low/2 } # when no big components, then adjust 'low'
  lowie <- (ee<low) ; highie <- ee>low/8
  low.ones <- which(lowie & highie)
  others <- length(which(!lowie))
  if(length(low.ones)>0) {
    if(length(low.ones)==1) {
      elbow <- low.ones
    } else {
      set <- ee[low.ones]
      pc.drops <- abs(diff(set))/(set[1:(length(set)-1)])
      infz <- is.infinite(pc.drops)
      #print(pc.drops)
      elbow <- which(pc.drops==max(pc.drops[!infz],na.rm=T))[1]+others
    }
  } else {
    # if somehow there are no small eigenvalues, just choose the elbow as the second last
    warning("no eigenvalues were significantly smaller than the previous\n")
    elbow <- length(ee)
  }
  if(tail(cumsum(ee[1:elbow]),1)>max.pc) {
    elbow <- which(cumsum(ee)>max.pc)[1]-1
  }
  if(elbow<1) {
    warning("elbow calculation failed, return zero")
    return(0)
  }
  names(elbow) <- NULL
  return(elbow)
}





#' Identify temporally varying pathways using Tempora
#'
#' Identify temporally varying pathways using Tempora
#'
#' @param object Tempora Object
#' @param pval_threshold P-value threshold to determine the significance of pathway enrichment over time. Default is 0.05.
#' @param adjust.p Flag to adjust p value using BH correction. Default is TRUE.
#' @param pathway.filter Numeric specifying z score threshold used to pre-filter pathways. Default is 3.
#' @return List containing generalized additive model (GAMs) fits, plots and p-values
#' @author Gary Bader et. al
#' @name varyingPaths.Tempora
#' @examples
#'
varyingPaths.Tempora <- function (object, pval_threshold = 0.05, adjust.p = T, pathway.filter = 3) {

  # initiate results list
  tp.list <- list()

  if (class(object)[1] != "Tempora") {
    stop("Not a valid Tempora object")
  }
  if (is.null(object@n.pcs)) {
    stop("BuildTrajectory has not been run. See ?Tempora::BuildTrajectory for details")
  }
  if (is.null(object@cluster.pathways)) {
    stop("CalculatePWProfiles has not been run. See ?Tempora::CalculatePWProfiles for details")
  }
  gsva_bycluster <- object@cluster.pathways
  significant_pathways <- c()

  for (i in 1:object@n.pcs) {
    genes_scaled <- scale(object@cluster.pathways.dr$rotation[,
                                                              i])
    significant_pathways <- c(names(which(genes_scaled[, 1] > pathway.filter |
                                            genes_scaled[, 1] < -pathway.filter)), significant_pathways)
  }

  tp.list$genes.scaled <- genes_scaled
  tp.list$n.pathways.considered <- length(significant_pathways)
  tp.list$considered.pathways <- significant_pathways

  pca_pathways <- sub("%.*", "", significant_pathways)
  pca_pathways_cleaned <- gsub("[[:punct:]]", "", pca_pathways)
  themes <- pca_pathways_cleaned
  warning("Fitting GAM models...")
  p_vals <- gams <- list()

  for (i in 1:length(themes)) {
    if (length(grep(themes[i], rownames(gsva_bycluster))) > 1) {
      plot_df <- data.frame(cluster = colnames(gsva_bycluster[grep(themes[i],
                                                                   rownames(gsva_bycluster)), ]),
                            value = colMeans(gsva_bycluster[grep(themes[i],
                                                                 rownames(gsva_bycluster)), ], na.rm = T))
    }
    else if (length(grep(themes[i], rownames(gsva_bycluster))) == 1) {
      plot_df <- data.frame(cluster = names(gsva_bycluster[grep(themes[i],
                                                                rownames(gsva_bycluster)), ]),
                            value = gsva_bycluster[grep(themes[i],
                                                        rownames(gsva_bycluster)), ])
    } else {
      next
    }
    plot_df$time <- object@cluster.metadata$Cluster_time_score
    gams[[themes[i]]] <- mgcv::gam(value ~ s(time, k = 3, bs = "cr"),
                                   data = plot_df)
    temp_anova <- mgcv::anova.gam(gams[[themes[i]]])
    p_vals[[themes[i]]] <- temp_anova$s.pv
  }

  if (adjust.p){
    p_vals_adj <- p.adjust(unlist(p_vals[which(unlist(p_vals) >
                                                 0)]), method = "BH")
    p.val.label <- "P-value = "
  } else {
    p_vals_adj <- unlist(p_vals[which(unlist(p_vals) > 0)])
    p.val.label <- "Adjusted p-value = "
  }

  varying_pathways <- p_vals_adj[which(p_vals_adj < pval_threshold)]
  varying_pathways <- varying_pathways[!duplicated(names(varying_pathways))]

  tp.list$p.vals <-  p_vals
  tp.list$p.vals.adj <- p_vals_adj
  tp.list$gams <- gams
  tp.list$varying.pathways <- varying_pathways

  plt.varying.pathway <- list()
  if (length(varying_pathways) == 0){

  } else {


    warning("Generating time-dependent pathways plots...")

    for (i in 1:length(varying_pathways)) {

      pathway.name <- varying_pathways[i]

      if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) > 1) {
        plot_df <- data.frame(cluster = colnames(gsva_bycluster[grep(names(varying_pathways)[i],
                                                                     rownames(gsva_bycluster)), ]),
                              value = colMeans(gsva_bycluster[grep(names(varying_pathways)[i],
                                                                   rownames(gsva_bycluster)), ]))
        plot_df$time <- object@cluster.metadata$Cluster_time_score
      } else if (length(grep(names(varying_pathways)[i], rownames(gsva_bycluster))) ==  1) {
        plot_df <- data.frame(cluster = names(gsva_bycluster[grep(names(varying_pathways)[i],
                                                                  rownames(gsva_bycluster)), ]),
                              value = gsva_bycluster[grep(names(varying_pathways)[i],
                                                          rownames(gsva_bycluster)), ])
        plot_df$time <- object@cluster.metadata$Cluster_time_score
      }

      id <- which(names(gams) == names(varying_pathways)[i])

      # get plotting data
      plot.list <- NULL
      plot.list <- R.devices::suppressGraphics({mgcv::plot.gam(gams[[id[1]]], main = paste0(names(varying_pathways)[i]),
                                                               xlab = "Inferred Time", ylab = "Pathway Expression Level",
                                                               bty = "l", cex.main = 1,  shade = F,
                                                               se = 3, scheme = 1)})

      # prep data as dataframe
      df.plot <- data.frame(x = plot.list[[1]][["x"]],
                            y = plot.list[[1]][["fit"]],
                            se = plot.list[[1]][["se"]],
                            ci = 1.96*plot.list[[1]][["se"]])

      # generate ggplot

      x.min <- min(df.plot$x, na.rm = T); x.max <- max(df.plot$x, na.rm = T)
      plt.varying.pathway[[i]] <- ggplot() +
        xlab("Inferred Time") +
        ylab("Pathway Activity") +
        theme_miko() +
        geom_ribbon(data = df.plot, aes(x, y, ymin = y-ci, ymax = y+ci), alpha = 0.3) +
        geom_line(data = df.plot, aes(x, y), size = 1, alpha = 0.8) +
        geom_text(data = plot_df, aes(x = time, y = value, label = cluster), size = 5) +
        labs(title = plot.list[[1]][["main"]], subtitle = paste0(p.val.label, round(varying_pathways[[i]], 5))) +
        scale_x_continuous(breaks=c(x.min, x.max), labels=c("Early", "Late"))

      # print(plt.varying.pathway[[i]])
    }
  }

  tp.list$plots <-  plt.varying.pathway

  return(tp.list)
}


#' Calcualte pathway enrichment scores in Tempora Framework
#'
#' Calcualte pathway enrichment scores in Tempora Framework. Uses GSVA backend.
#'
#' @param object Tempora Object
#' @param gmt_path Local path to database of pathways or genesets organized as a .gmt file. Genesets files in GMT format can be downloaded at http://baderlab.org/GeneSets.
#' @param method Method used to estimate pathway enrichment profile per cluster. Can be "gsva", "ssgsea", "zscore" or "plage", default to "gsva". See ?gsva for more information.
#' @param min.sz Minimum size of the genesets used in enrichment estimation, set to 5 genes by default.
#' @param max.sz Maximum size of the genesets used in enrichment estimation, set to 200 genes by default.
#' @param parallel.sz Type of cluster architecture when using snow. If 1, no parallelization will be used. If 0, all available cores will be used.
#' @param verbose Flag for reporting progress
#' @param do.plot Flag for ploting scree plot
#' @return Tempora object with internally calculated pathway scores
#' @author Gary Bader et. al
#' @name pathActivity.Tempora
#' @examples
#'
pathActivity.Tempora <- function (object, gmt_path, method = "gsva", min.sz = 5, max.sz = 200,
                                  parallel.sz = 1, verbose = F, do.plot = F)
{
  if (class(object)[1] != "Tempora") {
    stop("Not a valid Tempora object")
  }
  warning("Calculating cluster average gene expression profile...")
  exprMatrix <- object@data
  exprMatrix_bycluster <- list()
  pathwaygmt <- GSEABase::getGmt(gmt_path)

  for (i in sort(unique(object@meta.data$Clusters))) {
    exprMatrix_bycluster[[i]] <- rowMeans(exprMatrix[, which(colnames(exprMatrix) %in%
                                                               rownames(object@meta.data)[which(object@meta.data$Clusters ==
                                                                                                  i)])])
  }

  names(exprMatrix_bycluster) <- sort(unique(object@meta.data$Clusters))
  exprMatrix_bycluster <- do.call(cbind, exprMatrix_bycluster)
  colnames(exprMatrix_bycluster) <- sort(unique(object@meta.data$Clusters))
  rownames(exprMatrix_bycluster) <- rownames(exprMatrix)

  warning("\nCalculating cluster pathway enrichment profiles...\n")
  gsva_bycluster <- GSVA::gsva(as.matrix(exprMatrix_bycluster),
                               pathwaygmt, method = method, min.sz = min.sz, max.sz = max.sz,
                               parallel.sz = parallel.sz, verbose = verbose)

  colnames(gsva_bycluster) <- colnames(exprMatrix_bycluster)
  object@cluster.pathways <- gsva_bycluster

  gsva_bycluster_pca <- prcomp(t(gsva_bycluster), scale = T,
                               center = T)

  if (do.plot){
    screeplot(gsva_bycluster_pca, npcs = 25, type = "lines",
              main = "PCA on pathway enrichment analysis result")
  }

  object@cluster.pathways.dr <- gsva_bycluster_pca
  validObject(object)
  return(object)
}


#' Row-wise matrix binning
#'
#' Bin matrix row-wise by taking averages across regular row intervals.
#'
#' @param m matrix
#' @param bin.size Numeric indicating number of bins. Resulting matrix will have bin.size number of rows.
#' @return binned matrix
#' @author Nicholas Mikolajewicz
#' @name binMatrix
#' @examples
#'
binMatrix <- function(m, bin.size) {

  if (bin.size > nrow(m)) stop("bin.size cannot exceed number of rows in matrix")
  bm <- matrix(nrow = round((nrow(m)/bin.size)), ncol = ncol(m))
  for (i in 1:round((nrow(m)/bin.size))) {
    er <- (i*bin.size)
    sr <- er-(bin.size-1)
    if (er > nrow(m)) er <- nrow(m)
    bm[i,] <- colMeans(m[sr:er ,],na.rm=T)
  }
  colnames(bm) <- colnames(m)
  bm
}


#' Vector binning
#'
#' Bin vector by taking averages across regular intervals.
#'
#' @param v numeric vector
#' @param bin.size Numeric indicating number of bins. Resulting vector will have bin.size length.
#' @return binned vector
#' @author Nicholas Mikolajewicz
#' @name binVector
#' @examples
#'
binVector <- function(v, bin.size) {

  if (bin.size > length(v)) stop("bin.size cannot exceed length of vector")
  bv <- c()
  for (i in 1:round((length(v)/bin.size))) {
    er <- (i*bin.size)
    sr <- er-(bin.size-1)
    if (er > length(v)) er <- length(v)
    bv <- c(bv, mean(v[sr:er],na.rm=T))
  }

  bv
}



#' Returns top module genes from NMF feature loading matrix
#'
#' Returns top module genes from NMF feature loading matrix
#'
#' @param feature.loading Matrix; m genes x n modules loading
#' @param norm.cutoff Numeric [0,1]; cutoff to get top genes. Default is 0.5.
#' @author Nicholas Mikolajewicz
#' @name getNMFGenes
#' @examples
#'
#' # get feature loading matrix (W matrix)
#' nmf.kme <- as.matrix(so.query@misc[["nmf"]][["nmf"]][["W"]])
#'
#' # get module genes
#' nmf.module.genes <- getNMFGenes((nmf.kme), norm.cutoff = (1/ncol(nmf.kme)) + .25)
#'
getNMFGenes <- function(feature.loading, norm.cutoff = 0.5){

  if (!("matrix" %in% class(feature.loading))) stop("feature.loading input is not a matrix")

  nmf.kme <- t(feature.loading)
  nmf.kme <- (apply(nmf.kme, 2, function(x) (x/sum(x))))

  # get module genes
  module.genes <-  apply(nmf.kme, 1, function(x) colnames(nmf.kme)[x>norm.cutoff])

  module.size <- unlist(lapply(module.genes, length))
  nmf.module.genes <- module.genes[module.size > 0]

  return(nmf.module.genes)
}


#' Reinstall scMiko package
#'
#' Reinstall scMiko package from private repository (https://github.com/NMikolajewicz/scMiko). Package is reloaded after update.
#'
#' @param token autherization key to for private git repository
#' @param ... Additional arguments passed to devtools::install_github()
#' @author Nicholas Mikolajewicz
#' @name scMikoUpdate
#' @examples
#'
#' # reinstall scMiko
#' scMikoUpdate()
#'
scMikoUpdate <- function(token = "a3c1c9b15c496991c952d1fe3ccc52db770f22fa", ...){

  try({detach("package:scMiko", unload = T)}, silent = T)

  devtools::install_github(
    repo = "NMikolajewicz/scMiko",
    ref = "master",
    auth_token = token,
    ...)

  try({scMiko::scMikoReload()})
}


#' Update central log
#'
#' Update central log. Central log keeps track of every scPipeline run and includes information about the module run, data used, user, date.
#'
#' @param Module Module Name/ID. e.g., "M01" is used to specify Module 1.
#' @param clog.file Central log file name. A character.
#' @param log.path Path to central log file. Set to data.path specified in .RProfile by default.
#' @param user.id User name. Set to user specified in .RProfile by default.
#' @param run.notes Additional notes to describe current run.
#' @param pdf.flag Logical indicating whether PDF figures were generated.
#' @author Nicholas Mikolajewicz
#' @name updateCentralLog
#' @examples
#'
#' # update central log (usually at the end of a successful module run)
#' updateCentralLog(Module = "M01", pdf.flag = F)
#'
updateCentralLog <- function(Module, clog.file = "moduleLog.csv", log.path = if(exists("data.path")) data.path else "", user.id = if(exists("user")) user else "guest", run.notes = NA, pdf.flag = save.pdf){

  warning("Updating central log...\n")
  df.clog <- read.csv(paste0(data.path, clog.file))
  colnames(df.clog) <- rmvCSVprefix(colnames(df.clog))

  if (nrow(df.clog) > 0) if (max(df.clog$ID) != df.clog$ID[nrow(df.clog)]) warning("Last run ID is not assigned the highest valued ID. Check log.")

  if (nrow(df.clog) == 0){
    last.run <- 0
  } else {
    last.run <- df.clog$ID[nrow(df.clog)]
  }

  current.run <- last.run + 1

  current.clog <- data.frame(
    ID = current.run,
    Identifier = paste0(user, "_R", current.run),
    Module = Module,
    User = user.id,
    Date = format(Sys.time(), '%d %B, %Y'),
    Input = paste(which.data, collapse = ", "),
    Subset = paste(which.strata, collapse = ", "),
    HTML = T,
    PDF = pdf.flag,
    Notes = run.notes
  )

  if (nrow(df.clog) == 0){
    df.clog <-current.clog
  } else {
    df.clog <- bind_rows(df.clog, current.clog)
  }


  # override previous log
  write.csv(df.clog,
            file = paste0(data.path, clog.file),
            row.names = F)
  warning("Central log update succesfull!\n")
}


