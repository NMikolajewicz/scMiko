

#' Initiate analysis log
#'
#' Create analysis log data.frame template
#'
#' @param module.name Character specifying module name
#' @name initiateLog
#' @seealso \code{\link{addLogEntry}}
#' @return data.frame
#' @examples
#'
#' # intiate data log
#' df.log <- initiateLog("log name")
#'
#' # add new entry to log
#' df.log <- addLogEntry("Query File (.Rdata)", (input.file), df.log, "input.file")
#'
initiateLog <- function (module.name = ""){

  # Module
  df.log <- data.frame()
  df.log[nrow(df.log)+1, 1] <- as.character("Module")
  df.log[nrow(df.log), 2] <- as.character("")
  df.log[nrow(df.log), 3] <- as.character(module.name)
  colnames(df.log) <- c("Description", "Variable Name", "Value")

  # Date
  df.log[nrow(df.log)+1, 1] <- as.character("Date")
  df.log[nrow(df.log), 2] <- as.character("Sys.time()")
  df.log[nrow(df.log), 3] <- as.character(Sys.time())

  return(df.log)
}


#' Add entry to analysis log
#'
#' Add entry to analysis log
#'
#' @param entry.name name of entry. A character.
#' @param entry entry.
#' @param df.log existing analysis log to add entry to. A data.frame.
#' @param entry.variable.name name of variable storing entry (optional). A character.
#' @name addLogEntry
#' @seealso \code{\link{initiateLog}}
#' @author Nicholas Mikolajewicz
#' @return datalog data.frame
#' @examples
#'
#' # intiate data log
#' df.log <- initiateLog("log name")
#'
#' # add new entry to log
#' df.log <- addLogEntry("Query File (.Rdata)", (input.file), df.log, "input.file")
#'
addLogEntry <- function(entry.name, entry, df.log, entry.variable.name = ""){

  df.log[nrow(df.log)+1, 1] <- as.character(entry.name)
  df.log[nrow(df.log), 2] <- as.character(entry.variable.name)
  df.log[nrow(df.log), 3] <- paste(entry, collapse=", ")

  return(df.log)
}




#' prep Seurat (Extended adaptation of prepSeurat)
#'
#' Extended version of original prep Seurat function. In addition to running standard seurat object checks, prepSeurat2 ensures correct species is specified, sets the cluster resolution and subsets and downsamples data. If data are subset, then expression values are normalized and (optionally) scaled again.
#'
#' @param object Seurat objects
#' @param e2s ensemble to gene symbol mapping vector. Is a named character vector where names are ENSEMBL and entries are SYMBOLs.
#' @param species Species. One of "Mm" or "Hs".
#' @param resolution cluster resolution. Numeric [0,Inf] that specifies resolution for data clustering. If requested resolution exists, no new clustering is performed.
#' @param subset.data Data.frame specifying how to subset seurat object. Data.frame must contain two columns: 'field' and 'subgroups'. 'field' column specifies which meta.data field to subset on, and 'subgroups' column specifies which subgroup to include within specified field. See scMiko::subsetSeurat() for details.
#' @param subsample Numeric [0,1] specifying what fraction of cells to include for analysis. Default is 1. See scMiko::downsampleSeurat() for details.
#' @param terms2drop Reduce memory footprint of seurat object by omitting terms that will not be used for current analysis. Supported terms for omission include: "pca", "umap", "ica", "tsne", "nmf", "corr", "gsva", "deg", "counts", "data", "scale", "rna", "sct", "integrated", "graphs", "integration.anchors".
#' @param rmv.pattern Provided as input into scMiko::clearGlobalEnv(pattern = rmv.pattern). Character specifying name of variables to remove from global environment. Useful if object is large.
#' @param reprocess.n.var Number of variable genes to use if data is reprocessed. Default is 3000. Note that if integrated assay is available, variable features are first identified for reprocessed data set, and subsequently merged with the variable features present in the integrated assay, thus allows for potentially more variable features than specified by this parameter.
#' @param neighbors.reprocessed Specifies whether to compute new neighborhood graph if graphs are missing. Note that if data are subset, graphs are inherently removed. If downstream clustering is anticipated, set as TRUE. Default is FALSE.
#' @param use.integrated If TRUE, sets default assay to "integrated" if present within seurat object.
#' @param scale.reprocessed if reprocessing data (i.e., normalizing), specify whether scaling should also be performed. Default is FALSE.
#' @param keep.default.assay.only Specify whether to omit assays that are not default. Default is FALSE.
#' @param coerce.assay.used.to.default Specify whether to coerce assay used to default assay. Necessary if omitting assays (e.g., integrated). Default is TRUE.
#' @param barcode.recode List specifying how to recode barcodes. Default is NULL. See recodeBarcode() for details.
#' @name prepSeurat2
#' @author Nicholas Mikolajewicz
#' @return list containing prepped Seurat object, default assay, and number of cells in seurat object.
#'
prepSeurat2 <- function (object, e2s = NULL, species = NULL, resolution= NULL, subset.data = NULL, subsample = 1,
                         terms2drop = NULL, rmv.pattern = NULL, reprocess.n.var = 3000, neighbors.reprocessed = F, scale.reprocessed = F, use.integrated = T,
                         keep.default.assay.only = F, coerce.assay.used.to.default = T, barcode.recode = NULL){

  miko_message("Checking seurat object...")
  # assertion
  if (class(object) != "Seurat") stop("input must be Seurat Object")


  if (use.integrated){
    if ("integrated" %in% names(object@assays)){
      DefaultAssay(object) <- "integrated"
    }
  }


  # remove object from global environment #####################################
  if (!is.null(rmv.pattern)) {
    scMiko::clearGlobalEnv(rmv.pattern)
    invisible({gc()})
  }

  # get species ################################################################
  if (is.null(species)){
    species <- detectSpecies(object)
  }

  # general object handling ####################################################

  object <- fixBarcodeLabel(object)

  try({
    object <- UpdateSeuratObject(object) # required after Seurat 3.1.2 update
  }, silent = T)

  try({
  object <- updateDimNames(object) # required after Seurat 3.1.2 update
  }, silent = T)

  # species
  if (!("Organism" %in% colnames(object@meta.data))){
    object@meta.data[["Organism"]] <- detectSpecies(object)
  }

  if (species != unique(object@meta.data[["Organism"]])) {
    if (length(unique(object@meta.data[["Organism"]])) == 1) {
      species <- unique(object@meta.data[["Organism"]])
      miko_message("Incorrect input species was provided, and was updated to reflect what was available in seurat object.")
    } else if (length(unique(object@meta.data[["Organism"]])) > 1)  {
      miko_message("Input species could not be verfied and is being used as-is.")
    }
  }

  # set resolution #############################################################
  if (!is.null(resolution) && is.numeric(resolution)){
    miko_message("Setting cluster resolution...")
    object <-   setResolution(object, resolution = resolution)
    invisible({gc()})
  }

  # hold graphs ################################################################
  if (neighbors.reprocessed) graph.holder <- object@graphs

  # subsample ##################################################################
  n.presubsample <- ncol(object)
  if (subsample < 1 && is.numeric(subsample)){
    miko_message("Downsampling data...")
    object <- tryCatch({
      object = subset(object, cells = sample(Cells(object), round(ncol(object)*subsample)))
    }, error = function(e){
      print(e);
      return(object)
    }
    )
  }
  n.postsubsample <- ncol(object)

  n.presubset <- ncol(object)
  if (!is.null(subset.data) && ("data.frame" %in% class(subset.data)) && !is.na(unique(subset.data[,1]))) {
    miko_message("Subsetting data...")
    object <- scMiko::subsetSeurat(object,subset.data)
    invisible({gc()})
  }
  n.postsubset <- ncol(object)

  # recompute neighbors if they were removed####################################
  if (neighbors.reprocessed){

    if ((n.presubset < n.postsubset) | (n.presubsample < n.postsubsample) | length(object@graphs) == 0){

      # try to use original neighbor graph. If error is encountered, compute new graph.
      object <- tryCatch({
        which.cells.remain <- colnames(object)
        for (i in 1:length(graph.holder)){
          graph.name <- names(graph.holder)[i]
          current.graph <- graph.holder[[graph.name]]
          graph.class <- class(current.graph)
          current.graph <- current.graph[colnames(current.graph) %in% which.cells.remain,
                                         colnames(current.graph) %in% which.cells.remain]
          if (canCoerce(current.graph, graph.class)) {
            current.graph <- as(current.graph, graph.class)
            current.graph@assay.used <- DefaultAssay(object)
          }
          graph.holder[[graph.name]] <- current.graph
        }

        object@graphs <- graph.holder
        rm(graph.holder)
        object <- suppressMessages({UpdateSeuratObject(object)})
      }, error = function(e){

        if ("pca" %in% names(object@reductions)){
          pca.prop <- propVarPCA(object)
        } else {
          if (length(object@assays[[DefaultAssay(object)]]@scale.data) == 0){
            object <- scNormScale(object, method = "NFS", assay = DefaultAssay(object), enable.parallelization = F)
            object <- RunPCA(object, verbose = F)
          }
          pca.prop <- propVarPCA(object)
        }


        target.pc <- max(pca.prop$pc.id[pca.prop$pc.cum_sum<0.9])+1
        object <- FindNeighbors(object, verbose = T, reduction = "pca", dims = 1:target.pc)
        object <- suppressMessages({UpdateSeuratObject(object)})
        return(object)
      })

    }
    rm(graph.holder); invisible({gc()})

  }

  # trim down seurat object ####################################################
  if (!is.null(terms2drop)){
    supported.terms <- c("pca", "umap", "ica", "tsne", "nmf", "corr", "gsva", "deg",
                         "counts", "data", "scale", "rna", "sct", "integrated", "graphs", "integration.anchors")
    terms2drop <- terms2drop[terms2drop %in% supported.terms]

    # empty sparse matrix
    null.dgcmat <- as(matrix(NA), "dgCMatrix")

    if (("pca" %in% terms2drop) & ("pca" %in% names(object@reductions))) try({object@reductions[["pca"]] <- NULL}, silent = T)
    if (("umap" %in% terms2drop) & ("umap" %in% names(object@reductions))) try({object@reductions[["umap"]] <- NULL}, silent = T)
    if (("tsne" %in% terms2drop) & ("tsne" %in% names(object@reductions))) try({object@reductions[["tsne"]] <- NULL}, silent = T)
    if (("ica" %in% terms2drop) & ("ica" %in% names(object@reductions))) try({object@reductions[["ica"]] <- NULL}, silent = T)
    if (("nmf" %in% terms2drop) & ("nmf" %in% names(object@misc))) try({object@misc[["nmf"]] <- NULL}, silent = T)
    if (("corr" %in% terms2drop) & ("similarity.scale" %in% names(object@misc))) try({object@misc[["similarity.scale"]] <- NULL}, silent = T)
    if (("gsva" %in% terms2drop) & ("gsva" %in% names(object@misc))) try({object@misc[["gsva"]] <- NULL}, silent = T)
    if (("deg" %in% terms2drop) & ("deg" %in% names(object@misc))) try({object@misc[["deg"]] <- NULL}, silent = T)
    if ("counts" %in% terms2drop) {
      try({object@assays[["RNA"]]@counts <- null.dgcmat}, silent = T)
      try({object@assays[["SCT"]]@counts <- null.dgcmat}, silent = T)
      try({object@assays[["integrated"]]@counts <- null.dgcmat}, silent = T)
    }
    if ("data" %in% terms2drop) {
      try({object@assays[["RNA"]]@data <- null.dgcmat}, silent = T)
      try({object@assays[["SCT"]]@data <- null.dgcmat}, silent = T)
      try({object@assays[["integrated"]]@data <- null.dgcmat}, silent = T)
    }
    if ("scale" %in% terms2drop) {
      try({object@assays[["RNA"]]@scale.data <- matrix(NA)}, silent = T)
      try({object@assays[["SCT"]]@scale.data <- matrix(NA)}, silent = T)
      try({object@assays[["integrated"]]@scale.data <- matrix(NA)}, silent = T)
    }
    if ("rna" %in% terms2drop) try({object@assays[["RNA"]] <- NULL}, silent = T)
    if ("sct" %in% terms2drop) try({object@assays[["SCT"]] <- NULL}, silent = T)
    if ("integrated" %in% terms2drop) try({object@assays[["integrated"]] <- NULL}, silent = T)
    if ("graphs" %in% terms2drop) try({object@graphs <- NULL}, silent = T)
    if ("integration.anchors" %in% terms2drop) try({object@tools[["Integration"]]@anchors <- NULL}, silent = T)

    object <- suppressMessages({UpdateSeuratObject(object)})
  }

  # convert ENSEBLE to GENE names in Seurat object #############################

  if (!is.null(e2s)){

    gene.rep <-  checkGeneRep(e2s, as.vector(rownames(object)))

    if (gene.rep == "ensembl"){
      miko_message("Converting ENSEMBL to SYMBOL...")
      # filter species-specific genes
      if (species == "Mm"){
        object <- subset(object, features = unique(rownames(object)[grepl("MUSG", rownames(object))]))
      } else if (species == "Hs"){
        object <- subset(object, features =  unique(rownames(object)[!(grepl("MUSG", rownames(object)))]))
      }

      # convert ensembl to symbol
      object <- ens2sym.so(object, gNames.list = e2s)
      gene.rep <-  checkGeneRep(e2s, as.vector(rownames(object)))
      object <- suppressMessages({UpdateSeuratObject(object)})

      invisible({gc()})
    }

  }

  # remove duplicate genes #####################################################
  object <- rmDuplicateGenes(object)
  invisible({gc()})

  # integrated data object handling ############################################
  all.assays <- names(object@assays)
  all.commands <- names(object@commands)
  n.var.genes <- reprocess.n.var
  data.reprocessed <- F


  if (object@version < 4){

    if (("integrated" %in% all.assays) & ("NormalizeData.RNA" %in% all.commands) & ("ScaleData.RNA" %in% all.commands)){
      miko_message("Ensuring correct assays are set...")
      if (DefaultAssay(object) != "RNA") {
        miko_message("Setting default assay to 'RNA'...")
        DefaultAssay(object) <- "RNA"
      }
      miko_message("Finding variable features...")
      object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = n.var.genes)
      if (length(object@assays[["integrated"]]@var.features) > 0){
        object@assays[["RNA"]]@var.features <- unique(c(object@assays[["RNA"]]@var.features, object@assays[["integrated"]]@var.features))
      }

    } else if (("integrated" %in% all.assays) & (!("NormalizeData.RNA" %in% all.commands) | !("ScaleData.RNA" %in% all.commands))){
      miko_message("Setting default assay to 'RNA'...")
      DefaultAssay(object) <- "RNA"

      if (!("NormalizeData.RNA" %in% all.commands)){
        if ("counts" %in% terms2drop) stop("Cannot normalize and scale data because counts were omitted from Seurat Object. Remove 'counts' from terms2drop and try again.")
        miko_message("Normalizing data...")
        object <-NormalizeData(object, verbose = FALSE)
        invisible({gc()})
      }
      if (!("ScaleData.RNA" %in% all.commands)){
        if (scale.reprocessed & ("integrated" %in% all.assays)){
          miko_message(paste0("Scaling ", length(rownames(object)), " genes in data..."))
          object <- ScaleData(object, verbose = FALSE, features = rownames(object))
          invisible({gc()})
        }
      }

      miko_message("Finding variable features...")
      object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = n.var.genes)
      if (length(object@assays[["integrated"]]@var.features) > 0){
        object@assays[["RNA"]]@var.features <- unique(c(object@assays[["RNA"]]@var.features, object@assays[["integrated"]]@var.features))
      }
      invisible({gc()})

      data.reprocessed <- T
    }


  } else {

    if (("integrated" %in% all.assays) & ("SCT" %in% all.assays) ){
      miko_message("Setting default assay to 'SCT'...")
      DefaultAssay(object) <- "SCT"
      try({
        object@assays[["SCT"]]@var.features <- object@assays[["integrated"]]@var.features
      }, silent = T)

    }

  }


  # rescale if data was subset #################################################
  if ((n.presubset != n.postsubset) && (data.reprocessed == F)){
    try({
      if (length(object@assays[[DefaultAssay(object)]]@var.features) > 0) nVar <- length(object@assays[[DefaultAssay(object)]]@var.features) else nVar <- reprocess.n.var
    }, silent = T)
    if (!exists("nVar")) nVar <- reprocess.n.var
    if ("counts" %in% terms2drop) stop("Cannot normalize and scale data because counts were omitted from Seurat Object. Remove 'counts' from terms2drop and try again.")
    miko_message("Normalizing subset data...")
    object <-NormalizeData(object, verbose = FALSE)
    invisible({gc()})
    if (scale.reprocessed & ("integrated" %in% all.assays)){
      miko_message(paste0("Scaling ", length(rownames(object)), " genes in subset data..."))
      object <- ScaleData(object, verbose = FALSE, features = rownames(object))
      invisible({gc()})
    }
    miko_message(paste0("Finding top ", nVar, " variable genes in subset data..."))
    object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = nVar)
    invisible({gc()})

    data.reprocessed <- T
  }

  # Coerce 'used assays' to default ############################################
  if (coerce.assay.used.to.default){
    try({object@reductions[["pca"]]@assay.used <- DefaultAssay(object)}, silent = T)
    try({object@reductions[["tsne"]]@assay.used <- DefaultAssay(object)}, silent = T)
    try({object@reductions[["umap"]]@assay.used <- DefaultAssay(object)}, silent = T)
    try({object@reductions[["ica"]]@assay.used <- DefaultAssay(object)}, silent = T)
  }


  # Omit assays ################################################################
  if (keep.default.assay.only){

    which.default <- DefaultAssay(object)
    all.assays <- names(object@assays)
    which.omit <- all.assays[!(all.assays %in% which.default)]

    if (length(which.omit) > 0){
      miko_message(paste0("Omitting the following assays: ", paste(which.omit, collapse = ", "), "..."))
      for (i in 1:length(which.omit)){
        try({object@assays[[which.omit[i]]] <- NULL}, silent = T)
      }
    }
    invisible({gc()})
  }

  # Recode Barcodes ############################################################
  if (!is.null(barcode.recode)){
    object <- scMiko::recodeBarcode(object, barcode.recode)
  }

  # enforce ordered clusters  ##################################################
  try({
    if (all(Idents(object) == object@meta.data$seurat_clusters) ){
      Idents(object) <- orderedFactor(Idents(object))
      object@meta.data$seurat_clusters <- orderedFactor(object@meta.data$seurat_clusters)
    } else {
      object@meta.data$seurat_clusters <- orderedFactor(object@meta.data$seurat_clusters)
    }
  }, silent = T)


  # generate UMAP if missing ###################################################
  if (!("umap" %in% names(object@reductions))){
    object <- RunUMAP(object, dims = 1:30)
  }

  # Return results #############################################################

  return(list(
    so = object,
    assay = DefaultAssay(object),
    n.cell = ncol(object),
    rescaled = data.reprocessed,
    species = species
  ))
}

#' prep Seurat
#'
#' Preprocess using fixBarcodeLabel(), UpdateSeuratObject() and updateDimNames() Functions.
#'
#' @param object Seurat objects
#' @name prepSeurat
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#'
prepSeurat <- function (object){

  if (class(object) != "Seurat") stop("input must be Seurat Object")

  object <- fixBarcodeLabel(object)

  try({
    object <- UpdateSeuratObject(object) # required after Seurat 3.1.2 update
  }, silent = T)

  object <- updateDimNames(object) # required after Seurat 3.1.2 update

  return(object)
}


#' Remove variables from global environment
#'
#' Remove variables in global environment that match a given pattern.
#'
#' @param pattern Character specifying name of variable(s) to remove from global environment.
#' @param exact.match If TRUE (default), exact pattern matching is enforced. If FALSE, partial pattern matching is facilitated with grep().
#' @name clearGlobalEnv
#' @author Nicholas Mikolajewicz
#' @return NULL
#'
clearGlobalEnv <- function(pattern, exact.match = T){

  pattern <- as.character(pattern)

  if (!exact.match){
    which.match <- ls(envir=globalenv())[grep(pattern, ls(envir=globalenv()))]
    rm(list = which.match, envir = globalenv())
    warning(paste0("The following variable were removed from the global enviroment: ", paste(which.match, collapse = ", "), "\n"))
  } else {
    which.match <- ls(envir=globalenv())[which(ls(envir=globalenv()) %in% pattern)]
    rm(list = which.match, envir = globalenv())
    warning(paste0("The following variable were removed from the global enviroment: ", paste(which.match, collapse = ", "), "\n"))
  }
  invisible({gc()})

}


#' Recode (i.e., relabel) metadata in Seurat object
#'
#' Recode (i.e., relabel) metadata in Seurat object. If old label is specified as NA, treated as wildcard that captures all remaining barcodes that have not been relabeled.
#'
#' @param object Seurat Object
#' @param barcode.recode Named list where names are new labels, entries are old labels (can be vector of characters). E.g., barcode.recode <- list(in.vitro = c("invitro"), in.vivo = NA)
#' @param old.field.name name of field in which meta dat for recoding is stored. Default is "Barcode".
#' @param new.field.name name of new field in meta data where recoded data will be stored. Default is "Barcode".
#' @name recodeBarcode
#' @author Nicholas Mikolajewicz
#' @return Seurat Object with relabeled metadata
#'
recodeBarcode <- function(object, barcode.recode, old.field.name = "Barcode", new.field.name = "Barcode"){

  # input options ##############################################################
  # named list where names are new barcodes, entries are old barcodes (can be vector of characters)
  # e.g., barcode.recode <- list(in.vitro = c("invitro"), in.vivo = NA)
  #
  # Note that if NA is specified as old name, this is equivalent to getting all remainder cells that aren't specified by other old barcodes.
  #
  ##############################################################################

  stopifnot("'object' is not a Seurat object" = ("Seurat" %in% class(object)))
  stopifnot("'old.field.name' is not found in meta data of Seurat object" = (old.field.name %in% colnames(object@meta.data)))

  if (exists("barcode.recode") && (((!is.null(barcode.recode)) | (!is.na(barcode.recode))) & is.list(barcode.recode))){
    barcodes.old <- as.character(object@meta.data[[old.field.name]])
    barcodes.new <- barcodes.old

    recode.log <- c()
    recode.rest <- F
    for (i in 1:length(barcode.recode)){

      if (all(!is.na(barcode.recode[[i]]))){
        barcode.pattern <- paste(barcode.recode[[i]], collapse = "|")
        which.match <- grep(barcode.pattern, barcodes.old)
        recode.log <- c(recode.log, which.match)
        barcodes.new[which.match] <- names(barcode.recode)[i]
      } else {
        recode.rest <- T
        rest.ind <- i
      }
    }

    if (recode.rest){
      all.ind <- c(1:length(barcodes.old))
      which.rest <- all.ind[!(all.ind %in% recode.log)]
      barcodes.new[which.rest] <- names(barcode.recode)[rest.ind]
    }


    object@meta.data[[new.field.name]] <- barcodes.new
  }

  return(object)

}


