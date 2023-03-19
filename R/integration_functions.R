#' Integrate scRNA-seq data using Scanorama
#'
#' This is a wrapper for Scanorama (https://github.com/brianhie/scanorama). Scanorama enables batch-correction and integration of heterogeneous scRNA-seq datasets.
#'
#' @param object Seurat object or list of seurat objects.
#' @param batch name of meta data field specifying batch groupings.
#' @param do.umap compute UMAP. Default is T.
#' @param nfeatures number of features to use for integration. Default is 2000.
#' @param vars2regress variable to regress out. Ignored if SCT assay is available in provided object.
#' @param dim number of dimensions to use for integration. Default is 50L.
#' @param verbose print progress. Default is T.
#' @name runScanorama
#' @return Seurat object.
#' @author Nicholas Mikolajewicz
#' @examples
#'  object <- runScanorama(object, batch = "sample")
#'
runScanorama <- function(object, batch,  do.umap = F,  nfeatures = 2000, vars2regress = "percent.mt", dim = 50L, verbose = T){

  miko_message("Checking python modules...")
  library(reticulate, quietly = T)
  if (!(py_module_available("scanorama"))){
    miko_message("Installing 'Scanorama' python module...")
    py_install("scanorama", pip = T)
  }
  scanorama <- import('scanorama')

  if (!("SCT" %in% names(object@assays))){
    object <- getMitoContent(object)
    object <- SCTransform(object, vars.to.regress = vars2regress, assay = "RNA", method = "glmGamPoi", verbose = F, vst.flavor = "v2")
  }

  merge.exists <- F
  if ("Seurat" %in% class(object)){

    object <- getMitoContent(object)
    so <- object;

    merge.exists <- T

    if (!("SCT" %in% names(so@assays))){
      miko_message("Running SCTransform...", verbose = verbose)
      so <- SCTransform(so, vars.to.regress = vars2regress, assay = "RNA", method = "glmGamPoi", verbose = F, vst.flavor = "v2")
    }


    object <- Seurat::SplitObject(so, split.by = batch)

  } else if ("list" %in% class(object)){

    object <- pbapply::pblapply(X = object, FUN = function(x){
      x <- getMitoContent(x)
      return(x)
    })

    for (i in 1:length(object)){
      if (!("SCT" %in% names(object[[i]]@assays))){
        miko_message(paste0("Running SCTransform for '", names(object)[i], "'..."), verbose = verbose)
        object[[i]] <- SCTransform(object[[i]], vars.to.regress = vars2regress, assay = "RNA", method = "glmGamPoi", verbose = F, vst.flavor = "v2")
      }
    }

  }

  miko_message("Selecting integration features...", verbose = verbose)
  so.features <- SelectIntegrationFeatures(object.list = object,
                                           assay = rep("SCT", length(object)),
                                           nfeatures = nfeatures,
                                           fvf.nfeatures = nfeatures)


  miko_message("Computing residuals for integration features...", verbose = verbose)
  for (i in 1:length(object)){
    object[[i]] <- GetResidual(object = object[[i]],features =  so.features, assay = "SCT", replace.value = F, verbose = F)
  }


  if (!merge.exists){
    miko_message("Merging seurat objects...", verbose = verbose)
    so <- tryCatch({
      merge(object[[1]], y = object[-1]) # 2-3x faster than alternative method
    }, error = function(e){
      return(mergeSeuratList(object))     # slower but robust
    }, silent = T)
  }

  datasets <- list()
  genes_list <- list()
  for (i in 1:length(object)) {
    datasets[[i]] <- as.matrix( Matrix::t((
      object[[i]]@assays$SCT@scale.data[so.features , ]
    )))
    genes_list[[i]] <- colnames(datasets[[i]])
  }

  # Integration.
  miko_message("Running Scanorama Integration...", verbose = verbose)
  suppressMessages({
    suppressWarnings({
      integrated.data <- scanorama$integrate(datasets, genes_list, dimred= dim)
    })
  })

  # df.embed <- rbind(integrated.data[[1]])
  df.embed <- NULL
  for (i in 1:length(integrated.data[[1]])){

    current_embed <- (integrated.data[[1]][[i]])

    if ("list" %in% class(object)){
      rownames(current_embed) <- colnames(object[[i]])
    }
    # rownames(current_embed) <- colnames(so.list[[i]])
    colnames(current_embed) <- paste0("PC_", seq(1, ncol(current_embed)))

    df.embed <- rbind(df.embed, current_embed)

  }

  if ("Seurat" %in% class(object)){
    rownames(current_embed) <- colnames(object)
  }

  so@assays[["integrated"]] <- CreateAssayObject(counts = so@assays[["SCT"]]@counts, min.cells = 0, min.features = 0)
  so@assays[["integrated"]]@key <- "integrated_"


  # set default assay
  DefaultAssay(so) <- "integrated"

  so[["pca"]]  <- CreateDimReducObject(embeddings = df.embed,
                                       key = "PC_",
                                       global  = F,
                                       loading = new(Class = "matrix"),
                                       assay = "integrated")
  # cluster.UMAP(so)

  if (do.umap){
    so <- RunUMAP(object = so, dims = 1:dim, verbose = verbose)
  }

  return(so)


}




#' Integrate scRNA-seq data using batch-balanced KNN (BBKNN)
#'
#' This is a wrapper for BBKNN (https://github.com/Teichlab/bbknn). BBKNN is a fast and intuitive batch effect removal tool. It generates a BBKNN cell x cell similarity graph that can be used for clustering and UMAP embedding.
#'
#' @param object Seurat object
#' @param batch name of meta data field specifying batch groupings.
#' @param reduction reduction used for batch correction. Default is "pca".
#' @param assay assay name. Default is DefaultAssay(object).
#' @param do.umap compute UMAP. Default is T.
#' @param verbose print progress. Default is T.
#' @name runBBKNN
#' @return Seurat object with bbKNN-corrected UMAP stored in "b" reduction slot (if do.umap = T).
#' @author Nicholas Mikolajewicz
#' @examples
#'  object <- runBBKNN(object, batch = "sample")
#'
runBBKNN <- function(object, batch, reduction = "pca", assay  = DefaultAssay(object), do.umap = T, verbose = T){

  require(reticulate, quietly = T)
  require(Matrix, quietly = T)
  miko_message("Checking inputs...", verbose = verbose)
  stopifnot("Seurat" %in% class(object))
  stopifnot(reduction %in% names(object@reductions))
  stopifnot(batch %in% colnames(object@meta.data))
  if (!py_module_available("anndata")) {
    ad.success <- F
    try({
      py_install("anndata")
      ad.success <- T
    }, silent = T)
    if (!(ad.success))
      try({
        py_install("anndata", pip = T)
      }, silent = T)
  }
  if (!py_module_available("scanpy")) {
    sp.success <- F
    try({
      py_install("scanpy")
      sp.success <- T
    }, silent = T)
    if (!(sp.success))
      try({
        py_install("scanpy", pip = T)
      }, silent = T)
  }
  if (!py_module_available("bbknn")) {
    bk.success <- F
    try({
      py_install("bbknn")
      bk.success <- T
    }, silent = T)
    if (!(bk.success))
      try({
        py_install("bbknn", pip = T)
      }, silent = T)
  }
  if (!py_module_available("anndata"))
    stop("Python module 'anndata' is not available.")
  if (!py_module_available("scanpy"))
    stop("Python module 'scanpy' is not available.")
  if (!py_module_available("bbknn"))
    stop("Python module 'bbknn' is not available.")
  if (ulength(object@meta.data[, batch]) > 1) {
    miko_message("Running BBKNN...", verbose = verbose)
    anndata = import("anndata", convert = FALSE)
    bbknn = import("bbknn", convert = FALSE)
    sc = import("scanpy", convert = FALSE)
    adata = anndata$AnnData(X = object@reductions[[reduction]]@cell.embeddings,#dtype = "float64", #dtype = "float32",
                            obs = object@meta.data[, batch])
    sc$tl$pca(adata)
    adata$obsm$X_pca = object@reductions[[reduction]]@cell.embeddings
    # adata$uns[["pca"]][["variance"]] = object@reductions[["pca"]]@stdev^2
    # adata$uns[["pca"]][["variance_ratio"]] = object@reductions[["pca"]]@stdev^2/sum(object@reductions[["pca"]]@stdev^2)
    bbknn$bbknn(adata, batch_key = 0, n_pcs  = 49L)
    miko_message("Getting BBKNN Graph...", verbose = verbose)
    snn.matrix <- py_to_r(adata$obsp[["connectivities"]])
    rownames(x = snn.matrix) <- colnames(x = snn.matrix) <- Cells(x = object)
    snn.matrix <- as(snn.matrix, "CsparseMatrix")
    snn.matrix <- as.Graph(x = snn.matrix)
    slot(object = snn.matrix, name = "assay.used") <- assay
    object[["bbknn"]] <- snn.matrix
    if (do.umap) {
      miko_message("Embedding UMAP...", verbose = verbose)
      sc$tl$umap(adata = adata)

      umap = py_to_r(adata$obsm[["X_umap"]])


      # umap <- a@cell.embeddings
      #
      # py_help(bbknn$bbknn)

      # nei <- py_to_r(adata$uns[["neighbors"]])
      # con <- py_to_r(adata$obsm[["X_pca"]])
      # pc = py_to_r(adata$uns[["neighbors"]])
      # pc = py_to_r(adata$uns[["neighbors"]])

      # umap <- RunUMAP(object@graphs[["bbknn"]])@cell.embeddings

      object@meta.data$b1 <- umap[, 1]
      object@meta.data$b2 <- umap[, 2]
      colnames(umap) <- c("b_1", "b_2")
      rownames(umap) <- colnames(object)
      miko_message("Storing results...", verbose = verbose)
      suppressMessages({
        suppressWarnings({
          object[["b"]] <- CreateDimReducObject(embeddings = as.matrix(umap),
                                                key = "B", global = T, loading = new(Class = "matrix"),
                                                assay = assay)
        })
      })
    }
    miko_message("Done.", verbose = verbose)
  }
  else {
    miko_message("Only one batch experiment provided. BBKNN was not run.")
  }
  return(object)

  #
  # require(reticulate, quietly = T)
  # require(Matrix, quietly = T)
  #
  #
  # miko_message("Checking inputs...", verbose = verbose)
  # stopifnot("Seurat" %in% class(object))
  # stopifnot(reduction %in% names(object@reductions))
  # stopifnot(batch %in% colnames(object@meta.data))
  #
  # # check that python modules are available
  # if (!py_module_available("anndata")){
  #   ad.success <- F
  #   try({py_install("anndata"); ad.success <- T}, silent = T)
  #   if (!(ad.success))  try({py_install("anndata", pip = T)}, silent = T)
  # }
  #
  # if (!py_module_available("scanpy")){
  #   sp.success <- F
  #   try({py_install("scanpy"); sp.success <- T}, silent = T)
  #   if (!(sp.success))  try({py_install("scanpy", pip = T)}, silent = T)
  # }
  #
  # if (!py_module_available("bbknn")){
  #   bk.success <- F
  #   try({py_install("bbknn"); bk.success <- T}, silent = T)
  #   if (!(bk.success))  try({py_install("bbknn", pip = T)}, silent = T)
  # }
  #
  # if (!py_module_available("anndata")) stop("Python module 'anndata' is not available.")
  # if (!py_module_available("scanpy")) stop("Python module 'scanpy' is not available.")
  # if (!py_module_available("bbknn")) stop("Python module 'bbknn' is not available.")
  #
  #
  # if (ulength(object@meta.data[ ,batch ]) > 1){
  #   # integrate data ########################################
  #   miko_message("Running BBKNN...", verbose = verbose)
  #   anndata = import("anndata",convert=FALSE)
  #   bbknn = import("bbknn", convert=FALSE)
  #   sc = import("scanpy",convert=FALSE)
  #
  #   adata = anndata$AnnData(X=object@reductions[[reduction]]@cell.embeddings, obs=object@meta.data[ ,batch ])
  #   sc$tl$pca(adata)
  #   adata$obsm$X_pca = object@reductions[[reduction]]@cell.embeddings
  #   bbknn$bbknn(adata,batch_key=0)
  #
  #   miko_message("Getting BBKNN Graph...", verbose = verbose)
  #   snn.matrix <- py_to_r(adata$obsp[["connectivities"]])
  #   rownames(x = snn.matrix) <- colnames(x = snn.matrix) <- Cells(x = object)
  #   snn.matrix <- as(snn.matrix, "CsparseMatrix")
  #   snn.matrix <- as.Graph(x = snn.matrix)
  #   slot(object = snn.matrix, name = "assay.used") <- assay
  #   object[["bbknn"]] <- snn.matrix
  #
  #   if (do.umap){
  #
  #     miko_message("Embedding UMAP...", verbose = verbose)
  #     sc$tl$umap(adata)
  #     umap = py_to_r(adata$obsm[["X_umap"]])
  #
  #     #########################################################
  #
  #     object@meta.data$b1 <- umap[ ,1]
  #     object@meta.data$b2 <- umap[ ,2]
  #
  #     colnames(umap) <- c("b_1", "b_2")
  #     rownames(umap) <- colnames(object)
  #
  #     miko_message("Storing results...", verbose = verbose)
  #     suppressMessages({
  #       suppressWarnings({
  #         object[["b"]]  <- CreateDimReducObject(embeddings = as.matrix(umap),
  #                                                key = "B",
  #                                                global  = T,
  #                                                loading = new(Class = "matrix"),
  #                                                assay = assay)
  #       })
  #     })
  #
  #   }
  #   miko_message("Done.", verbose = verbose)
  # } else {
  #   miko_message("Only one batch experiment provided. BBKNN was not run.")
  # }
  #
  #
  #
  # return(object)

}



#' scRNAseq integration wrapper
#'
#' scRNAseq normalization and integration wrapper. Given seurat object input, data are split, normalized and integrated.
#'
#' @param object Seurat object
#' @param split.by Meta data feature to split and integrate data.
#' @param min.cell Minimum number of cells permitted per object prior to integration.
#' @param k.anchor 	How many neighbors (k) to use when picking anchors.
#' @param k.weight Number of neighbors to consider when weighting anchors.
#' @param nfeatures Number of features to return (passed to SelectIntegrationFeatures)
#' @param split.prenorm Split data before (TRUE) or after (FALSE) normalization.
#' @param assay Assay to use for normalization.
#' @param variable.features.n Use this many features as variable features after ranking by residual variance; default is 3000.
#' @param verbose Print progress. Default is TRUE.
#' @param use.existing.sct If TRUE, existing SCT model is used. Default is FALSE (new SCT model is fit)
#' @param conserve.memory If set to TRUE the residual matrix for all genes is never created in full when running SCTransform; useful for large data sets, but will take longer to run; this will also set return.only.var.genes to TRUE; default is FALSE
#' @param vars.to.regress meta features to regress out. Default is "percent.mt". Set to NULL if unspecified.
#' @name miko_integrate
#' @seealso \code{\link{IntegrateData}}
#' @author Nicholas Mikolajewicz
#' @return Integrated seurat object
#' @examples
#'
miko_integrate <- function(object, split.by = "Barcode", min.cell = 50, k.anchor = 20, k.weight = 35, nfeatures = 3000, split.prenorm = F, assay = "RNA", variable.features.n = 3000, verbose = T, use.existing.sct = F, conserve.memory = F, vars.to.regress = "percent.mt"){

  verbose. <- verbose

  if(!("Seurat" %in% class(object))) stop("'object' must belong to Seurat class")

  # split data
  DefaultAssay(object) <- assay

  if (!is.null(vars.to.regress)){
    vars.to.regress <- vars.to.regress[vars.to.regress %in% colnames(object@meta.data)]
    if (length(vars.to.regress) == 0) vars.to.regress <- NULL
  }
  # if ("integrated" %in% names(object@assays)) object@assays$integrated <- NULL

  if (!use.existing.sct){

    if (split.prenorm){
      if (verbose.) miko_message("Spliting object...")
      object.list <- SplitObject(object, split.by = split.by)
      cell.counts <- (unlist(lapply(object.list, ncol)))
      object.list <- object.list[cell.counts>=min.cell]

      # renormalize
      if (verbose.) miko_message("Running SCTransform...", verbose = verbose.)

      object.list <- pbapply::pblapply(X = object.list, FUN = SCTransform, method = "glmGamPoi", vst.flavor = "v2", verbose = verbose., assay = assay, conserve.memory = conserve.memory,
                                       vars.to.regress = vars.to.regress, variable.features.rv.th = 1.3, variable.features.n = variable.features.n)

    } else {
      if (verbose.) miko_message("Running SCTransform...", verbose = verbose.)
      object <- SCTransform(object, method = "glmGamPoi", verbose = verbose., assay = assay, vst.flavor = "v2",
                            vars.to.regress = vars.to.regress, variable.features.rv.th = 1.3, variable.features.n = variable.features.n, conserve.memory = conserve.memory)
      if (verbose.)  miko_message("Spliting object...", verbose = verbose.)
      object.list <- SplitObject(object, split.by = split.by)
    }

  } else {
    object.list <- SplitObject(object, split.by = split.by)
  }

  # select features
  if (verbose.) miko_message("Selecting features for integration...", verbose = verbose.)
  cell.counts <- (unlist(lapply(object.list, ncol)))
  object.list <- object.list[cell.counts>=min.cell]
  object.features <- SelectIntegrationFeatures(object.list = object.list, nfeatures = nfeatures)
  object.list <- PrepSCTIntegration(object.list = object.list, anchor.features = object.features)

  # run PCA

  if (verbose.) miko_message("Running PCA...", verbose = verbose.)
  min.dim <- min(unlist(lapply(object.list, ncol)))
  if (min.dim > 50){pca.dim <- 50} else {pca.dim <- min.dim-1}
  object.list <- pbapply::pblapply(X = object.list, FUN = RunPCA, features = object.features, verbose = verbose., npcs = pca.dim)

  # find integration anchors
  if (verbose.) miko_message("Finding integration anchors...", verbose = verbose.)
  object.anchors <- FindIntegrationAnchors(object.list = object.list, normalization.method = "SCT", k.filter = 70,
                                           anchor.features = object.features, dims = 1:30, reduction = "rpca", k.anchor = k.anchor)

  # integrate data
  if (k.weight > min.dim) k.weight <- min.dim
  if (verbose.) miko_message("Integrating data...", verbose = verbose.)
  object2 <- IntegrateData(anchorset = object.anchors, normalization.method = "SCT", dims = 1:30, k.weight = k.weight)

  # dimensional reduction
  if (verbose.) miko_message("Running PCA on integrated object...", verbose = verbose.)
  var.threshold <- 0.9
  object2 <- RunPCA(object2, verbose = verbose.)
  df.object.var <- propVarPCA(object2)
  object.n.dim <- df.object.var$pc.id [min(which(df.object.var$pc.cum_sum > var.threshold))]
  miko_message("Embedding UMAP...", verbose = verbose)
  object2 <- RunUMAP(object2, reduction = "pca", dims = 1:object.n.dim, return.model = TRUE)

  # cluster data
  if (verbose.) miko_message("Running clustering algorithms...", verbose = verbose.)
  object2 <- FindNeighbors(object2, dims = 1:object.n.dim, verbose = verbose.)
  object2 <- FindClusters(object2, verbose = verbose.)

  # correct SCT counts
  miko_message("Adjusting counts to fixed sequencing depth...", verbose = verbose.)
  try({
    object2 <- PrepSCTFindMarkers(object2)
  }, silent = T)


  if (verbose.) miko_message("Complete!")

  return(object2)

}

