
#' Calculate Gini marker specificity
#'
#' Calculate Gini marker specificity.
#'
#' @param object Seurat object
#' @param features feature to run gini specificy scoring on.
#' @param group.by groups used for feature specificity assignment. E.g., Each gene is assigned to a single group in which it exhibits maximal specificity. Default is 'seurat_clusters'.
#' @param assay Assay to run  gini specificity scoring on. Default is DefaultAssay(object).
#' @param min.pct minimal expression fraction per group. Gini scoring performs poorly on genes expressed at extremely low levels. Default (recommended) is 0.1.
#' @param verbose print progress. Default is T.
#' @name findGiniMarkers
#' @seealso \code{\link{fdrtool}} for p values calculations, \code{\link{ineq}} for gini index calculations.
#' @author Nicholas Mikolajewicz
#' @return data.frame with feature-specific gini scores.
#' @examples
#'
findGiniMarkers <- function(object, features = NULL, group.by = "seurat_clusters", assay = NULL, min.pct = 0.1, verbose = T){

  miko_message("Running Gini specificity analysis...", verbose = T)


  getExpr <- function (object, assay = NULL, features,  group.by = NULL, verbose = T){

    idents <-  NULL
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    feature.groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
      feature.groups <- unlist(x = sapply(X = 1:length(features),
                                          FUN = function(x) {
                                            return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                          }))
      if (any(is.na(x = feature.groups))) {
        warning("Some feature groups are unnamed.", call. = FALSE,
                immediate. = TRUE)
      }
      features <- unlist(x = features)
      names(x = feature.groups) <- features
    }
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
    data.features <- FetchData(object = object, vars = features,
                               cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
      Idents(object = object)[cells, drop = TRUE]
    } else {
      object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    # [cells, drop = TRUE]
    # a <- object[[group.by]]
    if (!is.factor(x = data.features$id)) {
      data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)

    if (verbose){
      data.plot <- pbapply::pblapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
                                  1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
          return(mean(x = expm1(x = x)))
        })
        return(list(avg.exp = avg.exp))
      })
    } else {
      data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
        data.use <- data.features[data.features$id == ident,
                                  1:(ncol(x = data.features) - 1), drop = FALSE]
        avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
          return(mean(x = expm1(x = x)))
        })
        return(list(avg.exp = avg.exp))
      })
    }

    names(x = data.plot) <- unique(x = data.features$id)
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)

    if (!is.null(x = id.levels)) {
      data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }

    return(data.plot)

  }

  # if (length(features) == 0) stop("'features' are not specified'")
  if (is.null(assay))  assay <- DefaultAssay(object = object)
  if (is.null(features)) {
    features <- getExpressedGenes(object = object, min.pct = min.pct, group = group.by)
  } else {
    features <- features[features %in% rownames(object)]
  }

  df.dot <- getExpr(object = object, features = features, group.by = group.by, verbose = verbose, assay = assay)
  df.dot.sum <- df.dot %>%
    dplyr::group_by(features.plot) %>%
    dplyr::summarize(gini.expr = ineq::ineq(x = rescale(avg.exp), type = "Gini"),
                     group = id[which.max(avg.exp)])
  df.dot.sum$ngini <- df.dot.sum$gini.expr/((ulength(object@meta.data[ ,group.by]) -1)/ulength(object@meta.data[ ,group.by]))
  df.dot.sum$gini.expr <- signif(df.dot.sum$gini.expr , 3)
  df.dot.sum$ngini <- signif(df.dot.sum$ngini , 3)
  colnames(df.dot.sum) <- c("feature", "gini","group", "ngini")

  gini_p <- fdrtool::fdrtool((df.dot.sum$ngini - median(df.dot.sum$ngini))/sd(df.dot.sum$ngini), statistic=c("normal"),
                             plot=F, color.figure=F, verbose=F,
                             cutoff.method=c("fndr"),
                             pct0=0.75)
  df.dot.sum$p <- signif(gini_p[["pval"]], 3)

  df.dot.sum <- df.dot.sum %>% dplyr::arrange(-ngini)

  return(df.dot.sum)

}




#' Calculate feature dependency index
#'
#' Calculate feature dependency index (CDI).
#'
#' @param object Seurat object
#' @param features.x feature or meta feature. CDI between features.x and features.y are computed.
#' @param features.y feature or meta feature. CDI between features.x and features.y are computed.
#' @param ncell.subset max number of cells to run analysis on. Default is 5000. Computationally intensive for larger datasets.
#' @param geosketch.subset Use GeoSketch method to subsample scRNA-seq data while preserving rare cell states (https://doi.org/10.1016/j.cels.2019.05.003). Logical, T or F (Default F). Recommended if cell type representation is imbalanced.
#' @param assay Assay to run CDI scoring on. Default is DefaultAssay(object).
#' @param slot slot to run CDI scoring on. Default is data.
#' @param n.workers number of workers for parallel implementation. Default is 1 (no parallel).
#' @param verbose print progress. Default is T.
#' @name findCDIMarkers
#' @seealso \code{\link{binom.test}}
#' @author Nicholas Mikolajewicz
#' @return data.frame with CDI scores.
#' @examples
#'
findCDIMarkers <- function(object, features.x = NULL, features.y = rownames(object),
                           ncell.subset = 5000, geosketch.subset = F, assay =  DefaultAssay(object), slot = "data", n.workers = 1, verbose = T){


  miko_message("Running CDI specificity analysis...", verbose = verbose)
  if (verbose){
    mylapply <- pbapply::pblapply
    myapply <- pbapply::pbapply
    mysapply <- pbapply::pbsapply
  } else {
    mylapply <- lapply
    myapply <- apply
    mysapply <- sapply
  }

  if (!("Seurat" %in% class(object))) stop("'object' is not Seurat object")
  # get expression matrix
  if (is.null(ncell.subset)) ncell.subset <- ncol(object)
  emat <- getExpressionMatrix(object, which.assay =DefaultAssay(object), which.data =slot )
  if (ncell.subset >= ncol(object)){
    ncell.subset <- ncol(object)
  }

  # check available features
  if (is.null(features.x)) features.x <- rownames(object)
  if (is.null(features.y)) features.y <- rownames(object)
  x_av <- features.x[features.x %in% rownames(emat)]
  x_missing <- features.x[!(features.x %in% x_av)]
  y_av <- features.y[features.y %in% rownames(emat)]
  y_missing <- features.y[!(features.y %in% y_av)]

  # convert meta features to incidence matrix and append to expression matrix
  x_meta_names <- c()
  if (length(x_missing) > 0){
    x_meta <- x_missing[x_missing %in% colnames(object@meta.data)]
    for (i in 1:length(x_meta)){
      meta.list <- group2list(object = object, group = x_meta[i])
      names(meta.list) <- paste0(x_meta[i], "_", names(meta.list))
      meta.list <- meta.list[order(names(meta.list))]
      meta.mat <- matrix(nrow = length(meta.list), ncol = ncol(object), data = 0)
      rownames(meta.mat) <- names(meta.list)
      colnames(meta.mat) <- colnames(object)
      for (j in 1:length(meta.list)){
        meta.mat[names(meta.list)[j],colnames(meta.mat) %in% meta.list[[j]]] <- 1
      }
      x_av <- c(x_av, names(meta.list))
      x_meta_names <- c(x_meta_names, names(meta.list))
      emat <- rbind(emat, meta.mat)
    }
  }

  y_meta_names <- c()
  if (length(y_missing) > 0){
    y_meta <- y_missing[y_missing %in% colnames(object@meta.data)]
    for (i in 1:length(y_meta)){
      meta.list <- group2list(object = object, group = y_meta[i])
      names(meta.list) <- paste0(y_meta[i], "_", names(meta.list))
      meta.list <- meta.list[order(names(meta.list))]
      meta.mat <- matrix(nrow = length(meta.list), ncol = ncol(object), data = 0)
      rownames(meta.mat) <- names(meta.list)
      colnames(meta.mat) <- colnames(object)
      for (j in 1:length(meta.list)){
        meta.mat[names(meta.list)[j],colnames(meta.mat) %in% meta.list[[j]]] <- 1
      }
      y_av <- c(y_av, names(meta.list))
      y_meta_names <- c(x_meta_names, names(meta.list))
      emat <- rbind(emat, meta.mat)
    }
  }

  # subsample
  if (ncell.subset < ncol(object)){
    set.seed(1023)

    if (geosketch.subset){
      library(reticulate, quietly = T)

      if (!py_module_available("geosketch")){
        ad.success <- F
        try({py_install("geosketch"); ad.success <- T}, silent = T)
        if (!(ad.success))  try({py_install("geosketch", pip = T)}, silent = T)
      }

      geosketch = import("geosketch",convert=FALSE)
      X_dimred <- object@reductions[["pca"]]@cell.embeddings

      miko_message(paste0("Subsampling ", ncell.subset, " cells using geometric sketching..."), verbose = verbose)
      sketch_index = geosketch$gs(X_dimred, N = as.integer(ncell.subset), replace= F, one_indexed = T)

      sub_ind <- unlist(py_to_r(sketch_index))
      emat <- emat[ ,sub_ind]

    } else {
      miko_message(paste0("Subsampling ", ncell.subset, " cells using uniform sampling..."), verbose = verbose)
      emat <- emat[ ,sample(seq(1, ncol(emat)), ncell.subset)]
    }
  }

  all.av <- unique(c(x_av, y_av))
  emat <- emat[rownames(emat) %in% all.av, ]
  emat <- (emat > 0)

  # which cells express
  miko_message("Identifying cells with non-zero gene expression...", verbose = verbose)
  # which.cells2 <- myapply(emat, 1, which) # slow

  nonzero <- function(x){
    ## function to get a two-column matrix containing the indices of the
    ### non-zero elements in a "dgCMatrix" class matrix

    stopifnot(inherits(x, "dgCMatrix"))
    if (all(x@p == 0))
      return(matrix(0, nrow=0, ncol=2,
                    dimnames=list(character(0), c("row","col"))))
    res <- cbind(x@i+1, rep(seq(dim(x)[2]), diff(x@p)))
    colnames(res) <- c("row", "col")
    res <- res[x@x != 0, , drop = FALSE]
    return(res)
  }

  split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])

  emat.sparse <- Matrix::t( SeuratObject::as.sparse(emat))
  which.cells2 <- as.data.frame(nonzero(x = emat.sparse))
  which.cells2$gene <- rownames(emat)[which.cells2$col]
  which.cells2 <- as_tibble(which.cells2)
  which.cells2 <- mylapply(split_tibble(which.cells2, col = "gene"), function(x){
    c(as.integer(unlist(x$row)))
  })
  n_cell <- unlist(lapply(which.cells2, length))
  which.cells2 <- which.cells2[n_cell != 0]
  which.zero <- rownames(emat)[!(rownames(emat)  %in% names(which.cells2))]
  which.cells2 <- c(which.cells2, sapply(which.zero,function(x) c()))
  which.cells2 <- which.cells2[rownames(emat)]

  if (!is.null(x_av)){
    which.cells <- which.cells2[names(which.cells2) %in% x_av]
  } else {
    which.cells <- which.cells2
  }

  # get co-expression probabilities
  miko_message("Quantifying co-expression probabilities...", verbose = verbose)
  pmat <- mysapply(which.cells2, function(x) sapply(which.cells, function(y) length(x) * length(y)))
  fmat <- pmat/(ncol(emat)^2)

  # get co-expression occurences
  miko_message("Counting co-expression occurences...", verbose = verbose)
  imat <- mysapply(which.cells2, function(x) sapply(which.cells, function(y) length(intersect(x, y))))

  # m2d helper function
  mat2df <- function(M) {
    if (!methods::is(M, "matrix")) stop("M must be a square matrix. (M is not a matrix).")
    if (nrow(M)!=ncol(M))   stop("M must be a square matrix. (M is not square).")
    if (is.null(colnames(M))) colnames(M) <- 1:ncol(M)
    if (is.null(rownames(M))) rownames(M) <- 1:ncol(M)
    if (!identical(rownames(M), colnames(M))) stop("rownames(M) != colnames(M)")
    xy <- rbind(admisc::permutations(colnames(M)),
                matrix(nrow = length(colnames(M)), ncol = 2, data = c(colnames(M), colnames(M))))
    data.frame(id1=xy[,1], id2=xy[,2], value=M[xy], stringsAsFactors = FALSE)
  }

  # specify null matrices
  if (is.null(dim(imat))){
    all.col <- names(imat)
    imat <- t(as.matrix(as.numeric(imat)));
    if (length(names(which.cells)) == dim(imat)[1]) {
      rownames(imat) <- names(which.cells);
      colnames(imat) <- names(which.cells2);
    }
    if (length(names(which.cells2)) == dim(imat)[1]) {
      rownames(imat) <- names(which.cells2);
      colnames(imat) <- names(which.cells);
    }

    fmat <- t(as.matrix(as.numeric(fmat)));
    if (length(names(which.cells)) == dim(fmat)[1]) {
      rownames(fmat) <- names(which.cells);
      colnames(fmat) <- names(which.cells2);
    }
    if (length(names(which.cells2)) == dim(fmat)[1]) {
      rownames(fmat) <- names(which.cells2);
      colnames(fmat) <- names(which.cells);
    }
  }

  miko_message("Computing co-dependency indices...", verbose = verbose)

  if (n.workers > 1){
    cl <- parallel::makeCluster(n.workers)
    doParallel::registerDoParallel(cl)

    bin.mat.vec <- list()
    bin.mat.vec <- unlist(foreach(i = 1:nrow(imat), .packages = c("dplyr", "fgsea", "plyr"))  %dopar% {
      return(unlist(lapply(1:ncol(imat), function(j){
        binom.test(x = imat[i,j], n = ncol(emat), p = fmat[i,j], alternative = "greater")[["p.value"]]
      })))
    })
    parallel::stopCluster(cl)
  } else {
    bin.mat.vec <- unlist(mylapply(1:nrow(imat), function(i){
      unlist(lapply(1:ncol(imat), function(j){
        binom.test(x = imat[i,j], n = ncol(emat), p = fmat[i,j], alternative = "greater")[["p.value"]]
      }))
    }))
  }


  # bin.mat.vec <- unlist(mylapply(1:nrow(imat), function(i){
  #   unlist(lapply(1:ncol(imat), function(j){
  #     binom.test(x = imat[i,j], n = ncol(emat), p = fmat[i,j], alternative = "greater")[["p.value"]]
  #   }))
  # }))
  bin.mat <- matrix(ncol = ncol(imat), nrow = nrow(imat), data = bin.mat.vec, byrow = T)

  bin.mat.log <- -log10(bin.mat+ 1e-300) # 1e-300 offset included to avoid Inf values
  rownames(bin.mat.log) <- rownames(imat); colnames(bin.mat.log) <- colnames(imat)

  df.cdi <- tryCatch(mat2df(bin.mat.log), error = function(e){
    bin.mat.log.df <- as.data.frame(bin.mat.log)
    all.col <- colnames(bin.mat.log.df )
    bin.mat.log.df$gene <- rownames(bin.mat.log.df)
    bin.mat.log.df.long <- tidyr::pivot_longer(data = bin.mat.log.df, cols = all.col)
    return(bin.mat.log.df.long)
  })
  colnames(df.cdi) <- c("feature.x", "feature.y", "cdi")
  miko_message("Normalizing CDI...", verbose = verbose)
  df.cdi <- df.cdi %>%
    dplyr::group_by(feature.x) %>%
    dplyr::mutate(ncdi = cdi/max(cdi, na.rm = T),
                  denominator = max(cdi, na.rm = T))

  df.cdi <- df.cdi[!(df.cdi$feature.x == df.cdi$feature.y), ]
  df.cdi$p <- 10^(-1*df.cdi$cdi)
  df.cdi <- df.cdi %>% dplyr::group_by(feature.x) %>% dplyr::mutate(fdr = p.adjust(p, method = "BH"))

  df.cdi <- df.cdi %>% dplyr::filter(feature.x %in% x_av, feature.y %in% y_av)

  df.cdi$ncdi[df.cdi$cdi == 0] <- 0
  df.cdi <- df.cdi %>% dplyr::arrange(-ncdi)

  return(df.cdi)


}

