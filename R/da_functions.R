#' Run differential abundance analysis.
#'
#' Run differential abundance analysis. Adapted Milo method in a wrapper.
#'
#' @param object Seurat object
#' @param condition.group Character specifying which metadata column to perform differential analysis on.
#' @param sample.group Character specifying metadata column with replicate groupings. If not specified (Default), pseudoreplicates are generated (see scMiko::pseudoReplicates()).
#' @param design.groups Character vector specifying model design. Default: c(condition.group, sample.group).
#' @param balance.samples Logical to balance sample sizes between condition.group. Call to scMiko::balanceSamples()
#' @param balance.samples Numeric specifying target sample size when balancing samples. 'balance.size' argument for scMiko::balanceSamples.
#' @param assay Character specifying which Seurat assay to use.
#' @param reference.group Character specifying which condition.group is reference.
#' @param fdr.correction Logical to perform FDR correction.
#' @name da_Run
#' @author Nicholas Mikolajewicz
#' @return list of results
da_Run <- function(object, condition.group, sample.group = NA, design.groups = c(condition.group, sample.group),
                   balance.samples = T, balance.size = NA, assay = DefaultAssay(object), reference.group = NULL, fdr.correction = F){


  require(SingleCellExperiment)
  require(igraph)
  require(miloR)
  require(Matrix)
  require(scater)
  require(Seurat)

  milo_setup_knn_datafunction <- function(x, k, d = 50, get.distance = FALSE, BNPARAM, BSPARAM,
                                          BPPARAM) {
    BiocNeighbors::findKNN(x[, c(1:d)], k = k, BNPARAM = BNPARAM, BPPARAM = BPPARAM,
                           get.distance = get.distance)
  }


  milo_neighborsToKNNGraph <- function (nn, directed = FALSE){
    start <- as.vector(row(nn))
    end <- as.vector(nn)
    interleaved <- as.vector(rbind(start, end))
    if (directed) {
      g <- igraph::make_graph(interleaved, directed = TRUE)
    }
    else {
      g <- igraph::make_graph(interleaved, directed = FALSE)
      g <- igraph::simplify(g, edge.attr.comb = "first")
    }
    g
  }

  miko_buildGraph <- function (x, k = 10, d = 50, get.distance = FALSE, reduced.dim = "PCA",
                               BNPARAM = BiocNeighbors::KmknnParam(), BSPARAM = bsparam(), BPPARAM = BiocParallel::SerialParam())
  {
    nn.out <- milo_setup_knn_datafunction(x = reducedDim(x, reduced.dim),
                                          d = d, k = k, BNPARAM = BNPARAM, BSPARAM = BSPARAM,
                                          BPPARAM = BPPARAM)
    # sink(file = "/dev/null")
    # gc()
    # sink(file = NULL)
    miko_message(paste0("Constructing kNN graph with k:", k))
    zee.graph <- milo_neighborsToKNNGraph(nn.out$index, directed = FALSE)
    graph(x) <- zee.graph
    if (isTRUE(get.distance)) {
      miko_message(paste0("Retrieving distances from ", k, " nearest neighbours"))
      old.dist <- matrix(0L, ncol = ncol(x), nrow = ncol(x))
      n.idx <- ncol(x)
      for (i in seq_along(1:n.idx)) {
        i.knn <- nn.out$index[i, ]
        i.dists <- nn.out$distance[i, ]
        old.dist[i, i.knn] <- i.dists
        old.dist[i.knn, i] <- i.dists
      }


      old.dist <- as(old.dist, "CsparseMatrix")
      # old.dist <- tryCatch({
      #   as(old.dist, "dgCMatrix")
      # },
      # error = function(e){
      #   return( as(old.dist, "CsparseMatrix"))
      # })

      nhoodDistances(x) <- old.dist
      # sink(file = "/dev/null")
      # gc()
      # sink(file = NULL)
    }
    x@.k <- k
    x
  }


  miko_graph_nhoods_from_adjacency <- function (nhs, nhood.adj, da.res, is.da, merge.discord = FALSE,
                                                lfc.threshold = NULL, overlap = 1, subset.nhoods = NULL) {
    if (is.null(colnames(nhs))) {
      warning("No names attributed to nhoods. Converting indices to names")
      colnames(nhs) <- as.character(c(1:ncol(nhs)))
    }
    if (!is.null(subset.nhoods)) {
      if (mode(subset.nhoods) %in% c("character", "logical",
                                     "numeric")) {
        if (mode(subset.nhoods) %in% c("character", "numeric")) {
          sub.log <- colnames(nhs) %in% subset.nhoods
        }
        else {
          sub.log <- subset.nhoods
        }
        nhood.adj <- nhood.adj[sub.log, sub.log]
        if (length(is.da) == ncol(nhs)) {
          nhs <- nhs[sub.log]
          is.da <- is.da[sub.log]
          da.res <- da.res[sub.log, ]
        }
        else {
          stop("Subsetting `is.da` vector length does not equal nhoods length")
        }
      }
      else {
        stop(paste0("Incorrect subsetting vector provided:",
                    class(subset.nhoods)))
      }
    }
    else {
      if (length(is.da) != ncol(nhood.adj)) {
        stop("Subsetting `is.da` vector length is not the same dimension as adjacency")
      }
    }
    if (isFALSE(merge.discord)) {
      discord.sign <- sign(da.res[, "logFC"] %*% t(da.res[,
                                                          "logFC"])) < 0
      nhood.adj[discord.sign] <- 0
    }
    if (overlap > 1) {
      nhood.adj[nhood.adj < overlap] <- 0
    }
    if (!is.null(lfc.threshold)) {
      nhood.adj[, which(da.res$logFC < lfc.threshold)] <- 0
      nhood.adj[which(da.res$logFC < lfc.threshold), ] <- 0
    }
    nhood.adj <- as.matrix((nhood.adj > 0) + 0)
    n.dim <- ncol(nhood.adj)
    if (!isSymmetric(nhood.adj)) {
      stop("Overlap matrix is not symmetric")
    }
    if (nrow(nhood.adj) != ncol(nhood.adj)) {
      stop("Non-square distance matrix - check nhood subsetting")
    }
    g <- igraph::graph_from_adjacency_matrix(nhood.adj, mode = "undirected",
                                             diag = FALSE)
    groups <- igraph::cluster_louvain(g)$membership
    names(groups) <- colnames(nhood.adj)
    keep.groups <- intersect(unique(groups[is.da]), unique(groups))
    return(groups[groups %in% keep.groups])
  }



  miko_perform_lognormal_dge <- function (exprs.data, test.model, gene.offset = gene.offset,
                                          model.contrasts = NULL, n.coef = NULL) {
    if (isTRUE(gene.offset)) {
      n.gene <- apply(exprs.data, 2, function(X) sum(X > 0))
      old.col <- colnames(test.model)
      if (all(test.model[, 1] == 1)) {
        test.model <- cbind(test.model[, 1], n.gene, test.model[,
                                                                c(2:ncol(test.model))])
        colnames(test.model) <- c(old.col[1], "NGenes",
                                  old.col[c(2:length(old.col))])
      }
      else {
        test.model <- cbind(n.gene, test.model)
        colnames(test.model) <- c("NGenes", old.col)
      }
    }
    i.fit <- lmFit(exprs.data, test.model)
    if (!is.null(model.contrasts)) {
      mod.constrast <- makeContrasts(contrasts = model.contrasts,
                                     levels = test.model)
      i.fit <- contrasts.fit(i.fit, contrasts = mod.constrast)
      i.fit <- eBayes(i.fit, trend = TRUE)
      i.res <- as.data.frame(topTreat(i.fit, number = Inf,
                                      sort.by = "p", p.value = 1))
    }
    else {
      i.fit <- eBayes(i.fit, trend = TRUE)
      if (is.null(n.coef)) {
        n.coef <- ncol(test.model)
      }
      i.res <- as.data.frame(topTreat(i.fit, coef = ncol(test.model),
                                      number = Inf, sort.by = "p", p.value = 1))
    }
    return(i.res)
  }


  miko_group_nhoods_by_overlap <- function (nhs, da.res, is.da, overlap = 1, lfc.threshold = NULL,
                                            merge.discord = FALSE, subset.nhoods = NULL, cells = NULL)
  {
    nhood.adj <- miko_build_nhood_adjacency(nhs)
    groups <- miko_graph_nhoods_from_adjacency(nhs = nhs, nhood.adj = nhood.adj,
                                               is.da = is.da, da.res = da.res, subset.nhoods = subset.nhoods,
                                               overlap = overlap, lfc.threshold = lfc.threshold, merge.discord = merge.discord)
    return(groups)
  }

  miko_build_nhood_adjacency <- function (nhoods, overlap = 1)
  {
    nh_intersect_mat <- Matrix::crossprod(nhoods)
    nh_intersect_mat[nh_intersect_mat < overlap] <- 0
    rownames(nh_intersect_mat) <- colnames(nhoods)
    colnames(nh_intersect_mat) <- colnames(nhoods)
    return(nh_intersect_mat)
  }

  miko_findNhoodMarkers <- function (x, da.res, da.fdr = 0.1, assay = "logcounts", aggregate.samples = FALSE,
                                     sample_col = NULL, overlap = 1, lfc.threshold = NULL, merge.discord = FALSE,
                                     subset.row = NULL, gene.offset = TRUE, return.groups = FALSE,
                                     subset.nhoods = NULL, na.function = "na.pass", compute.new = FALSE) {
    if (!is(x, "Milo")) {
      stop("Unrecognised input type - must be of class Milo")
    }   else if (any(!assay %in% assayNames(x))) {
      stop(paste0("Unrecognised assay slot: ", assay))
    }


    if (is.null(na.function)) {
      warning("NULL passed to na.function, using na.pass")
      na.func <- get("na.pass")
    }   else {
      tryCatch({
        na.func <- get(na.function)
      }, warning = function(warn) {
        warning(warn)
      }, error = function(err) {
        stop(paste0("NA function ", na.function, " not recognised"))
      }, finally = {
      })
    }

    if (isTRUE(aggregate.samples) & is.null(sample_col)) {
      stop("if aggregate.samples is TRUE, the column storing sample information must be specified by setting 'sample_col'")
    }

    n.da <- sum(na.func(da.res$SpatialFDR < da.fdr))


    if (!is.na(n.da) & n.da == 0) {
      stop("No DA neighbourhoods found")
    }
    if (any(is.na(da.res$SpatialFDR))) {
      warning("NA values found in SpatialFDR vector")
    }
    miko_message(paste0("Found ", n.da, " DA neighbourhoods at FDR ",
                   da.fdr * 100, "%"))

    if ((ncol(nhoodAdjacency(x)) == ncol(nhoods(x))) & isFALSE(compute.new)) {
      miko_message("nhoodAdjacency found - using for nhood grouping")
      nhs.da.gr <- miko_graph_nhoods_from_adjacency(nhoods(x),
                                                    nhood.adj = nhoodAdjacency(x), da.res = da.res,
                                                    is.da = da.res$SpatialFDR < da.fdr, merge.discord = merge.discord,
                                                    lfc.threshold = lfc.threshold, overlap = overlap,
                                                    subset.nhoods = subset.nhoods)
    } else {
      miko_message("Computing nhood adjacency")
      nhs.da.gr <- miko_group_nhoods_by_overlap(nhoods(x), da.res = da.res,
                                                is.da = da.res$SpatialFDR < da.fdr, merge.discord = merge.discord,
                                                lfc.threshold = lfc.threshold, overlap = overlap,
                                                cells = seq_len(ncol(x)), subset.nhoods = subset.nhoods)
    }
    nhood.gr <- unique(nhs.da.gr)
    miko_message(paste0("Nhoods aggregated into ", length(nhood.gr),
                   " groups"))
    fake.meta <- data.frame(CellID = colnames(x), Nhood.Group = rep(NA,
                                                                    ncol(x)))
    rownames(fake.meta) <- fake.meta$CellID
    for (i in seq_along(nhood.gr)) {
      nhood.x <- names(which(nhs.da.gr == nhood.gr[i]))
      nhs <- nhoods(x)
      if (!is.null(subset.nhoods)) {
        nhs <- nhs[, subset.nhoods]
      }
      nhood.gr.cells <- rowSums(nhs[, nhood.x, drop = FALSE]) >
        0
      fake.meta[nhood.gr.cells, "Nhood.Group"] <- ifelse(is.na(fake.meta[nhood.gr.cells,
                                                                         "Nhood.Group"]), nhood.gr[i], NA)
    }
    x <- x[, !is.na(fake.meta$Nhood.Group)]
    fake.meta <- fake.meta[!is.na(fake.meta$Nhood.Group), ]
    if (!is.null(subset.row)) {
      x <- x[subset.row, , drop = FALSE]
    }
    exprs <- assay(x, assay)
    marker.list <- list()
    i.contrast <- c("TestTest - TestRef")
    if (length(nhood.gr) == 1) {
      if (sum(fake.meta$Nhood.Group == nhood.gr[1]) == nrow(fake.meta)) {
        warning("All graph neighbourhoods are in the same group - cannot perform DGE testing. Returning NULL")
        return(NULL)
      }
    }
    if (isTRUE(return.groups)) {
      group.meta <- fake.meta
    }
    if (isTRUE(aggregate.samples)) {
      fake.meta[, "sample_id"] <- colData(x)[[sample_col]]
      fake.meta[, "sample_group"] <- paste(fake.meta[, "sample_id"],
                                           fake.meta[, "Nhood.Group"], sep = "_")
      sample_gr_mat <- matrix(0, nrow = nrow(fake.meta), ncol = length(unique(fake.meta$sample_group)))
      colnames(sample_gr_mat) <- unique(fake.meta$sample_group)
      rownames(sample_gr_mat) <- rownames(fake.meta)
      for (s in colnames(sample_gr_mat)) {
        sample_gr_mat[which(fake.meta$sample_group == s),
                      s] <- 1
      }
      exprs_smp <- matrix(0, nrow = nrow(exprs), ncol = ncol(sample_gr_mat))
      if (assay == "counts") {
        summFunc <- rowSums
      }
      else {
        summFunc <- rowMeans
      }
      for (i in 1:ncol(sample_gr_mat)) {
        if (sum(sample_gr_mat[, i]) > 1) {
          exprs_smp[, i] <- summFunc(exprs[, which(sample_gr_mat[,
                                                                 i] > 0)])
        }
        else {
          exprs_smp[, i] <- exprs[, which(sample_gr_mat[,
                                                        i] > 0)]
        }
      }
      rownames(exprs_smp) <- rownames(exprs)
      colnames(exprs_smp) <- colnames(sample_gr_mat)
      smp_meta <- unique(fake.meta[, c("sample_group", "Nhood.Group")])
      rownames(smp_meta) <- smp_meta[, "sample_group"]
      fake.meta <- smp_meta
      exprs <- exprs_smp
    }
    for (i in seq_along(nhood.gr)) {
      i.meta <- fake.meta
      i.meta$Test <- "Ref"
      i.meta$Test[fake.meta$Nhood.Group == nhood.gr[i]] <- "Test"
      if (ncol(exprs) > 1 & nrow(i.meta) > 1) {
        i.design <- as.formula(" ~ 0 + Test")
        i.model <- model.matrix(i.design, data = i.meta)
        rownames(i.model) <- rownames(i.meta)
      }
      # sink(file = "/dev/null")
      # gc()
      # sink(file = NULL)
      if (assay == "logcounts") {
        i.res <- miko_perform_lognormal_dge(exprs, i.model,
                                            model.contrasts = i.contrast, gene.offset = gene.offset)
      } else if (assay == "counts") {
        i.res <- miko_perform_lognormal_dge(exprs, i.model, model.contrasts = i.contrast,
                                            gene.offset = gene.offset)
        colnames(i.res)[ncol(i.res)] <- "adj.P.Val"
      } else {
        warning("Assay type is not counts or logcounts - assuming (log)-normal distribution. Use these results at your peril")
        i.res <- miko_perform_lognormal_dge(exprs, i.model,
                                            model.contrasts = i.contrast, gene.offset = gene.offset)
      }
      i.res$adj.P.Val[is.na(i.res$adj.P.Val)] <- 1
      i.res$logFC[is.infinite(i.res$logFC)] <- 0
      i.res <- i.res[, c("logFC", "adj.P.Val")]
      colnames(i.res) <- paste(colnames(i.res), nhood.gr[i],
                               sep = "_")
      marker.list[[paste0(nhood.gr[i])]] <- i.res
      # sink(file = "/dev/null")
      # gc()
      # sink(file = NULL)
    }
    marker.df <- do.call(cbind.data.frame, marker.list)
    colnames(marker.df) <- gsub(colnames(marker.df), pattern = "^[0-9]+\\.",
                                replacement = "")
    marker.df$GeneID <- rownames(i.res)
    if (isTRUE(return.groups)) {
      out.list <- list(groups = group.meta, dge = marker.df)
      return(out.list)
    }   else {
      return(marker.df)
    }
  }



  miko_calcNhoodDistance <- function(x, d, reduced.dim = NULL, use.assay = "logcounts"){
    non.zero.nhoods <- which(as.matrix(nhoods(x) != 0), arr.ind = TRUE)
    if (any(names(reducedDims(x)) %in% c("PCA"))) {
      nhood.dists <- sapply(1:ncol(nhoods(x)), function(X) .calc_distance(reducedDim(x,
                                                                                     "PCA")[non.zero.nhoods[non.zero.nhoods[, "col"] ==
                                                                                                              X, "row"], c(1:d), drop = FALSE]))
      names(nhood.dists) <- nhoodIndex(x)
    }
    else if (is.character(reduced.dim)) {
      nhood.dists <- sapply(1:ncol(nhoods(x)), function(X) .calc_distance(reducedDim(x,
                                                                                     reduced.dim)[non.zero.nhoods[non.zero.nhoods[, "col"] ==
                                                                                                                    X, "row"], c(1:d), drop = FALSE]))
      names(nhood.dists) <- nhoodIndex(x)
    }
    # sink(file = "/dev/null")
    # gc()
    # sink(file = NULL)
    nhoodDistances(x) <- nhood.dists
    return(x)
  }


  miko_testDiffExp <- function (x, da.res, design, meta.data, da.fdr = 0.1, model.contrasts = NULL,
                                overlap = 1, lfc.threshold = NULL, assay = "logcounts",
                                subset.row = NULL, gene.offset = TRUE, n.coef = NULL, merge.discord = FALSE,
                                na.function = "na.pass", compute.new = FALSE)
  {
    if (!is(x, "Milo")) {
      stop("Unrecognised input type - must be of class Milo")
    }
    else if (any(!assay %in% assayNames(x))) {
      stop(paste0("Unrecognised assay slot: ", assay))
    }
    if (is.null(na.function)) {
      warning("NULL passed to na.function, using na.pass")
      na.func <- get("na.pass")
    }
    else {
      tryCatch({
        na.func <- get(na.function)
      }, warning = function(warn) {
        warning(warn)
      }, error = function(err) {
        stop(paste0("NA function ", na.function, " not recognised"))
      }, finally = {
      })
    }
    n.da <- sum(na.func(da.res$SpatialFDR < da.fdr))
    if (!is.na(n.da) & n.da == 0) {
      stop("No DA neighbourhoods found")
    }
    if (any(is.na(da.res$SpatialFDR))) {
      warning("NA values found in SpatialFDR vector")
    }
    miko_message(paste0("Found ", n.da, " DA neighbourhoods at FDR ",
                   da.fdr * 100, "%"))
    if (!is.null(nhoodAdjacency(x)) & isFALSE(compute.new)) {
      miko_message("nhoodAdjacency found - using for nhood grouping")
      nhs.da.gr <- miko_group_nhoods_from_adjacency(nhoods(x),
                                                    nhood.adj = nhoodAdjacency(x), da.res = da.res,
                                                    is.da = na.func(da.res$SpatialFDR < da.fdr), overlap = overlap)
    }
    else {
      nhs.da.gr <- miko_group_nhoods_by_overlap(nhoods(x), da.res = da.res,
                                                is.da = na.func(da.res$SpatialFDR < da.fdr), overlap = overlap)
    }
    copy.meta <- meta.data
    if (!all(rownames(copy.meta) == colnames(x))) {
      warning("Column names of x are not the same as meta-data rownames")
    }
    copy.meta$Nhood.Group <- NA
    nhood.gr <- unique(nhs.da.gr)
    for (i in seq_along(nhood.gr)) {
      nhood.x <- nhs.da.gr %in% nhood.gr[i]
      copy.meta[rowSums(nhoods(x)[, names(nhs.da.gr)])[nhood.x] >
                  0, ]$Nhood.Group <- nhood.gr[i]
    }
    subset.dims <- !is.na(copy.meta$Nhood.Group)
    x <- x[, subset.dims]
    copy.meta <- copy.meta[subset.dims, ]
    if (is(design, "formula")) {
      model <- model.matrix(design, data = copy.meta)
      rownames(model) <- rownames(copy.meta)
    }
    else if (is.matrix(design)) {
      model <- design
      if (nrow(model) != nrow(copy.meta)) {
        miko_message("Subsetting input design matrix to DA neighbourhood cells")
        if (length(subset.dims) == nrow(model)) {
          model <- model[subset.dims, ]
        }
        else {
          stop(paste0("Cannot subset model matrix, subsetting vector is wrong length:",
                      length(subset.dims)))
        }
      }
      if (any(rownames(model) != rownames(copy.meta))) {
        warning("Design matrix and meta-data dimnames are not the same")
      }
    }
    if (ncol(x) != nrow(model)) {
      stop(paste0("Design matrix (", nrow(model), ") and milo objects (",
                  ncol(x), ") are not the same dimension"))
    }
    if (!is.null(subset.row)) {
      x <- x[subset.row, , drop = FALSE]
    }
    # sink(file = "/dev/null")
    # gc()
    # sink(file = NULL)
    miko_message(paste0("Nhoods aggregated into ", length(nhood.gr),
                   " groups"))
    dge.list <- list()
    for (i in seq_along(nhood.gr)) {
      i.meta <- copy.meta[copy.meta$Nhood.Group == nhood.gr[i],
                          , drop = FALSE]
      i.exprs <- assay(x[, copy.meta$Nhood.Group == nhood.gr[i],
                         drop = FALSE], assay)
      i.model <- model[copy.meta$Nhood.Group == nhood.gr[i],
                       , drop = FALSE]
      if (!is.null(ncol(i.exprs))) {
        if (ncol(i.exprs) > (ncol(i.model) + 1) & nrow(i.meta) >
            (ncol(i.model) + 1)) {
          if (assay == "logcounts") {
            i.res <- miko_perform_lognormal_dge(i.exprs, i.model,
                                                model.contrasts = model.contrasts, gene.offset = gene.offset,
                                                n.coef = n.coef)
          }
          else if (assay == "counts") {
            i.res <- miko_perform_counts_dge(i.exprs, i.model,
                                             model.contrasts = model.contrasts, gene.offset = gene.offset,
                                             n.coef = n.coef)
          }
          else {
            warning("Assay type is not counts or logcounts - assuming (log)-normal distribution. Use these results at your peril")
            i.res <- .perform_lognormal_dge(i.exprs, i.model,
                                            model.contrasts = model.contrasts, gene.offset = gene.offset,
                                            n.coef = n.coef)
          }
          i.res$adj.P.Val[is.na(i.res$adj.P.Val)] <- 1
          i.res$logFC[is.infinite(i.res$logFC)] <- 0
          i.res$Nhood.Group <- nhood.gr[i]
          dge.list[[paste0(nhood.gr[i])]] <- i.res
        }
        else if (ncol(i.exprs) <= ncol(i.model)) {
          warning("Not enough cells to perform DE testing in this neighbourhood")
        }
      }
      else {
        warning("Not enough cells to perform DE testing in this neighbourhood")
      }
    }
    return(dge.list)
  }




  # balance sample sizes
  if (balance.samples){
    object <- balanceSamples(object = object, group = condition.group, balance.size = balance.size)
  }

  # generate pseudoreplicates
  if (is.na(sample.group)){
    object <- pseudoReplicates(object = object, split.by = condition.group, n = 2)
    sample.group <- "pseudo_replicates"
    design.groups <- c(condition.group, sample.group)
  }

  # object <- FindNeighbors(object)
  object <- UpdateSeuratObject(object)
  sce.object <- tryCatch({
    sce.object <- as.SingleCellExperiment(object, assay = assay)
  }, error = function(e){
    object_diet <- DietSeurat(object = object, assays = assay, graphs = c("pca", "umap"))
    sce.object <- as.SingleCellExperiment(object_diet, assay = assay)
    return(sce.object)
  })

  milo.object <- miloR::Milo(sce.object)
  milo.object <- miko_buildGraph(milo.object, k = 10, d = 30)
  # milo.object <- miloR::makeNhoods(milo.object, prop = 0.5, k = 10, d = 30, refined = TRUE)
  milo.object <- miloR::makeNhoods(milo.object, prop = 0.5, k = 10, d = 30,
                                   refined = TRUE, refinement_scheme = "graph") # 040623

  milo.object <- countCells(milo.object, meta.data = data.frame(colData(milo.object)), sample=sample.group)

  milo_design <- data.frame(data.frame(colData(milo.object))[,design.groups])
  milo_design <- dplyr::distinct(milo_design)
  rownames(milo_design) <-milo_design[ ,sample.group]

  # milo.object <- miko_calcNhoodDistance(x = milo.object, d=30, reduced.dim = "PCA") # 040623


  form.cur <- c()

  for (i in 1:length(design.groups)){
    if (design.groups[i] %in% sample.group) next
    if (length(form.cur) == 0){
      form.cur <- paste0("~ ", design.groups[i])
    } else {
      form.cur <- paste0(form.cur , "+", design.groups[i])
    }
  }


  if (!is.null(reference.group)){
    all.var <- unique(as.character(milo_design[,condition.group]))
    all.var <- all.var[all.var != reference.group]
    all.var <- c(reference.group, all.var)
    milo_design[,condition.group] <- factor(milo_design[,condition.group], levels = all.var)
  }

  da.res <- testNhoods(milo.object, design = as.formula(form.cur) , design.df = milo_design, fdr.weighting = "graph-overlap") # 040623
  # da.res <- testNhoods(milo.object, design = as.formula(form.cur) , design.df = milo_design)

  if (!fdr.correction){
    da.res$SpatialFDR <- da.res$PValue
  }

  da.res$sP <- sign(da.res$logFC) * (-log10(da.res$PValue))
  milo.object <- buildNhoodGraph(milo.object, overlap = 1)

  signif_res <- da.res
  signif_res[signif_res$SpatialFDR > 0.1, "logFC"] <- 0
  colData(milo.object)["logFC"] <- NA

  suppressMessages({suppressWarnings({
    colData(milo.object)[unlist(nhoodIndex(milo.object)[signif_res$Nhood]), ] <- signif_res$sP
  })})

  nh_graph <- nhoodGraph(milo.object)
  umap.layout <- reducedDim(milo.object, "UMAP")[as.numeric(vertex_attr(nh_graph)$name), ]
  umap.layout <- as(umap.layout, "matrix")

  col_vals <- colData(milo.object)[as.numeric(vertex_attr(nh_graph)$name),  "logFC"]

  V(nh_graph)$colour_by <- col_vals

  # plot.g <- simplify(nh_graph)
  pl <- ggraph::ggraph(simplify(nh_graph), layout = umap.layout) +
    ggraph::geom_node_point(aes(fill = colour_by, size = size), shape = 21, color = "grey66") +
    scale_size(range = c(0.5, 5), name = "Nhood size") +
    ggraph::scale_edge_width(range = c(0.2, 3), name = "overlap size") +
    theme_classic(base_size = 14) +
    theme(axis.line = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank()) +
    scale_fill_gradient2(name = "signed(logP)",high = muted("red"),  mid = "white",  low = muted("blue"))

  if (!is.null(reference.group)){
    pl <- pl + labs(caption = paste0("Reference group: ", reference.group))
  }


  da.res <- annotateNhoods(milo.object, da.res, coldata_col = "seurat_clusters")

  return(list(
    milo.object = milo.object,
    seurat.object = object,
    da.umap = pl,
    results = da.res
  ))

}


#' Identify correlated genes/pathways associated with differential abundance.
#'
#' Identify correlated genes/pathways associated with differential abundance.
#'
#' @param da.res.list da_Run output.
#' @param reduction Character specifying which reduction to use for differential abundance association analysis.
#' @param species Character specifying species. Either "Mm" or "Hs". Default is "Mm".
#' @param red.gene.percentile Numeric, top nth percentile genes/pathways.
#' @param do.enrich Logical to perform enrichment of DA-associated genesets.
#' @param pathway.db Character pathway database to use for enrichment. Ignored if do.enrich = F. Default: "Bader".
#' @name da_DEG
#' @author Nicholas Mikolajewicz
#' @return updated da_Run output.
da_DEG <- function(da.res.list, reduction = "pca", species = "Mm", red.gene.percentile = 0.99, do.enrich = F,  pathway.db = "Bader"){

  milo.object <- da.res.list$milo.object
  seurat.object <- da.res.list$seurat.object

  if (is.na(reduction) | is.null(reduction)){
    # do nothing...
  } else {
    milo.exp <- calcNhoodExpression(x = nhoods(milo.object), subset.row = NULL, exprs = t(seurat.object[[reduction]]@cell.embeddings))
    milo.object@nhoodExpression <- milo.exp@nhoodExpression
  }

  hood.exp <- t(milo.object@nhoodExpression)

  hood.cor <- cor(hood.exp, da.res.list[["results"]]$logFC, method = "spearman")


  hood.z <- as.data.frame((hood.cor)/mad(hood.cor, na.rm = T))
  colnames(hood.z) <- "z"
  hood.z$r <- hood.cor[,1]
  hood.z$gene <- rownames(hood.z)
  hood.z <- hood.z[complete.cases(hood.z), ]


  hood.z.sig <- hood.z[abs(hood.z$z) > 2, ]
  hood.z.sig$gene2 <- gsub("_", "", hood.z.sig$gene)



  # getReductionGenes


  if (is.na(reduction) | is.null(reduction)){

  } else {

    try({
      feature.loading <- as.matrix(seurat.object@reductions[[reduction]]@feature.loadings)
      red.res <- getReductionGenes(feature.loading = feature.loading,  sig.threshold = red.gene.percentile)
      red.genes <- red.res$module.genes.pos

      if (do.enrich){
        hg.res.red <- runHG(gene.list =red.genes, gene.universe = rownames(seurat.object),
                            species = species,   pathway.db = pathway.db)
        hg.sum.red <- summarizeHG(hg.res.red)
        names(hg.sum.red$plots) <- gsub("_", "",names(hg.sum.red$plots)  )

        da.res.list[["hg.results"]] <- hg.res.red
        da.res.list[["hg.summary"]] <- hg.sum.red
      }


      da.res.list[["reduction.genes"]] <- red.genes
    }, silent = T)


  }


  da.res.list[["milo.object"]] <- milo.object
  da.res.list[["seurat.object"]] <- seurat.object

  da.res.list[["da.deg.results"]] <- list(all.results = hood.z, sig.results = hood.z.sig)

  # hood.z

  return(da.res.list)


}

