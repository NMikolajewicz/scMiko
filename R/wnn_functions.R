
#' Run WNN Multi-Modal Integration.
#'
#' Run WNN Multi-Modal Integration. Modified wrapper for seurat WNN workflow.
#'
#' @param object Seurat object or list of expression matrices. If seurat object, expression matrices are extracted. If list, assumes that expression matrix entries have column-wise genes and row-wise cells.
#' @param wnn.knn the number of multimodal neighbors to compute. 20 by default
#' @param umap.knn This determines the number of neighboring points used in local approximations of manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50. Default: 20
#' @param umap.min.dist This controls how tightly the embedding is allowed compress points together. Larger values ensure embedded points are moreevenly distributed, while smaller values allow the algorithm to optimise more accurately with regard to local structure. Sensible values are in the range 0.001 to 0.5. Default: 0.1
#' @param do.scale Logical to scale expression. Default is F.
#' @param do.center Logical to center expression. Default is F.
#' @param normalize.margin If specified, normalize across rows/cells (1) or columns/genes (2)
#' @param pca.thres Variance explained threshold for PC component inclusion. Default is 0.9.
#' @param cluster.resolution Cluster resolution for integrated network. Default is 1.
#' @param cluster.algorithm Algorithm for modularity optimization (1 = original Louvain algorithm; 2 = Louvain algorithm with multilevel refinement; 3 = SLM algorithm; 4 = Leiden algorithm). See Seurat:FindClusters() for details. Default: 3.
#' @param min.pct Minimum expression fraction for inclusion in network integration. Default is 0.25. Ignored if object is list.
#' @param split.var Grouping variable for expression fraction filter. Default is 'seurat_clusters'. Ignored if object is list.
#' @param neighborhood.membership Logical whether to return list of local neighborhoods. Default: T.
#' @name wnn_Run
#' @author Nicholas Mikolajewicz
#' @return list of integrated results
wnn_Run <- function (object, wnn.knn = 20, umap.knn = 20, umap.min.dist = 0.1, do.scale = F, do.center = F, normalize.margin = NA, pca.thres = 0.9, cluster.resolution = 1,  cluster.algorithm = 3, min.pct = 0.25, split.var = "seurat_clusters", neighborhood.membership = T){

  suppressMessages({
    suppressWarnings({


  require(Seurat)
  require(pbapply)
  require(future)
  require(future.apply)
      require(Matrix)

  wnn_Sweep <- function (x, MARGIN, STATS, FUN = "-", check.margin = TRUE,  ...){
    if (any(grepl(pattern = "X", x = names(x = formals(fun = sweep))))) {
      return(sweep(X = x, MARGIN = MARGIN, STATS = STATS,
                   FUN = FUN, check.margin = check.margin, ...))
    }
    else {
      return(sweep(x = x, MARGIN = MARGIN, STATS = STATS,
                   FUN = FUN, check.margin = check.margin, ...))
    }
  }

  wnn_L2Norm <-  function (mat, MARGIN = 1) {
    normalized <- wnn_Sweep(x = mat, MARGIN = MARGIN, STATS = apply(X = mat,
                                                                    MARGIN = MARGIN, FUN = function(x) {
                                                                      sqrt(x = sum(x^2))
                                                                    }), FUN = "/")
    normalized[!is.finite(x = normalized)] <- 0
    return(normalized)
  }

  wnn_Neighbor <- new("classGeneratorFunction", .Data = function (...)
    new("Neighbor", ...), className = "Neighbor", package = "Seurat")

  wnn_NNHelper <- function (data, query = data, k, method, cache.index = FALSE,
                            ...){
    args <- as.list(x = sys.frame(which = sys.nframe()))
    args <- c(args, list(...))
    results <- (switch(EXPR = method, rann = {
      args <- args[intersect(x = names(x = args), y = names(x = formals(fun = nn2)))]
      do.call(what = "nn2", args = args)
    }, annoy = {
      args <- args[intersect(x = names(x = args), y = names(x = formals(fun = wnn_AnnoyNN)))]
      do.call(what = "wnn_AnnoyNN", args = args)
    }, stop("Invalid method. Please choose one of 'rann', 'annoy'")))

    #   mat <- Matrix::rsparsematrix(nrow = 10, ncol = 10, density = 0.1)
    # rownames(x = mat) <- paste0("feature_", 1:10)
    # colnames(x = mat) <- paste0("cell_", 1:10)

    n.ob <- wnn_Neighbor(nn.idx = results$nn.idx, nn.dist = results$nn.dists,
                         alg.info = results$alg.info %||% list(), cell.names = rownames(x = query))
    if (isTRUE(x = cache.index) && !is.null(x = results$idx)) {
      slot(object = n.ob, name = "alg.idx") <- results$idx
    }
    return(n.ob)
  }

  wnn_AnnoyNN <- function (data, query = data, metric = "euclidean", n.trees = 50,
                           k, search.k = -1, include.distance = TRUE, index = NULL){
    idx <- index %||% wnn_AnnoyBuildIndex(data = data, metric = metric,
                                          n.trees = n.trees)
    nn <- wnn_AnnoySearch(index = idx, query = query, k = k, search.k = search.k,
                          include.distance = include.distance)
    nn$idx <- idx
    nn$alg.info <- list(metric = metric, ndim = ncol(x = data))
    return(nn)
  }

  wnn_AnnoyBuildIndex <- function (data, metric = "euclidean", n.trees = 50) {
    f <- ncol(x = data)
    a <- switch(EXPR = metric, euclidean = new(Class = RcppAnnoy::AnnoyEuclidean,
                                               f), cosine = new(Class = RcppAnnoy::AnnoyAngular, f),
                manhattan = new(Class = RcppAnnoy::AnnoyManhattan, f),
                hamming = new(Class = RcppAnnoy::AnnoyHamming, f), stop("Invalid metric"))
    for (ii in seq(nrow(x = data))) {
      a$addItem(ii - 1, data[ii, ])
    }
    a$build(n.trees)
    return(a)
  }


  wnn_AnnoySearch <- function (index, query, k, search.k = -1, include.distance = TRUE) {
    n <- nrow(x = query)
    idx <- matrix(nrow = n, ncol = k)
    dist <- matrix(nrow = n, ncol = k)
    convert <- methods::is(index, "Rcpp_AnnoyAngular")
    if (!inherits(x = plan(), what = "multicore")) {
      oplan <- plan(strategy = "sequential")
      on.exit(plan(oplan), add = TRUE)
    }
    res <- future_lapply(X = 1:n, FUN = function(x) {
      res <- index$getNNsByVectorList(query[x, ], k, search.k,
                                      include.distance)
      if (convert) {
        res$dist <- 0.5 * (res$dist * res$dist)
      }
      list(res$item + 1, res$distance)
    })
    for (i in 1:n) {
      idx[i, ] <- res[[i]][[1]]
      if (include.distance) {
        dist[i, ] <- res[[i]][[2]]
      }
    }
    return(list(nn.idx = idx, nn.dists = dist))
  }

  wnn_ComputeSNN <- function (nn_ranked, prune) {
    .Call("_Seurat_ComputeSNN", PACKAGE = "Seurat", nn_ranked,
          prune)
  }

  wnn_ComputeSNNwidth <- function (snn.graph, embeddings, k.nn, l2.norm = TRUE, nearest.dist = NULL) {
    if (l2.norm) {
      embeddings <- wnn_L2Norm(mat = embeddings)
    }
    nearest.dist <- nearest.dist %||% rep(x = 0, times = ncol(x = snn.graph))
    if (length(x = nearest.dist) != ncol(x = snn.graph)) {
      stop("Please provide a vector for nearest.dist that has as many elements as",
           " there are columns in the snn.graph (", ncol(x = snn.graph),
           ").")
    }
    snn.width <- SNN_SmallestNonzero_Dist(snn = snn.graph, mat = embeddings,
                                          n = k.nn, nearest_dist = nearest.dist)
    return(snn.width)
  }


  SNN_SmallestNonzero_Dist <- function (snn, mat, n, nearest_dist) {
    .Call("_Seurat_SNN_SmallestNonzero_Dist", PACKAGE = "Seurat",
          snn, mat, n, nearest_dist)
  }



  wnn_NNdist <- function (nn.idx, embeddings, metric = "euclidean", query.embeddings = NULL,
                          nearest.dist = NULL) {
    if (!is.list(x = nn.idx)) {
      nn.idx <- lapply(X = 1:nrow(x = nn.idx), FUN = function(x) nn.idx[x,
      ])
    }
    query.embeddings <- query.embeddings %||% embeddings
    nn.dist <- wnn_fast_dist(x = query.embeddings, y = embeddings,
                             n = nn.idx)
    if (!is.null(x = nearest.dist)) {
      nn.dist <- lapply(X = 1:nrow(x = query.embeddings),
                        FUN = function(x) {
                          r_dist = nn.dist[[x]] - nearest.dist[x]
                          r_dist[r_dist < 0] <- 0
                          return(r_dist)
                        })
    }
    return(nn.dist)
  }

  wnn_fast_dist <- function (x, y, n) {
    .Call("_Seurat_fast_dist", PACKAGE = "Seurat", x, y, n)
  }


  wnn_FindModalityWeights <- function (object, reduction.list, dims.list, k.nn = 20, snn.far.nn = TRUE,
                                       s.nn = k.nn, prune.SNN = 0, l2.norm = TRUE, sd.scale = 1,
                                       query = NULL, cross.contant.list = as.list(rep(1e-04, length(reduction.list))), sigma.idx = k.nn,
                                       smooth = FALSE, verbose = TRUE){
    my.lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    reduction.set <- unlist(x = reduction.list)
    names(x = reduction.list) <- names(x = dims.list) <- names(x = cross.contant.list) <- reduction.set
    embeddings.list <- lapply(X = reduction.list, FUN = function(r) Embeddings(object = object,
                                                                               reduction = r)[, dims.list[[r]]])
    if (l2.norm) {
      embeddings.list.norm <- lapply(X = embeddings.list,
                                     FUN = function(embeddings) wnn_L2Norm(mat = embeddings))
    } else {
      embeddings.list.norm <- embeddings.list
    }
    if (is.null(x = query)) {
      query.embeddings.list.norm <- embeddings.list.norm
      query <- object
    } else {
      if (snn.far.nn) {
        stop("query does not support to use snn to find distant wnn_Neighbors")
      }
      query.embeddings.list <- lapply(X = reduction.list,
                                      FUN = function(r) {
                                        Embeddings(object = query, reduction = r)[,
                                                                                  dims.list[[r]]]
                                      })
      if (l2.norm) {
        query.embeddings.list <- lapply(X = query.embeddings.list,
                                        FUN = function(embeddings) wnn_L2Norm(mat = embeddings))
      }
      query.embeddings.list.norm <- query.embeddings.list
    }
    if (verbose) {
      message("Finding ", k.nn, " nearest wnn_Neighbors for each modality.")
    }
    nn.list <- my.lapply(X = reduction.list, FUN = function(r) {
      nn.r <- wnn_NNHelper(data = embeddings.list.norm[[r]], query = query.embeddings.list.norm[[r]],
                           k = max(k.nn, sigma.idx, s.nn), method = "annoy",
                           metric = "euclidean")
      return(nn.r)
    })
    sigma.nn.list <- nn.list
    if (sigma.idx > k.nn || s.nn > k.nn) {
      nn.list <- lapply(X = nn.list, FUN = function(nn) {
        slot(object = nn, name = "nn.idx") <- Indices(object = nn)[,
                                                                   1:k.nn]
        slot(object = nn, name = "nn.dists") <- Distances(object = nn)[,
                                                                       1:k.nn]
        return(nn)
      })
    }
    nearest_dist <- lapply(X = reduction.list, FUN = function(r) Distances(object = nn.list[[r]])[,
                                                                                                  2])
    within_impute <- list()
    cross_impute <- list()
    for (r in reduction.set) {
      reduction.norm <- paste0(r, ".norm")
      object[[reduction.norm]] <- CreateDimReducObject(embeddings = embeddings.list.norm[[r]],
                                                       key = paste0("norm", Key(object = object[[r]])),
                                                       assay = DefaultAssay(object = object[[r]]))
      within_impute[[r]] <- PredictAssay(object = object,
                                         nn.idx = Indices(object = nn.list[[r]]), reduction = reduction.norm,
                                         dims = 1:ncol(x = embeddings.list.norm[[r]]), verbose = FALSE,
                                         return.assay = FALSE)

      cross.impute.groups <- setdiff(x = reduction.set, y = r)

      if (length(cross.impute.groups) == 1){
        cross_impute[[r]] <- PredictAssay(object = object,
                                          nn.idx = Indices(object =
                                                             nn.list[[setdiff(x = reduction.set, y = r)]]),
                                          reduction = reduction.norm, dims = 1:ncol(x = embeddings.list.norm[[r]]),
                                          verbose = FALSE, return.assay = FALSE)
      } else {
        cross_imput <- NULL
        for (kk in 1:length(cross.impute.groups)){
          cross_imput.cur <-  PredictAssay(object = object,
                                           nn.idx = Indices(object =
                                                              nn.list[[cross.impute.groups[[kk]]]]),
                                           reduction = reduction.norm, dims = 1:ncol(x = embeddings.list.norm[[r]]),
                                           verbose = FALSE, return.assay = FALSE)

          if (is.null(cross_imput)){
            cross_imput <- cross_imput.cur
          } else {
            cross_imput <- cross_imput + cross_imput.cur
          }

        }
        cross_impute[[r]] <- cross_imput/length(cross.impute.groups)
      }

    }

    within_impute_dist <- lapply(X = reduction.list, FUN = function(r) {
      r_dist <- sqrt(x = rowSums(x = (query.embeddings.list.norm[[r]] -
                                        t(x = within_impute[[r]]))^2))
      r_dist <- r_dist - nearest_dist[[r]]
      r_dist[r_dist < 0] <- 0
      return(r_dist)
    })
    cross_impute_dist <- lapply(X = reduction.list, FUN = function(r) {
      r_dist <- sqrt(x = rowSums(x = (query.embeddings.list.norm[[r]] -
                                        t(x = cross_impute[[r]]))^2))
      r_dist <- r_dist - nearest_dist[[r]]
      r_dist[r_dist < 0] <- 0
      return(r_dist)
    })
    if (snn.far.nn) {
      if (verbose) {
        message("Calculating kernel bandwidths")
      }
      snn.graph.list <- lapply(X = sigma.nn.list, FUN = function(nn) {
        snn.matrix <- wnn_ComputeSNN(nn_ranked = Indices(object = nn)[,
                                                                      1:s.nn], prune = prune.SNN)
        colnames(x = snn.matrix) <- rownames(x = snn.matrix) <- Cells(x = object)
        return(snn.matrix)
      })
      farthest_nn_dist <- my.lapply(X = 1:length(x = snn.graph.list),
                                    FUN = function(s) {
                                      distant_nn <- wnn_ComputeSNNwidth(snn.graph = snn.graph.list[[s]],
                                                                        k.nn = k.nn, l2.norm = FALSE, embeddings = embeddings.list.norm[[s]],
                                                                        nearest.dist = nearest_dist[[s]])
                                      return(distant_nn)
                                    })
      names(x = farthest_nn_dist) <- unlist(x = reduction.list)
      modality_sd.list <- lapply(X = farthest_nn_dist, FUN = function(sd) sd *
                                   sd.scale)
    } else {
      if (verbose) {
        message("Calculating sigma by ", sigma.idx, "th wnn_Neighbor")
      }
      modality_sd.list <- lapply(X = reduction.list, FUN = function(r) {
        rdist <- Distances(object = sigma.nn.list[[r]])[,
                                                        sigma.idx] - nearest_dist[[r]]
        rdist <- rdist * sd.scale
        return(rdist)
      })
    }
    within_impute_kernel <- lapply(X = reduction.list, FUN = function(r) {
      exp(-1 * (within_impute_dist[[r]]/modality_sd.list[[r]]))
    })
    cross_impute_kernel <- lapply(X = reduction.list, FUN = function(r) {
      exp(-1 * (cross_impute_dist[[r]]/modality_sd.list[[r]]))
    })
    params <- list(reduction.list = reduction.list, dims.list = dims.list,
                   l2.norm = l2.norm, k.nn = k.nn, sigma.idx = sigma.idx,
                   snn.far.nn = snn.far.nn, sigma.list = modality_sd.list,
                   nearest.dist = nearest_dist)
    modality_score <- lapply(X = reduction.list, FUN = function(r) {
      score <- within_impute_kernel[[r]]/(cross_impute_kernel[[r]] +
                                            cross.contant.list[[r]])
      score <- MinMax(data = score, min = 0, max = 200)
    })
    if (smooth) {
      modality_score <- lapply(X = reduction.list, FUN = function(r) {
        apply(X = Indices(object = nn.list[[r]]), MARGIN = 1,
              FUN = function(nn) mean(x = modality_score[[r]][nn[-1]]))
      })
    }


    # modality weights ###########################################################
    mod.mat <- matrix(nrow = length(modality_score[[1]]), ncol = length(reduction.set))
    for (i in 1:length(reduction.set)){
      mod.mat[ ,i] <- modality_score[[i]]
    }

    mod.mat.sum <- rowSums(exp(x = mod.mat))

    modality.weights.all <- list()
    for (i in 1:length(reduction.set)){
      modality.weights.all[[reduction.set[i]]] <- exp(x = modality_score[[i]])/mod.mat.sum
      names(modality.weights.all[[reduction.set[i]]]) <- names(modality_score[[i]])
    }

    score.mat <- cbind(Reduce(f = cbind, x = within_impute_dist),
                       Reduce(f = cbind, x = cross_impute_dist), Reduce(f = cbind,
                                                                        x = within_impute_kernel), Reduce(f = cbind, x = cross_impute_kernel),
                       Reduce(f = cbind, x = modality_score))


    cheaders <- (c("modalityX_nnZ",  "modalityX_nn", "modalityX_nnZ_kernel", "modalityX_nn_kernel", "modalityX_score"))


    cheaders.all <- c()
    for (i in 1:length(cheaders)){
      for (j in 1:length(reduction.set)){
        if (i == 1){
          cheaders.all <- c(cheaders.all, gsub( "Z", j, gsub("X", j,  cheaders[i])))
        } else if (i == 2){
          cheaders.all <- c(cheaders.all, gsub("X", j,  cheaders[i]))
        } else if (i == 3){
          cheaders.all <- c(cheaders.all, gsub( "Z", j, gsub("X", j,  cheaders[i])))
        } else  if (i == 4){
          cheaders.all <- c(cheaders.all, gsub("X", j,  cheaders[i]))
        } else  if (i == 5){
          cheaders.all <- c(cheaders.all, gsub("X", j,  cheaders[i]))
        }

      }
    }

    colnames(x = score.mat) <- cheaders.all

    score.mat <- as.data.frame(x = score.mat)

    modality.assay <- sapply(X = reduction.list, FUN = function(r) slot(object[[r]],
                                                                        name = "assay.used"))
    modality.weights <- new(Class = "ModalityWeights", first.modality.weight = modality.weights.all, # modality1.weight
                            modality.assay = modality.assay, params = params, score.matrix = score.mat) # , score.matrix = score.mat
    return(modality.weights)
  }

  wnn_MultiModalNN <- function (object, query = NULL, modality.weight = NULL, k.nn = NULL,
                                reduction.list = NULL, dims.list = NULL, knn.range = 200,
                                kernel.power = 1, nearest.dist = NULL, sigma.list = NULL,
                                l2.norm = NULL, verbose = TRUE) {
    my.lapply <- ifelse(test = verbose, yes = pblapply, no = lapply)
    k.nn <- k.nn %||% slot(object = modality.weight, name = "params")$k.nn
    reduction.list <- reduction.list %||% slot(object = modality.weight,
                                               name = "params")$reduction.list
    dims.list = dims.list %||% slot(object = modality.weight,
                                    name = "params")$dims.list
    nearest.dist = nearest.dist %||% slot(object = modality.weight,
                                          name = "params")$nearest.dist
    sigma.list = sigma.list %||% slot(object = modality.weight,
                                      name = "params")$sigma.list
    l2.norm = l2.norm %||% slot(object = modality.weight, name = "params")$l2.norm
    modality.weight.value <- modality.weight@first.modality.weight
    names(x = modality.weight.value) <- unlist(x = reduction.list)
    if (inherits(x = object, what = "Seurat")) {
      reduction_embedding <- lapply(X = 1:length(x = reduction.list),
                                    FUN = function(x) {
                                      Embeddings(object = object, reduction = reduction.list[[x]])[,
                                                                                                   dims.list[[x]]]
                                    })
    } else {
      reduction_embedding <- object
    }
    if (is.null(x = query)) {
      query.reduction_embedding <- reduction_embedding
      query <- object
    } else {
      if (inherits(x = object, what = "Seurat")) {
        query.reduction_embedding <- lapply(X = 1:length(x = reduction.list),
                                            FUN = function(x) {
                                              Embeddings(object = query, reduction = reduction.list[[x]])[,
                                                                                                          dims.list[[x]]]
                                            })
      }
      else {
        query.reduction_embedding <- query
      }
    }
    if (l2.norm) {
      query.reduction_embedding <- lapply(X = query.reduction_embedding,
                                          FUN = function(x) wnn_L2Norm(mat = x))
      reduction_embedding <- lapply(X = reduction_embedding,
                                    FUN = function(x) wnn_L2Norm(mat = x))
    }
    query.cell.num <- nrow(x = query.reduction_embedding[[1]])
    reduction.num <- length(x = query.reduction_embedding)
    if (verbose) {
      message("Finding multimodal wnn_Neighbors")
    }
    redunction_nn <- my.lapply(X = 1:reduction.num, FUN = function(x) {
      nn_x <- wnn_NNHelper(data = reduction_embedding[[x]], query = query.reduction_embedding[[x]],
                           k = knn.range, method = "annoy", metric = "euclidean")
      return(nn_x)
    })
    redunction_nn <- lapply(X = redunction_nn, FUN = function(x) Indices(object = x)[,
                                                                                     -1])
    nn_idx <- lapply(X = 1:query.cell.num, FUN = function(x) {
      Reduce(f = union, x = lapply(X = redunction_nn, FUN = function(y) y[x,
      ]))
    })
    nn_dist <- my.lapply(X = 1:reduction.num, FUN = function(r) {
      wnn_NNdist <- wnn_NNdist(nn.idx = nn_idx, embeddings = reduction_embedding[[r]],
                               query.embeddings = query.reduction_embedding[[r]],
                               nearest.dist = nearest.dist[[r]])
      return(wnn_NNdist)
    })
    if (length(x = sigma.list[[1]]) == 1) {
      sigma.list <- lapply(X = sigma.list, FUN = function(x) rep(x = x,
                                                                 ncol(x = object)))
    }
    nn_weighted_dist <- lapply(X = 1:reduction.num, FUN = function(r) {
      lapply(X = 1:query.cell.num, FUN = function(x) {
        exp(-1 * (nn_dist[[r]][[x]]/sigma.list[[r]][x])^kernel.power) *
          modality.weight.value[[r]][x]
      })
    })
    nn_weighted_dist <- sapply(X = 1:query.cell.num, FUN = function(x) {
      Reduce(f = "+", x = lapply(X = 1:reduction.num, FUN = function(r) nn_weighted_dist[[r]][[x]]))
    })
    select_order <- lapply(X = nn_weighted_dist, FUN = function(dist) {
      order(dist, decreasing = TRUE)
    })
    select_nn <- t(x = sapply(X = 1:query.cell.num, FUN = function(x) nn_idx[[x]][select_order[[x]]][1:k.nn]))
    select_dist <- t(x = sapply(X = 1:query.cell.num, FUN = function(x) nn_weighted_dist[[x]][select_order[[x]]][1:k.nn]))

    select_dist.inter <- (1 - select_dist)/2
    select_dist.inter[select_dist.inter < 0] <- 0
    select_dist <- sqrt(x = select_dist.inter)
    weighted.nn <- wnn_Neighbor(nn.idx = select_nn, nn.dist = select_dist,
                                alg.info = list(), cell.names = Cells(x = query))
    return(weighted.nn)
  }


  wnn_FindMultiModalNeighbors <- function (object, reduction.list, dims.list, k.nn = 20, l2.norm = TRUE,
                                           sd.scale = 1, cross.contant.list = as.list(rep(1e-04, length(reduction.list))), smooth = FALSE,
                                           knn.range = 200, knn.graph.name = "wknn", snn.graph.name = "wsnn",
                                           weighted.nn.name = "weighted.nn", modality.weight = NULL,
                                           prune.SNN = 1/15, weighted.graph = FALSE, return.intermediate = FALSE,
                                           modality.weight.name = NULL, verbose = TRUE){
    if (is.null(x = modality.weight)) {
      if (verbose) {
        message("Calculating cell-specific modality weights")
      }
      modality.weight <- wnn_FindModalityWeights(object = object,
                                                 reduction.list = reduction.list, dims.list = dims.list,
                                                 k.nn = k.nn, sd.scale = sd.scale, l2.norm = l2.norm,
                                                 cross.contant.list = cross.contant.list, smooth = smooth,
                                                 verbose = verbose)
    }
    modality.weight.name <- modality.weight.name %||% paste0(DefaultAssay(object = object[[reduction.list[[1]]]]),
                                                             ".weight")
    k.nn <- k.nn %||% slot(object = modality.weight, name = "params")$k.nn
    first.assay <- slot(object = modality.weight, name = "modality.assay")[1]
    weighted.nn <- wnn_MultiModalNN(object = object, k.nn = k.nn,
                                    modality.weight = modality.weight, knn.range = knn.range,
                                    verbose = verbose)
    select_nn <- Indices(object = weighted.nn)
    select_nn_dist <- Distances(object = weighted.nn)


    if (weighted.graph) {
      if (verbose) {
        message("Constructing multimodal KNN graph")
      }
      select_nn_dist <- t(x = apply(X = select_nn_dist, MARGIN = 1,
                                    FUN = function(x) log2(k.nn) * x/sum(x)))
      nn.matrix <- sparseMatrix(i = 1:ncol(x = object), j = 1:ncol(x = object),
                                x = 1)
      for (i in 1:ncol(x = object)) {
        nn.matrix[i, select_nn[i, ]] <- select_nn_dist[i,
        ]
        slot(object = weighted.nn, name = "nn.dist") <- select_nn_dist
      }
    } else {
      if (verbose) {
        message("Constructing multimodal KNN graph")
      }
      j <- as.numeric(x = t(x = select_nn))
      i <- ((1:length(x = j)) - 1)%/%k.nn + 1
      nn.matrix <- sparseMatrix(i = i, j = j, x = 1, dims = c(ncol(x = object),
                                                              ncol(x = object)))
      diag(x = nn.matrix) <- 1
    }
    rownames(x = nn.matrix) <- colnames(x = nn.matrix) <- colnames(x = object)
    nn.matrix <- nn.matrix + Matrix::t(x = nn.matrix) - Matrix::t(x = nn.matrix) *
      nn.matrix
    nn.matrix <- as.Graph(x = nn.matrix)
    slot(object = nn.matrix, name = "assay.used") <- first.assay
    object[[knn.graph.name]] <- nn.matrix
    if (verbose) {
      message("Constructing multimodal SNN graph")
    }
    snn.matrix <- wnn_ComputeSNN(nn_ranked = select_nn, prune = prune.SNN)
    rownames(x = snn.matrix) <- colnames(x = snn.matrix) <- Cells(x = object)
    snn.matrix <- as.Graph(x = snn.matrix)
    slot(object = snn.matrix, name = "assay.used") <- first.assay
    object[[snn.graph.name]] <- snn.matrix
    object[[weighted.nn.name]] <- weighted.nn


    for (i in 1:length(reduction.list)){
      object[[paste0(unlist(reduction.list)[i], "_", "weights")]] <- modality.weight@first.modality.weight[[i]]
    }

    return(object)
  }

    })
  }
  )

  if (class(object) == "Seurat"){

    message("Preparing expression matrices...")
    which.rep <- getExpressedGenes(object = object, min.pct = min.pct, group = split.var, group.boolean = "OR")

    # split seurat objects #########################################################
    object.list <- SplitObject(object[rownames(object) %in% which.rep, ], split.by = split.var)
    sample.n <- unlist(lapply(object.list, function(x) ncol(x)))

    # get sample specific e.mat
    exp.list <- list()

    for (i in 1:length(object.list)){
      exp.list[[names(object.list)[i]]] <- t(getExpressionMatrix(object.list[[i]], which.data = "data"))
    }

    rm(object)
    rm(object.list)


  } else if (class(object) == "list"){
    exp.list <- object
    rm(object)
  }

  invisible({gc()})

  # assemble seurat object #######################################################
  message("Preparing integration object...")
  for (i in 1:length(exp.list)){
    set.name <- names(exp.list)[i]
    if (i == 1){
      so.gene <- CreateSeuratObject(counts = (exp.list[[set.name]]), project = "multi.mod.integration", assay = set.name)
    } else {
      so.gene@assays[[set.name]] <- CreateAssayObject(counts = (exp.list[[set.name]]), min.cells = 0, min.features = 0)
    }
  }



  # normalize, scale and dimensionally reduce data ###############################
  message("Running PCA...")
  for (i in 1:length(exp.list)){
    set.name <- names(exp.list)[i]
    DefaultAssay(so.gene) <- set.name
    pca.name <- paste0("pc_", set.name)
    VariableFeatures(so.gene) <- rownames(so.gene[[set.name]])

    if (!is.na(normalize.margin)){
      so.gene <- Seurat::NormalizeData(so.gene, margin = normalize.margin, normalization.method = "CLR", verbose = F)
    }

    so.gene <- ScaleData(so.gene, do.scale =do.scale, do.center = do.center)

    nDim <- 50
    if (nrow(so.gene@assays[[set.name]]) < nDim) nDim <- nrow(so.gene@assays[[set.name]])
    so.gene <- RunPCA(so.gene, reduction.name = pca.name, npcs = nDim, verbose = F)

  }

  # get optimal PCA number #######################################################
  message("Identifying optimal number of PCs...")
  # pca.thres <- 0.9
  nDim.optimal <- c()
  for (i in 1:length(exp.list)){
    set.name <- names(exp.list)[i]
    DefaultAssay(so.gene) <- set.name
    pca.name <- paste0("pc_", set.name)
    pca.prop <- propVarPCA(so.gene, reduction.name = pca.name)
    nDim.optimal <- c(nDim.optimal, max(pca.prop$pc.id[pca.prop$pc.cum_sum<pca.thres])+1)
  }
  names(nDim.optimal) <- names(exp.list)

  dim.lists <- list()
  red.lists <- list()
  for (i in 1:length(nDim.optimal)){
    dim.lists[[names(nDim.optimal)[i]]] <- 1:nDim.optimal[i]
    set.name <- names(nDim.optimal)[i]
    pca.name <- paste0("pc_", set.name)
    red.lists[[i]] <- pca.name
  }


  # weighted nearest neighbor analysis ###########################################

  n.run <- length(red.lists)
  x.constant <- 1e-4 # 1e-4 default
  k.range <- 200 # 200 default

  message("Integrating object...")
  # Run 1 ######
  so.gene <- wnn_FindMultiModalNeighbors(
    object = so.gene, reduction.list = red.lists[1:n.run], k.nn = wnn.knn,  #  round(0.005 * ncol(so.gene))
    dims.list = dim.lists[1:n.run], modality.weight.name = "wnn.weight", knn.range = k.range,
    smooth = F,
    cross.contant.list =  as.list(rep(x.constant, length(red.lists[1:n.run]))),
    knn.graph.name = "wknn", snn.graph.name = "wsnn", weighted.nn.name = "weighted.nn"
  )


  do.filter = F

  if (do.filter){

  # Flag ill-connected genes and filter out
  graph.holder <- so.gene@graphs
  neighbor.holder <- so.gene@neighbors

  nn.dist.dat <- as.data.frame(so.gene@neighbors[["weighted.nn"]]@nn.dist)
  rownames(nn.dist.dat) <- so.gene@neighbors[["weighted.nn"]]@cell.names
  keep.which <- complete.cases(nn.dist.dat)
  so.gene2 <- so.gene[ ,keep.which]


  # try to use original neighbor graph. If error is encountered, compute new graph.

  which.cells.remain <- colnames(so.gene2)
  for (i in 1:length(graph.holder)){
    graph.name <- names(graph.holder)[i]
    current.graph <- graph.holder[[graph.name]]
    graph.class <- class(current.graph)
    current.graph <- current.graph[colnames(current.graph) %in% which.cells.remain,
                                   colnames(current.graph) %in% which.cells.remain]
    if (canCoerce(current.graph, graph.class)) {
      current.graph <- as(current.graph, graph.class)
      current.graph@assay.used <- DefaultAssay(so.gene2)
    }
    graph.holder[[graph.name]] <- current.graph
  }

  neighbor.holder2 <- list()

  keep.which2 <-  neighbor.holder[["weighted.nn"]]@cell.names %in% which.cells.remain
  neighbor.holder2[["weighted.nn"]] <-   wnn_Neighbor(nn.idx =  neighbor.holder[["weighted.nn"]]@nn.idx[keep.which2, ],
                                                      nn.dist =  neighbor.holder[["weighted.nn"]]@nn.dist[keep.which2, ],
                                                      alg.info = list(), cell.names = neighbor.holder[["weighted.nn"]]@cell.names[keep.which2])

  so.gene2@graphs <- graph.holder
  so.gene2@neighbors <- neighbor.holder2
  so.gene2 <- UpdateSeuratObject(so.gene2)


  } else {
    so.gene2 <- so.gene
    rm(so.gene)
  }

  message("Embedding UMAP...")

  so.gene2 <- tryCatch({
    so.gene2 <- RunUMAP(so.gene2, nn.name = "weighted.nn",
                        min.dist = umap.min.dist, n.neighbors = umap.knn,
                        reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
  }, error = function(e){
    so.gene2 <- RunUMAP(so.gene2, graph = "wknn", min.dist = umap.min.dist,
                        reduction.name = "wnn.umap", reduction.key = "wnnUMAP_" )
    return(so.gene2)
  })


  message("Finding clusters...")
  so.gene2 <- FindClusters(so.gene2, graph.name = "wsnn", algorithm = cluster.algorithm, resolution = cluster.resolution, verbose = T)

  wnnUMAP.list <- getUMAP(so.gene2, umap.key = "wnn.umap", node.type = "point", size = 0.01)

  df.wnn.umap <- wnnUMAP.list$df.umap
  plt.wnn.umap <- wnnUMAP.list$plt.umap
  plt.wnn.umap <- plt.wnn.umap +
    xlab("wnnUMAP 1") + ylab("wnnUMAP 2") +
    labs(title = "Gene Regulatory Network", subtitle = "Weighted nearest-neighbor analysis")

  # get nearest neighbor neighborhood ############################################

  gname <- so.gene2@neighbors[["weighted.nn"]]@cell.names
  wmat <- so.gene2@neighbors[["weighted.nn"]]@nn.dist
  nmat <- so.gene2@neighbors[["weighted.nn"]]@nn.idx

  neighborhood.list2 <- list()
  if (neighborhood.membership){
    message("Getting neighborhoods...")
    for (i in 1:length(gname)){
      neighborhood.list2[[gname[i]]] <- gname[nmat[i,]]
    }
  }

  return(
    list(
      so.gene = so.gene2,
      exp.list = exp.list,
      wnnUMAP.list = wnnUMAP.list,
      df.wnn.umap = df.wnn.umap,
      plt.wnn.umap = plt.wnn.umap,
      neighborhood.list = neighborhood.list2,
      dim.lists = dim.lists
    )
  )


}



#' Compute network component UMAPs and visualize component weights.
#'
#' Compute network component UMAPs and visualize component weights.
#'
#' @param object wnn integrated seurat object (wnn_Run output)
#' @param exp.list List of expression matrices (wnn_Run list input)
#' @param dim.lists List of PC numbers (wnn_Run output)
#' @param plt.wnn.umap ggplot handle for integrated network umap (wnn_Run output)
#' @param umap.min.dist min.dist for umap embedding. See RunUMAP(). Default: 0.1
#' @param umap.n.neighbors n.neighbor for umap embedding. See RunUMAP(). Default: 20
#' @name wnn_Components
#' @author Nicholas Mikolajewicz
#' @return list of ggplot handle and object
wnn_Components <- function(object, exp.list, dim.lists, plt.wnn.umap, umap.min.dist = 0.1, umap.n.neighbors = 20){


  # get component umaps ##########################################################
  # umap.min.dist <- 0.100
  # umap.n.neighbors <- 20

  which.node.type <- "point"
  # umap.lists <- list()
  plt.umap.list <- list()

  message("Embedding UMAPs...")
  for (i in 1:length(exp.list)){
    set.name <- names(exp.list)[i]

    pca.reduction.name <- paste0("pc_", set.name)
    # red.lists[[i]] <- pca.name

    # pca.reduction.name <- red.lists[[i]] #; analysis.parameters[[names(nDim.optimal)[i]]][["reduction.name"]]
    umap.reduction.name1 <- paste0(tolower(set.name), ".umap")
    umap.reduction.name2 <- paste0(tolower(set.name), "UMAP_")

    object <- RunUMAP(object, reduction = pca.reduction.name, dims = dim.lists[[set.name]], assay = set.name,
                        reduction.name = umap.reduction.name1, reduction.key = umap.reduction.name2,
                        min.dist = umap.min.dist, n.neighbors = umap.n.neighbors)

    UMAP.list.current <- getUMAP(object, umap.key = umap.reduction.name1, node.type = which.node.type, size = 0.01)

    plt.umap.current <- UMAP.list.current$plt.umap
    plt.umap.current <- plt.umap.current +
      xlab(paste0(tolower(set.name), "UMAP 1")) + ylab(paste0(tolower(set.name), "UMAP 2")) +
      labs(title = "Regulator network", subtitle = set.name) + theme_miko()
    # + scale_color_manual(values = div7.col)
    plt.umap.list[[set.name]] <- plt.umap.current
  }

  which.weights <- names(object@meta.data)[grepl("weights", names(object@meta.data))]

  message("Constructing plot...")
  plt.weight.vln <- list()
  for (i in which.weights){
    plt.weight.vln[[i]] <-   VlnPlot(object, features = i, group.by = 'seurat_clusters', sort = F, pt.size = 0.1) +
      NoLegend() + ylim(0, 1)  + theme_miko() +
      ylab("Weights") + xlab("Clusters") +
      labs(title = "WNN Weights")
  }

  plt.top <- cowplot::plot_grid(plotlist = plt.umap.list, ncol = length(plt.umap.list), align = "hv")
  plt.middle <- cowplot::plot_grid(plotlist = plt.weight.vln, ncol = length(plt.weight.vln), align = "hv")
  plt.merge <- cowplot::plot_grid(plt.top, plt.middle, ncol = 1)
  plt.all <- cowplot::plot_grid(plt.merge, plt.wnn.umap, ncol = 2, rel_widths = c(3,2))


  return(
    list(
      object = object,
      plot = plt.all
    )
  )

}

