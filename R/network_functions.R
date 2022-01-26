#' Identify features for gene program discovery.
#'
#' Identify features for gene program discovery.
#'
#' @param object Seurat object
#' @param method Feature selection method.
#' \itemize{
#' \item "expr" - Top expressed genes within `split_var` groups, together with variable genes.
#' \item "hvg" - Highly-variable genes (number specified by `n_features`)
#' \item "deviance" - Deviant genes (number specified by `n_features`)
#' }
#' @param n_features Number of features to use when calculating variable or deviant features. Default is 2000.
#' @param min_pct Minimum expression within spilt_var. Ignored if `method` is not "expr". Default is 0.5.
#' @param split_var Grouping variable to enforce `min_pct` criteria within. Default is 'seurat_clusters'.
#' @param batch_feature Variables to regress out. Ignored if `method` is not "deviance". Default is NULL.
#' @param verbose Show progress. Default is T.
#' @name findNetworkFeatures
#' @seealso \code{\link{runSSN}}
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @author Nicholas Mikolajewicz
#' @return vector of features
#' @examples
#'
findNetworkFeatures <- function(object, method = c("expr", "hvg", "deviance"), n_features = 2000, min_pct = 0.5, split_var = "seurat_clusters", batch_feature = NULL, verbose = T){

  miko_message("Selecting features for network analysis...", verbose = verbose)
  # get minimally expressed genes ################################################
  # default 0.2
  # split.var <- "seurat_clusters"


  if (!("Seurat") %in% class(object)) stop("'object' is not a Seurat object")
  if (!(split_var) %in% colnames(object@meta.data)) stop(paste0(split_var, " is not a meta feature in object"))
  if (is.na(n_features) | is.null(n_features)) n_features <- 2000

  if (method == "hvg"){
    miko_message("Getting highly variable genes...", verbose = verbose)
    try({
      object <- FindVariableFeatures(object = object, nfeatures = n_features)
    }, silent = T)

    which.rep.genes2 <- unique(VariableFeatures(object))
  } else if (method == "expr"){
    miko_message("Getting highly expressed and variable genes...", verbose = verbose)
    which.rep.genes <- getExpressedGenes(object = object, min.pct = min_pct, group = split_var, group.boolean = "OR")

    try({
      object <- FindVariableFeatures(object, nfeatures = n_features)
    }, silent = T)

    which.rep.genes2 <- unique(c(VariableFeatures(object), which.rep.genes))
  } else if (method == "deviance"){
    miko_message("Getting highly deviant genes...", verbose = verbose)
    set.seed(1023)
    so.sub <- downsampleSeurat(object = object, subsample.n = 20000, verbose = F)

    m <- GetAssayData(so.sub, slot = "data", assay = DefaultAssay(so.sub))
    m <- m[rownames(m) %in% rownames(so.sub), ]
    m <- as.matrix(m)

    if ((!is.null(batch_feature)) && batch_feature %in% colnames(object@meta.data)){
      batch.factor <- as.factor(so.sub@meta.data[ ,batch_feature])
      devs <- scry::devianceFeatureSelection(object = m, batch = batch.factor, fam = "binomial")
    } else {
      devs <- scry::devianceFeatureSelection(object = m, fam = "binomial")
    }

    df.dev <- data.frame(
      gene = rownames(m),
      d = devs) %>% dplyr::arrange(-d)

    df.dev$d.norm <- df.dev$d/sum(df.dev$d)
    df.dev$d.cumsum <- cumsum(df.dev$d.norm )

    dev.gene <- df.dev$gene[df.dev$d.cumsum < 0.8]
    if (length(dev.gene) > n_features){
      dev.gene <- df.dev$gene[1:n_features]
    }

    which.rep.genes2 <- dev.gene

    rm(m); rm(so.sub)
    invisible({gc()})

  } else {
    which.rep.genes <- getExpressedGenes(object = object, min.pct = min_pct, group = split_var, group.boolean = "OR")
    which.rep.genes2 <- unique(c(VariableFeatures(object), which.rep.genes))
  }

  which.rep.common <- which.rep.genes2

  return(which.rep.common)

}




#' Apply scale-free topology transform to shared-nearest neighbor graph
#'
#' Apply scale-free topology transform to shared-nearest neighbor graph
#'
#' @param object Seurat object
#' @param graph.name name of SNN graph in `object`. Default is "RNA_snn".
#' @param sf soft thresholds (sf) that are evaluated in search for optimal soft threshold (used to transform SNN to scale-free topology). Default is sf = c(seq(1, 5, by = 0.5), seq(6, 10)).
#' @param sf_threshold R2 threshold [0,1] for identifying optimal soft threshold. Default is 0.9 (recommended 0.8-0.9).
#' @param umap_knn This determines the number of neighboring points used in local approximations of UMAP manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50. default is 10.
#' @param n_dim Number of PC dimensions to use in generating UMAP. Default is PCA.
#' @param binary_threshold Numeric threshold used to binarize SNN graph. Default is 0.9.
#' @param n_workers Number of workers used for parallel implementation. Default is 1.
#' @param verbose Show progress. Default is T.
#' @name scaleFreeNet
#' @seealso \code{\link{runSSN}}
#' @author Nicholas Mikolajewicz
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @return list containing
#' \itemize{
#' \item "object" - scale-free SNN and binarized graph are stored in `object` graph slot, and umap embedding is stored in `object` reduction slot.
#' \item "st.res" - soft threshold optimization results.
#' }
#' @examples
#'
scaleFreeNet <- function(object, graph_name = "RNA_snn", sf = c(seq(1, 5, by = 0.5), seq(6, 10)), sf_threshold = 0.9, umap_knn = 10, n_dim = 30, binary_threshold = 0.9, n_workers = 1, verbose = T){

  suppressMessages({
    suppressWarnings({

      if (!(graph_name %in% names(object@graphs))) stop(paste0(graph_name, " is not a graph in 'object'"))

      # power.all <- c(seq(0.5, 5, by = 0.5), seq(6, 10))
      # n.work <-   parallel::detectCores()
      # if (length(power.all) < n.work) n.work <- length(power.all)
      miko_message("Identifying optimal power...", verbose = verbose)
      st.res <- getSoftThreshold2(
        s.mat = as.matrix(object@graphs[[graph_name]]),
        power = sf,
        network.type = "unsigned",
        nBreaks = 20,
        removeFirst = T,
        rescale.adjacency = F,
        n.cores = n_workers,
        r2.target = sf_threshold
      )


      new.graph <- paste0(graph_name, "_power")
      object@graphs[[new.graph]] <- object@graphs[[graph_name]]
      object@graphs[[new.graph]]@x <- object@graphs[[graph_name]]@x ^ st.res$powerEstimate
      miko_message(paste0("\tOptimal Power: ",  st.res$powerEstimate), time = F, verbose = verbose)

      miko_message("Running UMAP...", verbose = verbose)
      if (st.res$powerEstimate > 1){
        object <- tryCatch({
          object <- RunUMAP(object,  graph = new.graph, umap.method = "umap-learn", verbose = F) #, umap.method = "umap-learn"
          # object <- FindClusters(object, resolution = 1, verbose = verbose,  graph_name = new.graph)
        }, error = function(e){
          object <- RunUMAP(object, dims = 1:n_dim, n.neighbors = umap_knn, reduction  = "pca", verbose = F)
          # object <- FindClusters(object, resolution = 1, verbose = verbose)
          return(object)
        })
      } else {
        object <- RunUMAP(object, dims = 1:n_dim, n.neighbors = umap_knn, reduction  = "pca", verbose = F)
        # object <- FindClusters(object, resolution = 1, verbose = verbose)
      }


      snn.all <-  object@graphs[[new.graph]]@x
      snn.all[snn.all == 1] <- NA
      nn.threshold <- quantile(snn.all, binary_threshold, na.rm = T)

      new.nn <- paste0(new.graph, "_binary")
      object@graphs[[new.nn]] <- object@graphs[[new.graph]]
      object@graphs[[new.nn]]@x <- 0
      object@graphs[[new.nn]]@x[object@graphs[[graph_name]]@x > nn.threshold] <- 1


    })
  })

  return(
    list(
      object = object,
      st.res = st.res
    )
  )

}


#' Identify optimal cluster resolution of scale-free shared nearest neighbor network (SNN)
#'
#' Identify optimal cluster resolution of scale-free shared nearest neighbor network (SNN)
#'
#' @param object Seurat object
#' @param graph name of SNN graph in `object`.
#' @param target.purity Target purity to find resolution for.
#' @param start.res Initial resolution to begin optimization with. Default is 0.5.
#' @param step.size Step size between consecutive resolutions to test. Default is 0.05.
#' @param target.level Specify whether `target.purity` is achieved at "global" level [i.e., median(purity)] or at "cluster" level. Default is "global".
#' @param max.iter Maximum number of iterations. If `max.iter` is reached, non-optimal resolution may be returned.
#' @param verbose Show progress. Default is T.
#' @name SSNResolution
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{neighborPurity}}
#' @author Nicholas Mikolajewicz
#' @return optimal cluster resolution
#' @examples
#'
SSNResolution <- function(object, graph, target.purity = 0.7, start.res = 0.5, step.size = 0.05, target.level = "global", max.iter = 100, verbose = F){

  is.optimal <- F
  best.purity <- NA
  best.res <- NA

  i <- 0
  while(is.optimal == F){
    i <- i + 1
    initial.res <- start.res
    if (verbose)  miko_message(paste0("Current resolution: ", start.res))
    suppressWarnings({
      suppressMessages({
        object <- FindClusters(object, graph.name = graph, algorithm = 1, resolution = start.res, verbose = FALSE)
        object <- neighborPurity(object , graph = graph)
      })
    })

    if (target.level == "global"){
      m.purity <- signif(median(object$purity), 3)
    } else if (target.level == "cluster"){
      df.cluster.purity <- data.frame(cluster = object$seurat_clusters, purity = object$purity)
      df.cluster.purity.sum <- df.cluster.purity %>%
        dplyr::group_by(cluster) %>% dplyr::summarize(m.purity = median(purity, na.rm = T), .groups = "drop")
      m.purity <- signif(min(df.cluster.purity.sum$m.purity), 3)
    }


    if (m.purity > target.purity){
      if (verbose) miko_message(paste0("\tPurity: ", m.purity))
      best.purity <- m.purity
      best.res <- start.res
      start.res <- start.res + step.size
    } else if (m.purity == target.purity){
      best.purity <- m.purity
      best.res <- start.res
      is.optimal <- T
      if (verbose)  miko_message(paste0("Optimal Resolution: ", best.res , " (", scMiko::ulength(object$seurat_clusters), " clusters)")  )
    } else {
      if (verbose) miko_message("\tPurity: ", m.purity)

      if (is.na(best.purity)){
        start.res <-  0
      } else {
        if (verbose) message(paste0("\tToo low! Using resolution from prior iteration.") )
        is.optimal <- T
        if (verbose) message(paste0("Optimal Resolution: ", best.res , " (", scMiko::ulength(object$seurat_clusters), " clusters)") )
      }

    }


    if (i == max.iter){
      message("Failed to identify optimal resolution after ", max.iter, " iterations")
      is.optimal = T
    }

  }

  return(best.res)

}


#' Perform gene program discovery using SNN analysis
#'
#' Runs scale-free shared nearest neighbor network (SNN) analysis on subset of features specified in Seurat object.
#'
#' @param object Seurat object
#' @param features features to compute SNN on. If features are missing from scaled data, scaled data is recomputed.
#' @param scale_free Logical to enforce scale free topology. Default is T.
#' @param robust_pca Logical to run robust PCA (WARNING: computationally intensive, not recommended for large data). Default is F.
#' @param data_type Data type to compute SNN on.
#' \itemize{
#' \item "pearson" - pearson residuals for count data based on regularized negative binomial model.
#' \item "deviance" - deviance for count data based on multinomial null model (assumes each feature has constant rate).
#' }
#' @param reprocess_sct if `data_type` is "pearson", specify whether SCTransform is run (regardless whether features missing from existing scaled data or not). Default is F.
#' @param slot Slot to use.
#' \itemize{
#' \item "scale" - RECOMMENDED (default)
#' \item "data" - Not recommended and not tested extensively. Available for exploration. If specified, `data_type` is ignored.
#' }
#' @param batch_feature Variables to regress out. Default is NULL.
#' @param do_scale Whether to scale data (only if `slot` = "data")
#' @param do_center Whether to center data (only if `slot` = "data")
#' @param pca_var_explained Proportion of variance explained by PCA. Uses that top N PC components that explain `pca_var_explained` amount of variance. Default is 0.9.
#' @param weight_by_var Weight the feature embedding by the variance of each PC
#' @param umap_knn This determines the number of neighboring points used in local approximations of UMAP manifold structure. Larger values will result in more global structure being preserved at the loss of detailed local structure. In general this parameter should often be in the range 5 to 50. default is 10.
#' @param optimize_resolution Logical specifying whether to identify optimal clustering resolution. Optimal resolution identifying use target purity criteria. Default is T.
#' @param target_purity Target purity for identifying optimal cluster resolution. Default is 0.8.
#' @param step_size Step size between consecutive resolutions to test. Default is 0.05.
#' @param n_workers Number of workers for parallel implementation. Default is 1.
#' @param verbose Print progress. Default is T.
#' @name runSSN
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @import foreach parallel
#' @seealso \code{\link{findNetworkFeatures}} for finding features, \code{\link{SCTransform}} for gene count normalization and scaling, \code{\link{nullResiduals}} for deviance calculations, \code{\link{scaleFreeNet}} for scale-free topology transform.
#' @author Nicholas Mikolajewicz
#' @return Cell x Gene Seurat object, with gene-centric UMAP embedding and associated gene programs
#' @examples
#'
runSSN <- function(object, features, scale_free = T, robust_pca = F, data_type = c("pearson", "deviance"), reprocess_sct = F, slot = c("scale", "data"), batch_feature = NULL, do_scale = F, do_center = F, pca_var_explained = 0.9, weight_by_var = F, umap_knn = 10, optimize_resolution = T,
                   target_purity = 0.8, step_size = 0.05, n_workers = 1, verbose = T){

  require(parallel)
  require(foreach)

  if (!(class(object) %in% "Seurat")) stop("'object' is not Seurat object")
  which.rep.common <- features
  which.rep.common <- which.rep.common[which.rep.common %in% rownames(object)]
  which.data = slot

  if (robust_pca){
    generalPCA <- runRPCA
  } else {
    generalPCA <- RunPCA
  }


  # prep expression matrices
  miko_message("Preparing expression matrix...", verbose = verbose)
  if (which.data == "data"){
    expr.mat <- object@assays[[DefaultAssay(object)]]@data
  } else if (which.data == "scale"){

    all.required.genes <- unique(c(which.rep.common)) #which.TF


    # get SCTransform-scaled expression values
    if (data_type == "pearson"){

      if ((!("SCT" %in% names(object@assays))) | reprocess_sct){
        if ( "percent.mt" %in% colnames(object@meta.data)){
          v2r <- "percent.mt"
        } else {
          v2r <-NULL
        }

        if ((!is.null(batch_feature)) && batch_feature %in% colnames(object@meta.data)) v2r <- unique(c(v2r, batch_feature))

        object <- Seurat::SCTransform(object = object, method = "glmGamPoi", verbose = verbose,
                                      assay = "RNA", conserve.memory = T,
                                      vars.to.regress = v2r, variable.features.rv.th = 1.3,
                                      variable.features.n = 2000)
      }


      suppressMessages({
        suppressWarnings({


          try({
            miko_message("Calculating residuals...", verbose = verbose)
            available.genes <- rownames(object)
            available.genes <- available.genes[available.genes %in% all.required.genes]
            object@assays[[DefaultAssay(object)]]@var.features <- all.required.genes

            object <- GetResidual(object = object, features = all.required.genes, assay = "SCT", verbose = verbose,
                                  clip.range = c(-sqrt(x = ncol(x = object)/30), sqrt(x = ncol(x = object)/30)), replace.value = T)
          }, silent = T)

          if (!all(available.genes %in% rownames(object@assays[["SCT"]]@scale.data))){
            miko_message("Attempt to calculate residuals failed. Attempting alternative approach...", verbose = verbose)

            try({

              pmt.present <- "percent.mt" %in% colnames(object@meta.data)
              if (pmt.present){
                var2reg <- "percent.mt"
              } else {
                var2reg <- NULL
              }

              if (batch_feature %in% colnames(object@meta.data)) var2reg <- unique(c(var2reg, batch_feature))

              all.assays <- names(object@assays)
              if ("RNA" %in% all.assays){
                assay.use <- "RNA"
              } else {
                assay.use <- DefaultAssay(object)
              }

              so.query2 <- CreateSeuratObject(counts = object@assays[[assay.use]]@counts[rownames(object@assays[[assay.use]]@counts) %in%
                                                                                           all.required.genes, ],
                                              meta.data = object@meta.data)

              miko_message("Running SCTransform...", verbose = verbose)
              so.query2 <- SCTransform(object = so.query2, method = "glmGamPoi", verbose= verbose, return.only.var.genes	= F,
                                       vars.to.regress = var2reg, assay = assay.use)

              object@assays[["SCT"]]@scale.data <- so.query2@assays[["SCT"]]@scale.data
              object@assays[["SCT"]]@var.features <- rownames(object@assays[["SCT"]]@scale.data)

              rm(so.query2)
              invisible({gc()})

              if (all(available.genes %in% rownames(object@assays[["SCT"]]@scale.data))){
                miko_message(paste0("All ", length(available.genes), "/", length(available.genes), " residuals succesfully calculated. Proceeding with network generation."), verbose = verbose)
              } else {
                miko_message(paste0("Subset of ", length(rownames(object@assays[["SCT"]]@scale.data)), "/", length(available.genes), " residuals calculated. Proceeding with network generation."), verbose = verbose)
              }


            }, silent = T)

          }

        })
      })

      expr.mat <- object@assays[[DefaultAssay(object)]]@scale.data
      # parameter.list$general.do.scale <- parameter.list$general.do.center <-  F
      do_scale <- do_center <- F

    } else if (data_type == "deviance"){

      miko_message("Calculating deviance...", verbose = verbose)

      m <-  GetAssayData(object, slot = "data", assay = "SCT")
      m <- m[rownames(m) %in% all.required.genes, ]
      # so[["DEV"]] <- CreateAssayObject(counts = m)
      m <- as.matrix(m)

      # if ("batch.feature" %in% names(parameter.list)){
      if (batch_feature  %in% colnames(object@meta.data)){
        batch.factor <- as.factor(object@meta.data[ ,batch_feature])
      } else {
        batch.factor <- NULL
      }
      # } else {
      #   batch.factor <- NULL
      # }

      # system.time({
      m.deviance <- scry::nullResiduals(
        object = m,
        fam = c("binomial"),
        type = c("deviance"),
        batch = batch.factor
      )

      m.deviance[is.na(m.deviance)] <- 0
      rownames(m.deviance) <- rownames(m);
      colnames(m.deviance) <- colnames(m)
      invisible({gc()})

      # parameter.list$clip.range
      expr.mat <- object@assays[[DefaultAssay(object)]]@scale.data <- m.deviance
      rm(m.deviance);
      invisible({gc()})

      # parameter.list$general.do.scale <- parameter.list$general.do.center <-  F
      do_scale <- do_center <- F
    }


  }

  e.mat.success <- F
  try({
    # expr.mat.tf <- t(expr.mat[rownames(expr.mat) %in% which.TF, ])
    expr.mat <- t(expr.mat[rownames(expr.mat) %in% which.rep.common, ])
    e.mat.success <- T
  }, silent = T)

  if (!e.mat.success){
    # expr.mat.tf <- t(as.matrix(expr.mat[rownames(expr.mat) %in% which.TF, ]))
    expr.mat <- t(as.matrix(expr.mat[rownames(expr.mat) %in% which.rep.common, ]))
  }

  # general network
  miko_message(paste0("Generating transcriptional network using ", ncol(expr.mat), " genes..."))
  suppressMessages(suppressWarnings({
    so.gene <- CreateSeuratObject(expr.mat)
    so.gene <- so.gene[ ,colnames(so.gene) %in% which.rep.common]
    so.gene <- ScaleData(so.gene, do.scale = do_scale, do.center = do_center, verbose = verbose)
    if (robust_pca){
      so.gene <- runRPCA(so.gene, features = rownames(so.gene), verbose = F, reduction.key = "PC_", reduction.name = "pca")
    } else {


      so.gene <- RunPCA(so.gene, weight.by.var = weight_by_var, features = rownames(so.gene), verbose = F)
    }

    df.prop <- propVarPCA(so.gene)
    n.dim <- min(which(df.prop$pc.cum_sum > pca_var_explained))
    which.red <- "pca"
    so.gene <- FindNeighbors(so.gene, dims = 1:n.dim, reduction = which.red, verbose = verbose)


    if (scale_free){

      sfn.res <- scaleFreeNet(object = so.gene, graph_name = "RNA_snn", n_dim = n.dim, umap_knn = umap_knn, n_workers = n_workers)
      so.gene <- sfn.res$object
      sfn.res$object <- NULL
      plt.sft <- cowplot::plot_grid(sfn.res$st.res$optimization.plot,
                                    sfn.res$st.res$distribution.plot[[as.character(sfn.res$st.res$powerEstimate)]])

      so.gene@misc$scale_free = sfn.res$st.res

    } else {
      so.gene <- RunUMAP(so.gene, dims = 1:n.dim,  n.neighbors =umap_knn, reduction  = "pca")
      # so.gene <- FindClusters(so.gene, resolution = 1, verbose = verbose)

      plt.sft <- NULL
    }


    if (optimize_resolution){


      suppressMessages({
        suppressWarnings({

          # general network ################

          miko_message("Identifying optimal resolution...")

          # if (!is.null(so.gene@commands[["RunUMAP.RNA.pca"]]@params[["graph"]]) &&
          #     so.gene@commands[["RunUMAP.RNA.pca"]]@params[["graph"]] == "RNA_snn_power"){
          # so.gene@graphs[["RNA_snn_power"]]
          if ("RNA_snn_power" %in% names(so.gene@graphs)){
            # gn <- "RNA_snn_power_binary"
            gn <- "RNA_snn_power"
          } else {
            gn <- "RNA_nn"
          }
          best.res <- SSNResolution(object = so.gene, graph = gn,
                                    target.purity = target_purity, start.res = 1, verbose = verbose,
                                    step.size = step_size, target.level = "global", max.iter = 100)

          miko_message(paste0("Resolution of ", best.res, " achieves ", target_purity, " purity."), verbose = verbose)

          best.reg.general <- best.res
          # if (!is.null(so.gene@commands[["RunUMAP.RNA.pca"]]@params[["graph"]]) &&
          #     so.gene@commands[["RunUMAP.RNA.pca"]]@params[["graph"]] == "RNA_snn_power"){
          if ("RNA_snn_power" %in% names(so.gene@graphs)){
            so.gene <- FindClusters(so.gene, graph.name = gn,
                                    resolution = best.res, modularity.fxn = 1, algorithm =1, verbose = verbose)
          } else {
            so.gene <- FindClusters(so.gene,resolution = best.res,
                                    modularity.fxn = 1, algorithm =1, verbose = verbose)
          }


          # if (parameter.list$print.inline){
          #   cluster.UMAP(so.gene)
          # }

        })
      })


    }


  }))


  return(so.gene)

}


#' Get significant ICA genes
#'
#' Given feature loadings obtained from independent component analysis (ICA), significant genes for each component are identified.
#'
#' @param feature.loading ICA feature loadings
#' @param fdr.cutoff FDR threshold. Default is 0.0001.
#' @param local.fdr Whether to use local FDR. Default is T.
#' @param assert.positive.skew Whether to ensure that max(|feature.loading|) is positive. Default is T.
#' @param only.pos Whether significant genes are positively loaded only. Default is T.
#' @name getICAGenes
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{runICA}} for independent component analysis.
#' @author Nicholas Mikolajewicz
#' @return named list of ICA gene programs
#' @examples
#'
getICAGenes <- function(feature.loading, fdr.cutoff = 0.0001, local.fdr = T, assert.positive.skew = T, only.pos = T){

  ica.kme <- feature.loading
  ica.kme <- apply(ica.kme, 2, function(x) (x)/sd(x, na.rm = T))

  if (assert.positive.skew){
    median.loading <- abs(apply(ica.kme, 2, min))
    mean.loading <- abs(apply(ica.kme, 2, max))

    pos.skew <- mean.loading > median.loading

    reflection.vector <- pos.skew*1
    reflection.vector[reflection.vector == 0] <- -1

    for (i in 1:length(reflection.vector)){
      ica.kme[,i] <- ica.kme[,i]*reflection.vector[i]
    }
  }

  # compute FDR
  if (local.fdr){
    ica.p <- apply(ica.kme, 2,
                   function(x) fdrtool::fdrtool(x, statistic=c("normal"),
                                                plot=F, color.figure=F, verbose=F,
                                                cutoff.method=c("fndr"),
                                                pct0=0.75)$lfdr )

  } else {
    ica.p <- apply(ica.kme, 2,
                   function(x) fdrtool::fdrtool(x, statistic=c("normal"),
                                                plot=F, color.figure=F, verbose=F,
                                                cutoff.method=c("fndr"),
                                                pct0=0.75)$qval  )
  }

  rownames(ica.p) <- rownames(ica.kme)

  # count number of sig module genes
  if (only.pos){
    ica.p[ica.kme < 0] <- fdr.cutoff+1
  }

  # get module genes
  module.genes <- apply(ica.p, 2, function(x) rownames(ica.p)[x < fdr.cutoff])
  module.size <- unlist(lapply(module.genes, length))
  ica.module.genes <- module.genes[module.size > 0]

  return(ica.module.genes)
}


#' Run Independent Component Analysis on gene expression
#'
#' Run fastica algorithm from the ica package for ICA dimensionality reduction. Wrapper for Seurat's RunICA function, with additional significant gene identification step.
#'
#' @param object Seurat object
#' @param assay name of assay to use for ICA analysis. Expression data from 'scale.data' slot will be used.
#' @param features Features to compute ICA on. If not specified, all features present in `object` are used (not recommended, ICA is computationally expensive).
#' @param max_cells Max number of cells to run ICA on. If number of cells in `object` exceeded `max_cells`, object is subsampled prior to running ICA. Default is 20000.
#' @param verbose Print progress. Default is T.
#' @param ... Additional parameters passed to Seurat::RunICA(...)
#' @name runICA
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{RunICA}} for Seurat's independent component analysis, \code{\link{getICAGenes}} for significant ICA gene identification.
#' @author Nicholas Mikolajewicz
#' @return seurat object with significant genes stored in "misc" slot of ICA reduction slot.
#' @examples
#'
runICA <- function(object, assay = DefaultAssay(object), features = NULL, max_cells = 20000, verbose = T, ...){

  miko_message("Preparing expression data...", verbose = verbose)

  if (is.null(features)) features <- rownames(object)

  features_av <- features[features %in% rownames(object@assays[[assay]]@scale.data)]
  prop.av <- signif(100 * length(features_av)/length(features), 3)

  if (prop.av < 100){
    miko_message(length(features_av), "/", length(features) , " (", prop.av, "%) genes are available in 'scale.data'. Attempting to process remaining features...")
  } else {
    miko_message(length(features_av), "/", length(features) , " (", prop.av, "%) genes are available in 'scale.data'.")
  }


  if (prop.av < 100){
    # miko_message("Preprocessing expression data for missing features in 'scale.data'...")
    if (assay == "SCT"){
      miko_message("Calculating residuals using 'SCT' assay...")
      object <- GetResidual(object = object, features = features, verbose = verbose, assay = assay)
    } else if (assay == "RNA"){
      miko_message("Scaling features using 'RNA' assay...")
      object <- ScaleData(object = object, features = features, verbose = verbose, assay = assay)
    }

  }

  features_av <- features[features %in% rownames(object@assays[[assay]]@scale.data)]
  object@assays[[assay]]@var.features <- features_av

  miko_message(paste0("Running ICA on ", length(features_av), " features..."), verbose = verbose)

  if (ncol(object) > max_cells){
    set.seed(1023)
    so.ica <- RunICA(object[ ,sample(colnames(object), max_cells)],
                     features = features_av, verbose = F, ...)
  } else {
    so.ica <- RunICA(object, features = features_av, verbose = F, ...)
  }

  object@reductions$ica <- so.ica@reductions[["ica"]]
  object@reductions$ica@misc$meta.data <- so.ica@meta.data
  object@reductions$ica@misc$genes <- getICAGenes(object@reductions[["ica"]]@feature.loadings)
  object@reductions$ica@misc$commands <- so.ica@commands




  return(object)

}

#' Perform non-negative matrix factorization (NMF)
#'
#' Run non-negative matrix factorization(NMF or NNMF) using sequential coordinate-wise descent or multiplicative updates on Seurat object.
#'
#' @param object Seurat object
#' @param assay assay. Default is DefaultAssay(object).
#' @param k Number of NMF gene programs Default is 6.
#' @param raster rasterize output plot. Default is F.
#' @param n.threads Number of threads to use when running NMF. Default is 2.
#' @param features features to run NMF on.
#' @param feature.min.pct minimum expressing fraction for feature to run NMF on. Default is 0. Ignored if `features` are specified.
#' @param max.iter Maximum iteration of alternating NNLS solutions to H and W
#' @param gene.cutoff Feature loading cutoff threshold [0,1]. Default = 0.5. Ignored if gene.n is specified.
#' @param gene.n number of genes to return per gene program
#' @param sample.name sample name. Default is NULL.
#' @param show.top.n number of enrichment terms to show in summary plots. Ignoried if do.enrichment = F.
#' @param pathway.db pathway database to use for enrichment analysis. Options are "GO" or "Bader". Default is "Bader". Ignored if do.enrichment = F.
#' @param do.enrichment Whether to run enrichment analysis.
#' @param verbose Print progress. Default is TRUE.
#' @param ... additional arguments passed to NNLM::nnmf(...)
#' @name runNMF
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{nnmf}}
#' @author Nicholas Mikolajewicz
#' @return Seurat object with NMF results in reduction slot. Program genes are stored in "misc" slot of NMF reduction slot.
#' @examples
#'
runNMF <- function(object, assay = DefaultAssay(object), k = 6, raster = F, n.threads = 2, features = NULL, feature.min.pct = 0, max.iter = 50, gene.cutoff = 0.5, gene.n = 50,sample.name = NULL, show.top.n = 10,  pathway.db = "Bader", do.enrichment = T, verbose = T, ...){

  stopifnot("Seurat" %in% class(object))
  require("NNLM")

  miko_message("Preparing expression matrix...", verbose = verbose)

  DefaultAssay(object) <- assay

  if (!is.null(features)){
    expr.genes <- features
  } else if (is.null(features)){
    if (feature.min.pct >= 0 & feature.min.pct < 1){
      expr.genes <-  getExpressedGenes(object = object, min.pct = feature.min.pct)
    } else {
      expr.genes <- rownames(object)
    }
  }

  expr.genes <- unique(expr.genes)

  scale_success <- F
  try({

    expr.mat <- object@assays[[assay]]@scale.data
    if (!all(expr.genes %in% rownames(expr.mat))){
      object <- GetResidual(object = object, features = expr.genes, verbose = verbose, assay = assay)
    }
    scale_success <- T
  }, silent = T)

  if (scale_success){
    miko_message("Using 'scale.data' slot for NMF....", verbose = verbose)
    expr.mat <- object@assays[[assay]]@scale.data
    miko_message("Truncating negative values at zero....", verbose = verbose)

    expr.mat[expr.mat < 0] <- 0
  } else {
    miko_message("Using 'data' slot for NMF....", verbose = verbose)
    expr.mat <- object@assays[[assay]]@data
    expr.mat[expr.mat < 0] <- 0
  }

  if ("dgCMatrix" %in% class(expr.mat)) expr.mat <- as.matrix(expr.mat)
  expr.mat <- expr.mat[rownames(expr.mat) %in% expr.genes, ]


  miko_message("Running NMF...", verbose = verbose)
  set.seed(1023)
  nmf.opt  <- nnmf(expr.mat, k = k,
                   max.iter = max.iter,
                   # method = "scd",
                   loss = "mse",
                   rel.tol = 1e-4,
                   n.threads = n.threads,  #n.threads
                   verbose = F,
                   ...);

  nmf_name <- paste0("nmf_k", k)
  #
  object[[nmf_name]] <- CreateDimReducObject(loading = (nmf.opt[["W"]]),
                                             # key = "test",
                                             key = paste0(gsub("_", "", nmf_name), "_"),
                                             embeddings = t(nmf.opt[["H"]]),
                                             assay = DefaultAssay(object))


  getNMFGenes.dev <- function(feature.loading, norm.cutoff = 0.5, n.cutoff = NA){

    if (!("matrix" %in% class(feature.loading))) stop("feature.loading input is not a matrix")

    nmf.kme <- t(feature.loading)
    nmf.kme <- (apply(nmf.kme, 2, function(x) ((x^2)/sum(x^2))))

    # get module genes
    if (!is.na(n.cutoff)){
      module.genes <-  apply(nmf.kme, 1, function(x) colnames(nmf.kme)[order(-x)][1:n.cutoff])
      module.genes <-  wideDF2namedList(module.genes)
    } else {
      module.genes <-  apply(nmf.kme, 1, function(x) colnames(nmf.kme)[x>norm.cutoff])
    }

    module.size <- unlist(lapply(module.genes, length))
    nmf.module.genes <- module.genes[module.size > 0]

    return(nmf.module.genes)
  }

  # object[[nmf_name]]

  miko_message("Getting gene program features...", verbose = verbose)
  nmf.gene <- getNMFGenes.dev(feature.loading = nmf.opt$W, norm.cutoff = gene.cutoff, n.cutoff = gene.n)
  # df.cell.load <- data.frame(t(nmf.opt$H))
  df.cell.load <- data.frame(object@reductions[[nmf_name]]@cell.embeddings)
  colnames(df.cell.load) <- names(nmf.gene) <- paste0("NMF", seq(1, length(nmf.gene)))
  df.cell.load$var <- rownames(df.cell.load)

  umap_success <- F
  try({
    df.cell.umap <- getUMAP(object = object)[["df.umap"]]
    df.cell.umap <- merge(df.cell.umap, df.cell.load, by = "var")
    umap_success <- T
  }, silent = T)

  if (!umap_success) df.cell.umap <- NULL


  if (do.enrichment){

    miko_message("Running pathway enrichment...", verbose = verbose)
    nmf.res.hg <- runHG(gene.list =  nmf.gene, gene.universe = rownames(object), species = detectSpecies(object), pathway.db = pathway.db, verbose = verbose)
    nmf.sum.hg <- summarizeHG(nmf.res.hg, show.n = show.top.n)

    miko_message("Generating plots...", verbose = verbose)

    all.plots <- list()
    for (i in 1:length(nmf.gene)){

      df.cell.umap$score <- snip(x = df.cell.umap[ ,names(nmf.gene)[i]], lower.quantile = 0.001, upper.quantile = 0.999)
      if (raster){
        plt.nmf.load <- df.cell.umap %>%
          ggplot(aes(x = x, y = y, color = score)) +
          scattermore::geom_scattermore()
        # geom_point(size = autoPointSize(n.points = nrow(df.cell.umap))) +
        viridis::scale_color_viridis() +
          labs(title = paste0(names(nmf.gene)[i], " (", length(nmf.gene[[i]]), " genes)"),
               subtitle = sample.name, x = "UMAP 1", y = "UMAP 2", color = "Activity") +
          theme_miko(legend = T)
      } else {
        plt.nmf.load <- df.cell.umap %>%
          ggplot(aes(x = x, y = y, color = score)) +
          geom_point(size = autoPointSize(n.points = nrow(df.cell.umap))) +
          viridis::scale_color_viridis() +
          labs(title = paste0(names(nmf.gene)[i], " (", length(nmf.gene[[i]]), " genes)"),
               subtitle = sample.name, x = "UMAP 1", y = "UMAP 2", color = "Activity") +
          theme_miko(legend = T)
      }


      plt.nmf.load.combo <- cowplot::plot_grid(plt.nmf.load, nmf.sum.hg$plots[[names(nmf.gene)[i]]])

      # print(plt.nmf.load.combo)
      all.plots[[names(nmf.gene)[i]]] <- plt.nmf.load.combo
      # print(plt.nmf.load.combo)



    }


  } else {
    all.plots <- NULL
    nmf.res.hg <- NULL
    # df.cell.umap <- NULL
  }

  object@reductions[[nmf_name]]@misc$plots = all.plots
  object@reductions[[nmf_name]]@misc$genes = nmf.gene
  object@reductions[[nmf_name]]@misc$df.cell.umap = df.cell.umap
  object@reductions[[nmf_name]]@misc$hg.results = nmf.res.hg
  object@reductions[[nmf_name]]@misc$nnmf_output = nnmf

  return(object)

}

#' Consolidate several NMF reduction objects into single NMF reduction object within Seurat object.
#'
#' Consolidate several NMF reduction objects into single NMF reduction object within Seurat object.
#'
#' @param object Seurat object containing several NMF reductions (e.g., each NMF with different k parameter)
#' @param reduction_key Reduction key used to identify reduction objects for consolidation. Default is "nmf".
#' @name consolidateNMF
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{runNMF}}
#' @author Nicholas Mikolajewicz
#' @return Seurat object with with several NMF reductions consolidated into "nmf_all" reduction object.
#' @examples
#'
consolidateNMF <- function(object, reduction_key = "nmf"){

  knames <- names(object@reductions)[grepl( reduction_key, names(object@reductions))]

  df.umap.nmf <- NULL
  nmf.object <- NULL
  hmat_all <- wmat_all <- NULL
  genes_all <- list()
  umap_all <- NULL
  for (i in 1:length(knames)){

    k_cur <- gsub("nmf_", "", knames[i])
    res_cur <- object@reductions[[knames[i]]]
    nc <- ncol(res_cur@feature.loadings)
    # nmf_res_cur <- nmf_res[[i]]
    # nc <- ncol(nmf_res_cur[["nmf.results"]][["W"]])
    wmat <- res_cur@feature.loadings
    colnames(wmat) <- paste0(k_cur, "_NMF", seq(1, ncol(wmat)))
    wmat_all <- cbind(wmat_all, wmat)

    # hmat <- nmf_res_cur[["nmf.results"]][["H"]]
    hmat <- res_cur@cell.embeddings
    colnames(hmat) <- paste0(k_cur, "_NMF", seq(1, ncol(hmat)))
    hmat_all <- cbind(hmat_all, hmat)


    genes_cur <- res_cur@misc$genes
    names(genes_cur) <- paste0(k_cur, "_", names(genes_cur))
    genes_all <- c(genes_all, genes_cur)

    df.umap_cur <- res_cur@misc$df.cell.umap
    try({
      df.umap_cur <- df.umap_cur %>% dplyr::select(-c("score"))
    }, silent = T)

    colnames(df.umap_cur)[grepl("NMF", colnames(df.umap_cur))] <- paste0(k_cur, "_", colnames(df.umap_cur)[grepl("NMF", colnames(df.umap_cur))])
    if (is.null(umap_all)){
      umap_all <- df.umap_cur
    } else {
      umap_all <- merge(umap_all, df.umap_cur)
    }
  }

  object[["nmf_all"]] <- CreateDimReducObject(loading = wmat_all,
                                              # key = "test",
                                              key = paste0("NMF_"),
                                              embeddings = hmat_all,
                                              assay = DefaultAssay(object))

  colnames(object@reductions[["nmf_all"]]@feature.loadings) <- colnames(wmat_all)
  colnames(object@reductions[["nmf_all"]]@cell.embeddings) <- colnames(hmat_all)

  object@reductions[["nmf_all"]]@misc$genes <- genes_all
  object@reductions[["nmf_all"]]@misc$df.cell.umap <- umap_all
  object@reductions[["nmf_all"]]@misc$meta.data <- object@meta.data

  return(object)

}


#' Named list of cells grouped by meta feature from Seurat object
#'
#' Returns named list of cells grouped by meta feature from Seurat object.
#'
#' @param object Seurat object
#' @param group grouping feature (must be present in `object` meta data). Default is "seurat_clusters".
#' @param is.num whether `group` is a numerical feature.
#' @param prefix prefix added to each named entry in list. Default is "".
#' @name group2list
#' @author Nicholas Mikolajewicz
#' @return Returns named list of cells grouped by meta feature from Seurat object.
#' @examples
#'
group2list <- function(object, group = "seurat_clusters", is.num = F, prefix = ""){

  if (!("Seurat" %in% class(object))) stop("'object' is not a Seurat object")
  df.meta <- object@meta.data
  if (!(group %in% colnames(df.meta))) stop(paste0("'", group, "' is not present in meta data"))

  u.clust <- unique(df.meta[,group])
  if (is.num){
    u.clust <- as.numeric(as.character(u.clust))
    u.clust <- u.clust[order(u.clust)]
  } else {
    u.clust <- as.character(u.clust)
    u.clust <- u.clust[order(u.clust)]
  }

  group.list <- list()
  for (i in 1:length(u.clust)){
    group.list[[paste0(prefix, u.clust[i])]] <- rownames(df.meta)[df.meta[, group] %in% u.clust[i]]
  }

  return(group.list)
}


#' Identify and (optionally) prune gene program features in scale-free shared nearest neighbor network (SSN)
#'
#' Identify and (optionally) prune gene program features in scale-free shared nearest neighbor network (SSN)
#'
#' @param object Seurat object objected from SSN analysis.
#' @param graph name of shared nearest neighbor (SNN) graph used for SSN analysis.
#' @param prune.threshold threshold used for pruning features from SSN graph (features with low connectivity are preferentially pruned). Default is 0.1.
#' @param return.df whether to return data.frame instead of named list.
#' @name pruneSSN
#' @author Nicholas Mikolajewicz
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @return Returns named list of SSN gene program features.
#' @seealso \code{\link{runSSN}} for SSN analysis
#' @examples
#'
pruneSSN <- function(object, graph = "RNA_snn_power", prune.threshold = 0.1, return.df = F){

  valid.pt <- prune.threshold > 0 & prune.threshold< 1

  if (!(graph %in% names(object@graphs))) stop(paste0(graph, " is not in present in 'object'"))
  mod.list <- group2list(object, is.num = T, prefix = "m")
  if (valid.pt){

    pthresh <-prune.threshold
    mod.list <- mod.list[unlist(lapply(mod.list, function(x) length(x) > 0))]
    snn.graph <- as.matrix(object@graphs[[graph]])

    df.con.l2 <- NULL
    mod.list_v2 <- list()
    for (i in 1:length(mod.list)){
      snn_cur <- snn.graph[mod.list[[i]], mod.list[[i]]]
      df.connectivity <- getConnectivity(snn_cur, rownames(snn_cur))
      df.connectivity$wi_l2 <- (df.connectivity$wi^2)/sum(df.connectivity$wi^2)
      df.connectivity$wi_l2 <- df.connectivity$wi_l2 / max(df.connectivity$wi_l2)
      df.con.l2 <- bind_rows(df.con.l2, df.connectivity)
      mod.list_v2[[names(mod.list)[i]]] <- df.connectivity$genes[df.connectivity$wi_l2 > pthresh]
    }

    mod.list <- mod.list_v2
  }

  if (return.df){
    return(df.con.l2)
  } else {
    return(mod.list)
  }

}




#' SSN connectivity plot
#'
#' generates SSN connectivity plot used to visualize transcriptomic gene network.
#'
#' @param gene.object Seurat object (cell x gene) obtained from SSN analysis.
#' @param gene.list named list of gene program features/genes. If omitted, feature-annotated connectivity plot is not generated.
#' @param quantile_threshold Quantile threshold for visualized SSN graph edges. Default is 0.9.
#' @param raster_dpi prefix added to each named entry in list. Default is "".
#' @param edge.alpha alpha [0,1] parameter (i.e., transparency) for network edges. Default is 0.015. Use larger values for less dense networks.
#' @param edge.color color of network edges. Default is 'black'.
#' @param node.color color of node edges. Default is 'grey'.
#' @param node.alpha alpha [0,1] parameter for nodes. Default is 1 (obaque).
#' @param node.weights scale node size by connectivity weights. Default is F.
#' @param node.size.min minimal node size. Ignored if node.weights = F. Default is 1.
#' @param node.size.max maximal node size. Ignored if node.weights = F. Default is 5.
#' @param do.label Show gene program IDs  on network graph. Default is T.
#' @param label.size Size of gene program IDs on network graph. Default is 5.
#' @param verbose Print progress. Default is T.
#' @name SSNConnectivity
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @author Nicholas Mikolajewicz
#' @seealso \code{\link{runSSN}} for SSN analysis (gene.object), \code{\link{pruneSSN}} for gene program features (gene.list)
#' @return Returns 2 variants of SSN connectivity plot along with data.frame used to generate plots.
#' @examples
#'
SSNConnectivity <- function(gene.object, gene.list = NULL, quantile_threshold = 0.9, raster_dpi = 200,edge.alpha = 0.015,
                            edge.color = "black", node.color = "grey", node.alpha = 1,node.weights = F, node.size.min = 1, node.size.max = 5, do.label = T, label.size = 5,  verbose = T){

  # get connectivity ########################################
  stopifnot("'gene.object' is not a Seurat object" = ("Seurat" %in% class(gene.object)))
  if ("RNA_snn_power" %in% names(gene.object@graphs)){
    graph <- "RNA_snn_power"
  } else {
    graph <- "RNA_snn"
  }
  snn.graph <- as.matrix(gene.object@graphs[[graph]])

  df.connectivity <- getConnectivity(snn.graph, rownames(snn.graph))


  graph.threshold <- quantile(snn.graph[snn.graph > 0], quantile_threshold)
  miko_message(paste0("Edge threshold set at SNN = ", signif(graph.threshold,3), " (", quantile_threshold, " quantile)"), verbose = verbose)
  df.graph.edges <- NULL

  miko_message("Getting network edges...", verbose = verbose)
  knn.weights <- pbapply::pbapply(snn.graph, 2, function(x) x[x>=graph.threshold])
  df.kw <- bind_rows(knn.weights)
  df.kw <- as.data.frame(df.kw)
  rownames(df.kw) <- names(knn.weights)
  df.kw.col <- colnames(df.kw)
  df.kw$start <- names(knn.weights)
  df.graph.edges <- pivot_longer(df.kw, cols = df.kw.col)
  df.graph.edges <- df.graph.edges[complete.cases(df.graph.edges), ]
  colnames(df.graph.edges) <- c("end", "start", "weights")

  # remove duplicate entries
  df.graph.edges <- df.graph.edges %>% dplyr::filter(start != end)
  df.graph.edges$order <- df.graph.edges$start < df.graph.edges$end
  df.graph.edges$start_hold <- df.graph.edges$start
  df.graph.edges$end_hold <- df.graph.edges$end
  df.graph.edges$start[df.graph.edges$order] <- df.graph.edges$end_hold[df.graph.edges$order]
  df.graph.edges$end[df.graph.edges$order] <- df.graph.edges$start_hold[df.graph.edges$order]
  df.graph.edges$dup <- duplicated(paste0(df.graph.edges$start, "_", df.graph.edges$end))
  df.graph.edges <- df.graph.edges %>% dplyr::filter(!dup) %>% dplyr::select(c("end", "start", "weights"))

  snnUMAP.list <- getUMAP(gene.object, umap.key = "umap", node.type = "text")

  # snnUMAP.list <- getUMAP(so.gene.meso, umap.key = "wnn.umap", node.type = "text")
  df.snn.umap <- snnUMAP.list$df.umap
  df.snn.umap.x <- df.snn.umap
  colnames(df.snn.umap.x) <- c("x.start", "y.start", "start", "cluster")
  df.snn.umap.y <- df.snn.umap
  colnames(df.snn.umap.y) <- c("x.end", "y.end", "end", "cluster")
  df.ux <- merge(df.graph.edges, df.snn.umap.x, by = "start")
  df.uy <- merge(df.graph.edges, df.snn.umap.y, by = "end")
  df.umap.merge <- merge(df.ux, df.uy, by = c("start", "end"))

  df.snn.umap$genes <- df.snn.umap$var

  df.snn.umap <- merge(df.snn.umap, df.connectivity, by = "genes")

  plt.general.connectivitiy <- ggplot() +
    ggrastr::rasterise(geom_segment(data = df.umap.merge, aes(x = x.start, y = y.start,
                                                              xend = x.end, yend= y.end,
    ), color = edge.color, alpha =edge.alpha), dpi = raster_dpi)  +
    theme_miko(legend = T) +
    theme(
      panel.border = element_blank(),
      axis.text = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank()
    ) + xlab("") + ylab("")


  df.snn.umap$module <- df.snn.umap$seurat_clusters
  df.snn.umap$module <- paste0("m", df.snn.umap$module)
  df.snn.umap$module <- gsub("mother", "other", df.snn.umap$module)

  df.wnn.umap.sum <- df.snn.umap %>%
    dplyr::filter(seurat_clusters != "other") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarize(x.mean = median(x), y.mean = median(y), .groups = 'drop')

  # get color palettes for each gene module
  umod <- unique(as.numeric(as.character(gene.object@meta.data$seurat_clusters)))
  umod <- umod[!is.na(umod)]
  umod <- umod[order(umod)]
  col.pal <- categoricalColPal(labels = umod)
  col.pal <- c(col.pal, "lightgrey")
  names(col.pal) <- c(paste0("m", umod), "other")


  if (!is.null(gene.list)){

    if (node.weights){
      gp <- geom_point(aes(fill = module, label = gene, size = wi), pch = 21, color = node.color, alpha = node.alpha)
    } else {
      gp <- geom_point(aes(fill = module, label = gene), pch = 21, color = node.color, alpha = node.alpha)
    }

    df.snn.umap$gene <- df.snn.umap$genes
    plt.wnn.umap.connectivity.net.all <- df.snn.umap %>%
      dplyr::filter(var %in% unique(unlist(gene.list))) %>%
      ggplot(aes(x = x, y = y)) +
      ggrastr::rasterise(geom_segment(data = df.umap.merge, aes(x = x.start, y = y.start,
                                                                xend = x.end, yend= y.end,
      ), color = edge.color, alpha = edge.alpha), dpi = raster_dpi) +   #alpha = weights.x , alpha = 0.8
      gp +
      # geom_point(aes(fill = module, label = gene), pch = 21, color = node.color) +
      theme_miko(legend = T) +
      # geom_point(aes(fill = seurat_clusters), shape = 21, color = "grey", ...) + theme_miko(legend = T) +  #, size = wi
      scale_size_continuous(range = c(node.size.min, node.size.max)) +
      theme(
        panel.border = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_blank()
      ) + xlab("") + ylab("")  +
      labs(size = "Degree", alpha = "Weight", fill = "Gene Program") +
      scale_fill_manual(values = col.pal)

    if (do.label){
      # plt.wnn.umap.connectivity.net.all <- plt.wnn.umap.connectivity.net.all + ggrepel::geom_text_repel(data = df.wnn.umap.sum, aes(x = x.mean, y = y.mean, label = seurat_clusters), size = label.size)
      plt.wnn.umap.connectivity.net.all <- plt.wnn.umap.connectivity.net.all + geom_text(data = df.wnn.umap.sum, aes(x = x.mean, y = y.mean, label = seurat_clusters), size = label.size)
    }

    output.list <-   list(
      plot_edge = plt.general.connectivitiy,
      plot_network = plt.wnn.umap.connectivity.net.all,
      df.snn.umap = df.snn.umap,
      df.umap.merge = df.umap.merge
    )
  } else {
    output.list <-   list(
      plot_edge = plt.general.connectivitiy,
      # plot_network = plt.wnn.umap.connectivity.net.all,
      df.snn.umap = df.snn.umap,
      df.umap.merge = df.umap.merge
    )
  }


  return(
    output.list
  )
}


#' SSN connectivity plot
#'
#' generates SSN connectivity plot used to visualize transcriptomic gene network.
#'
#' @param cell.object Seurat object (gene x cell). Same as input to runSSN(object = cell.object, ...).
#' @param gene.object Seurat object (cell x gene) obtained from SSN analysis.
#' @param group_by meta feature used for grouping. Must be present in `cell.object`.
#' @param features feature to visualize expression.
#' @param max_cells maximum number of cells used. Default is 20000. If number of cells in `cell.object` exceeds `max_cells`, `cell.object` is downsampled.
#' @param connecitivity.plot "plot_edge" generated by SSNConnectivity function. If not specified, not edges are plotted in SSN graph plot.
#' @name SSNExpression
#' @author Nicholas Mikolajewicz
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{runSSN}} for SSN analysis (gene.object), \code{\link{pruneSSN}} for gene program features (gene.list), \code{\link{SSNConnectivity}} for connectivity plot.
#' @return Returns 2 variants of SSN connectivity plot along with data.frame used to generate plots.
#' @examples
#'
SSNExpression <- function(cell.object, gene.object, group_by = "seurat_clusters", features = NULL, max_cells = 20000, connectivity.plot = NULL){

  if (!("Seurat" %in% class(cell.object))) stop("'cell.object' is not a Seruat object")
  if (!("Seurat" %in% class(gene.object))) stop("'gene.object' is not a Seruat object")

  if (!is.null(features)){
    all.gene <- unique(unlist(features))
    all.gene <- all.gene[all.gene %in% colnames(gene.object)]
  } else {
    all.gene <- colnames(gene.object)
  }

  plt.dot <- DotPlot(cell.object,
                     features = all.gene, group.by = group_by)
  df.dot <- plt.dot[["data"]]

  if (group_by == "seurat_clusters"){
    u.clust <- as.numeric(as.character(unique(df.dot$id)))
    u.clust <- as.character(u.clust[order(u.clust)])
    prefix <- "C"
  } else {
    u.clust <- (as.character(unique(df.dot$id)))
    prefix <- ""
  }


  df.dot$gene <- df.dot$features.plot

  df.wnn.umap.cur <- getUMAP(gene.object)[["df.umap"]]
  df.wnn.umap.cur$gene <- df.wnn.umap.cur$var
  df.umap.sub <- getUMAP(cell.object, meta.features = group_by)[["df.umap"]]

  if (nrow(df.umap.sub) > max_cells){
    set.seed(1023)
    df.umap.sub <- df.umap.sub[sample(seq(1, nrow(df.umap.sub)), max_cells), ]
  }


  plt.clut.net.highlight.list <-  pbapply::pblapply(u.clust, function(x){

    uc_cur <- unique(unlist(unlist(x)))
    # print(uc_cur)
    df.dot.cur <- df.dot %>% dplyr::filter(id == uc_cur)
    df.wnn.umap.cur.cur <- df.wnn.umap.cur

    df.umap.sub.cur <- df.umap.sub
    df.umap.sub.cur$is.clust <- df.umap.sub.cur[ ,group_by] %in% uc_cur
    df.wnn.umap.cur.cur <- merge(df.wnn.umap.cur.cur, df.dot.cur, by = "gene")

    plt.clust.highlight <- df.umap.sub.cur %>%
      ggplot(aes(x = x, y = y, color = is.clust)) +
      geom_point(size = autoPointSize(nrow(df.umap.sub.cur))) +
      scale_color_manual(values= c("TRUE" = "tomato", "FALSE" = "grey")) +
      theme_miko()+ theme_void() + theme(legend.position = "none") +
      labs(x = "UMAP 1", y = "UMAP 2", title = paste0(prefix, uc_cur,  " Cell Population"), subtitle = " ")


    if (!is.null(connectivity.plot)){
      plt.clut.net <- connectivity.plot +

        geom_point(aes(x = x, y = y, color = avg.exp.scaled),
                   data = df.wnn.umap.cur.cur, size = autoPointSize(nrow(df.wnn.umap.cur.cur))) +
        scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
        theme_miko(legend = T) + theme_void()+
        scale_size(range = c(0, 5)) +
        labs(x = "UMAP 1", y = "UMAP 2", title = paste0(prefix, uc_cur,  " Network Expression"), subtitle = "Transcriptomic Network", color = "expr")

    } else {
      plt.clut.net <- df.wnn.umap.cur.cur %>%
        ggplot(aes(x = x, y = y, color = avg.exp.scaled)) +  #, size = pct.exp
        scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
        geom_point(size = autoPointSize(nrow(df.wnn.umap.cur.cur))) + theme_miko(legend = T) + theme_void()+
        scale_size(range = c(0, 5)) +
        labs(x = "UMAP 1", y = "UMAP 2", title = paste0("C", uc_cur,  " Network Expression"), subtitle = "Transcriptomic Network", color = "expr")
    }


    plt.clut.net.highlight <- cowplot::plot_grid(plt.clust.highlight, plt.clut.net)

    return(plt.clut.net.highlight)

  })

  names(plt.clut.net.highlight.list) <- paste0("m", u.clust)

  return(plt.clut.net.highlight.list)
}


#' Summarize SSN, ICA or NMF gene program/module analyses
#'
#' Summarize SSN, ICA or NMF gene program/module analyses. Expression heatmaps and gene program enrichments are generated to facilitate analysis and interpretation of gene programs
#'
#' @param cell.object Seurat object (gene x cell). Same as input to runSSN(object = cell.object, ...).
#' @param gene.object Seurat object (cell x gene) obtained from SSN analysis.
#' @param module.type Which gene program results to summarize
#' \itemize{
#' \item "ssn" - scale-free shared nearest neighbor network
#' \item "ica" - independent component analysis
#' \item "nmf" - non-negative matrix factorization
#' }
#' @param gene.list feature to visualize expression.
#' @param reduction specify which reduction slot use for summary. Ignored if `module.type` = "ssn".
#' @param connecitivity.plot "plot_edge" generated by SSNConnectivity function. If not specified, not edges are plotted in SSN graph plot.
#' @param group_by which meta feature to group by. Default is "seurat_clusters".
#' @param raster whether UMAPs are rasterized (recommended for large datasets).
#' @param show.n.pathways number of pathway annotations to show in each enrichment plot. Default is 10.
#' @param n.workers number of workers for parallel implementation of pathway enrichment analysis. Default is 1.
#' @param raster.threshold cell count threshold at which to switch to rasterized plots. Default is 10000.
#' @param verbose Print progress. Default is T.
#' @name summarizeModules
#' @author Nicholas Mikolajewicz
#' @references \url{https://nmikolajewicz.github.io/scMiko/articles/Module_Detection.html}
#' @seealso \code{\link{runSSN}} for SSN analysis (gene.object), \code{\link{pruneSSN}} for gene program features (gene.list), \code{\link{SSNConnectivity}} for connectivity plot.
#' @return List of summarize gene program results, including expression heatmaps and pathway enrichments
#' @examples
#'
summarizeModules <- function(cell.object, gene.object, module.type = c("ssn", "ica", "nmf"), gene.list = NULL,
                             reduction = NULL, connectivity_plot = NULL, group_by = "seurat_clusters", raster = F,
                             show.n.pathways = 10, n.workers = 1, raster.threshold = 10000, verbose = T){


  if (!(group_by %in% colnames(cell.object@meta.data))) stop("'group_by' is not in 'cell.object'")
  # cell to cluster mapping

  if (module.type == "ica"){
    miko_message("Summarizing ICA gene programs...", verbose = verbose)
    df.meta <- cell.object@reductions[[reduction]]@misc[["meta.data"]]
    cellClusterMap <- data.frame(cell = rownames(df.meta),
                                 cluster = df.meta[ ,group_by])
    cell.emb <- cell.object@reductions[[reduction]]@cell.embeddings
    feat.emb <- cell.object@reductions[[reduction]]@feature.loadings
    feat.emb <- data.frame(feat.emb)

    subtitle_label <- "x = IC programs; y = Cell Clusters; z = Cell Loading (Scaled)"
  } else if (module.type == "nmf"){
    miko_message("Summarizing NMF gene programs", verbose = verbose)
    df.meta <- cell.object@reductions[[reduction]]@misc[["meta.data"]]
    cellClusterMap <- data.frame(cell = rownames(df.meta),
                                 cluster = df.meta[ ,group_by])

    cell.emb <-  cell.object@reductions[[reduction]]@cell.embeddings
    feat.emb <- cell.object@reductions[[reduction]]@feature.loadings
    feat.emb <- data.frame(feat.emb)

    subtitle_label <- "x = NMF programs; y = Cell Clusters; z = Cell Loading (Scaled)"

  } else if (module.type == "ssn"){
    miko_message("Summarizing SSN gene programs...", verbose = verbose)
    if (is.null(gene.list)) stop("'gene.list' is required for summarizing SSN programs")

    miko_message("Scoring SSN gene programs", verbose = verbose)
    suppressMessages({
      suppressWarnings({

        df.meta <- cell.object@meta.data
        cellClusterMap <- data.frame(cell = rownames(df.meta),
                                     cluster = df.meta[ ,group_by])

        mod.res <-  runMS(
          object = cell.object,
          genelist = gene.list,
          assay = DefaultAssay(cell.object),
          size = autoPointSize(ncol(cell.object)),
          return.plots = F,
          verbose = F,
          winsorize.quantiles = c(0.01, 0.99)
        )

      })
    })

    if (!is.null(connectivity_plot)){
      df.con <- connectivity_plot[["plot_env"]][["df.snn.umap"]]
    } else {
      df.con <- NULL
    }

    subtitle_label <- "x = SSN programs; y = Cell Clusters; z = Gene Program Activity (Scaled)"

  }

  # cluster-level activity
  if (group_by == "seurat_clusters"){
    u.cluster <- as.numeric(as.character(unique(cellClusterMap$cluster)))
    u.cluster <- u.cluster[order(u.cluster)]
  } else {
    u.cluster <- (as.character(unique(cellClusterMap$cluster)))
  }

  # cluster-expression heatmaps
  if (module.type %in% c("ica", "nmf")){
    cluster.activity <- matrix(nrow = length(u.cluster), ncol = dim(cell.emb)[2])
    for (i in 1:(length(u.cluster))){
      which.cell <- cellClusterMap$cell[cellClusterMap$cluster == u.cluster[i]]
      for (j in 1:(dim(cell.emb)[2])){
        cluster.activity[i,j] <- mean(cell.emb[ rownames(cell.emb) %in% which.cell,j])
      }
    }
    rownames(cluster.activity) <- paste0("C", u.cluster)
    colnames(cluster.activity) <- colnames(cell.emb)

    plt.expr <- miko_heatmap(mat = cluster.activity,
                             color = colorRampPalette(rev(brewer.pal(n = 7, name =   "RdBu")))(100),
                             scale  = "column") +
      labs(
        title = "Gene Program Activities",
        subtitle = subtitle_label
      )

  } else if (module.type == "ssn"){

    mod.res[["data"]]$cluster <- cell.object@meta.data[, group_by]
    df.mod.wide <- mod.res[["data"]] %>%
      dplyr::select(grep("m", colnames(mod.res[["data"]]))) %>%
      dplyr::select(-c("class.ms"))

    cell.emb <-  df.mod.wide

    df.mod.wide$cluster <- cell.object@meta.data[, group_by]
    suppressMessages({
      suppressWarnings({
        df.mod.wide <- df.mod.wide %>%
          pivot_longer(cols = grep("m", colnames(df.mod.wide))) %>%
          dplyr::group_by(cluster, name) %>%
          dplyr::summarize(x = mean(value, na.rm = T)) %>%
          dplyr::mutate(cluster = paste0("c", cluster)) %>%
          pivot_wider(names_from = cluster,values_from = x)
        df.mod.wide <- col2rowname(df.mod.wide, "name")
        df.mod.wide <- as.data.frame(t(df.mod.wide))
      })
    })

    cluster.activity <- df.mod.wide
    plt.expr <- miko_heatmap(cluster.activity, scale ="column") +
      labs(title = "Gene Program Activities",
           subtitle = subtitle_label)

  }



  #
  # get ica genes
  if (module.type == "ica"){
    module.genes = cell.object@reductions[[reduction]]@misc$genes
    names(module.genes) <- gsub("_", "", names(module.genes))
  } else if (module.type == "nmf"){
    module.genes = cell.object@reductions[[reduction]]@misc$genes
  } else if (module.type == "ssn"){
    module.genes = gene.list
  }



  module.merge.genes <- module.genes

  miko_message("Performing functional annotation...", verbose = verbose)


  miko_message("\tGO Biological Processes...", verbose = verbose, time = F)
  suppressMessages({
    suppressWarnings({
      bp.res <- runHG(module.merge.genes, gene.universe = rownames(cell.object), n.workers = n.workers,
                      species = detectSpecies(cell.object), pathway.db = "GO", go.ontology = "BP", verbose = F)
    })
  })
  miko_message("\tGO Molecular Functions...", verbose = verbose, time = F)
  suppressMessages({
    suppressWarnings({
      mf.res <- runHG(module.merge.genes, gene.universe = rownames(cell.object),n.workers = n.workers,
                      species = detectSpecies(cell.object), pathway.db = "GO", go.ontology = "MF", verbose = F)
    })
  })
  miko_message("\tGO Cellular Compartments...", verbose = verbose, time = F)
  suppressMessages({
    suppressWarnings({
      cc.res <- runHG(module.merge.genes, gene.universe = rownames(cell.object), n.workers = n.workers,
                      species = detectSpecies(cell.object), pathway.db = "GO", go.ontology = "CC", verbose = F)    })
  })



  bp.sum <- summarizeHG(bp.res, show.n = show.n.pathways, pathway.name.size = 10)
  mf.sum <- summarizeHG(mf.res, show.n = show.n.pathways, pathway.name.size = 10)
  cc.sum <- summarizeHG(cc.res, show.n = show.n.pathways, pathway.name.size = 10)

  col.pal <- categoricalColPal(names(module.merge.genes))

  df.umap.gene <- getUMAP(gene.object)[["df.umap"]]
  df.umap.cell2 <- getUMAP(cell.object, meta.features = group_by)[["df.umap"]]
  df.umap.gene$gene <- df.umap.gene$var
  miko_message("Consolidating results...", verbose = verbose)

  all.names <- names(module.merge.genes)
  plt.list <- pbapply::pblapply(all.names,function(x){

    mod.name <- unlist(x[[1]])
    ngenes <- length(module.merge.genes[[mod.name]])

    if (module.type == "ica"){
      df.umap.cell <- df.umap.cell2

      ica.suc <- F
      try({
        df.umap.cell$cell.loading <- cell.emb[, gsub("IC", "IC_", mod.name)]
        ica.suc <- T
      }, silent = T)
      if (!ica.suc){
        df.cell.emb <- as.data.frame(cell.emb[, gsub("IC", "IC_", mod.name)])
        colnames(df.cell.emb) <- "cell.loading"
        df.cell.emb$cells <- rownames(df.cell.emb)
        df.umap.cell$cells <- rownames(df.umap.cell)
        df.umap.cell <- merge(df.umap.cell, df.cell.emb, all.y = T, by = "cells")
      }

      if (median(df.umap.cell$cell.loading) > mean(df.umap.cell$cell.loading)){
        df.umap.cell$cell.loading <- df.umap.cell$cell.loading*(-1)
      }
      color.label <- "ICA Loading"
      cell.label <- " cell loading"

    } else if (module.type == "nmf"){
      df.umap.cell <- df.umap.cell2

      nmf.suc <- F
      try({
        df.umap.cell$cell.loading <- cell.emb[, mod.name]
        nmf.suc <- T
      }, silent = T)
      if (!nmf.suc){
        df.cell.emb <- as.data.frame(cell.emb[, mod.name])
        colnames(df.cell.emb) <- "cell.loading"
        df.cell.emb$cells <- rownames(df.cell.emb)
        df.umap.cell$cells <- rownames(df.umap.cell)
        df.umap.cell <- merge(df.umap.cell, df.cell.emb, all.y = T, by = "cells")
      }

      color.label <- "NMF Loading"
      cell.label <- " cell loading"

    } else if (module.type == "ssn"){
      df.umap.cell <- df.umap.cell2
      df.umap.cell$cell.loading <- cell.emb[, mod.name]
      color.label <- "Gene Program Activity"
      cell.label <- " cell activity"

    }

    if (nrow(df.umap.cell) > raster.threshold) {
      gp <- scattermore::geom_scattermore
    } else {
      gp <- geom_point
    }

    if (raster) gp <- scattermore::geom_scattermore

    plt.cell <- df.umap.cell %>%
      ggplot(aes(x = x, y = y, color = cell.loading)) +
      gp(size = autoPointSize(nrow(df.umap.cell))) +
      theme_miko() + theme_void() + labs(title = "Cellular UMAP",
                                         subtitle = paste0(mod.name, cell.label)) +
      scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      theme(legend.position = "bottom") + labs(color = color.label)

    df.umap.cur <- df.umap.gene
    df.umap.cur$is.mod <- df.umap.cur$var %in% module.merge.genes[[mod.name]]
    n.mod.gene <- sum(df.umap.cur$is.mod)
    top.n <- 30
    if (top.n > n.mod.gene) top.n <- n.mod.gene

    # cell.emb
    # feat.emb
    if (module.type == "ica"){
      feat.emb.cur <- data.frame(gene = rownames(feat.emb),
                                 loading = feat.emb[,gsub("IC", "IC_", mod.name)])
    } else if (module.type == "nmf"){
      feat.emb.cur <- data.frame(gene = rownames(feat.emb),
                                 loading = feat.emb[,mod.name])
    } else if (module.type == "ssn"){

      if (!is.null(df.con)){
        feat.emb.cur <-       df.con[ ,c("genes", "wi")]
        colnames(feat.emb.cur) <- c("gene", "loading")
      } else {
        feat.emb.cur <- NULL
      }

    }

    if (!is.null(feat.emb.cur)){
      df.umap.cur <- merge(df.umap.cur, feat.emb.cur, by = "gene", all.x = T)
      min.size <- 0.1 * min(abs(df.umap.cur$loading[df.umap.cur$is.mod]), na.rm = T)
      df.umap.cur$loading[!df.umap.cur$is.mod] <- min.size
    } else {
      min.size <- 0.01
      df.umap.cur$loading <- 1
      df.umap.cur$loading[!df.umap.cur$is.mod] <- 0.01
    }

    df.umap.cur$loading <- abs(df.umap.cur$loading)
    df.umap.top <- df.umap.cur %>% dplyr::top_n(top.n, loading)
    df.umap.top <- df.umap.top %>% dplyr::filter(loading != min.size)
    if (nrow(df.umap.top) > top.n) df.umap.top <- df.umap.top[1:top.n, ][1:top.n, ]

    if (is.null(connectivity_plot)){
      plt.net <-  ggplot(data = df.umap.cur) +
        geom_point( pch = 21, aes(x = x, y = y, fill = is.mod, color = is.mod, size = loading)) + theme_miko() +
        theme_void() + theme(legend.position = "none") +
        labs(x = "geneUMAP 1", y = "geneUMAP 2",
             title = "Transcription network", subtitle = paste0(mod.name , " (", ngenes, " genes)")  ) +
        scale_fill_manual(values = c("TRUE" = as.vector(col.pal[mod.name]), "FALSE" = "grey98")) +
        scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "grey80")) +
        scale_size(range = c(0.5,5)) +
        ggrepel::geom_text_repel(data = df.umap.top, aes(x = x,  y =y, label = gene),
                                 inherit.aes = F, max.overlaps = Inf)
    } else {
      df.umap.cur <- df.umap.cur %>% dplyr::filter(is.mod)
      plt.net <-  connectivity_plot +
        geom_point(data = df.umap.cur, pch = 21, aes(x = x, y = y, fill = is.mod, color = is.mod, size = loading)) + theme_miko() +
        theme_void() + theme(legend.position = "none") +
        labs(x = "geneUMAP 1", y = "geneUMAP 2",
             title = "Transcription network", subtitle = paste0(mod.name , " (", ngenes, " genes)")  ) +
        scale_fill_manual(values = c("TRUE" = as.vector(col.pal[mod.name]), "FALSE" = "grey98")) +
        scale_color_manual(values = c("TRUE" = "grey", "FALSE" = "grey80")) +
        scale_size(range = c(0.5,5)) +
        ggrepel::geom_text_repel(data = df.umap.top, aes(x = x,  y =y, label = gene), inherit.aes = F, max.overlaps = Inf)
    }


    plt.bp <- bp.sum$plots[[mod.name]] + labs(subtitle = "BP")
    plt.mf <- mf.sum$plots[[mod.name]] + labs(subtitle = "MF")
    plt.cc <-cc.sum$plots[[mod.name]] + labs(subtitle = "CC")

    # plt.sum <- cowplot::plot_grid(plt.cell, plt.net, plt.bp, plt.mf, plt.cc, ncol = 5)
    plt.sum <- list(cell.umap = plt.cell,
                    gene.umap = plt.net,
                    bp.enrich = plt.bp,
                    mf.enrich = plt.mf,
                    cc.enrich = plt.cc)
    return(plt.sum)
  })

  names(plt.list) <- all.names

  return(
    list(
      plt.summary = plt.list,
      plt.heatmap = plt.expr,
      data.heatmap = cluster.activity,
      bp.enrichment = bp.res,
      mf.enrichment = mf.res,
      cc.enrichment = cc.res,
      module.genes = module.merge.genes
    )
  )
}
