
#' Run AUCell classification
#'
#' Run AUCell classification to identify active genesets in scRNAseq sample. Wrapper for AUCell package.
#'
#' @param object Seurat object
#' @param genelist Named list of genesets.
#' @param assay Assay used for expression matrix.
#' @param n.workers Number of workers for parallelized implementation. Default is 1.
#' @param n.repeats Number of fitting repeats (best fit taken). Default is 5.
#' @param n.iterations Number of fitting iterations per a repeat. Default is 1000.
#' @param posterior.p Posterior distribution membership threshold. Default is 0.9 (e.g., p > 0.9 is member).
#' @param mixture.analysis logical to perform mixture analysis. Default = T (computationally intensive).
#' @param size UMAP point size.
#' @param ... additional parameters passed to geom_point(...)
#' @name runAUC
#' @seealso \code{\link{AUCell_calcAUC }}
#' @author Nicholas Mikolajewicz
#' @return list of results along with ggplot handles visualizing class predictions overlaid on UMAP.
#' @examples
#'
#' # get genesets
#' verhaak.df <- geneSets[["Verhaak_CancerCell_2010"]]
#' verhaak.list <- wideDF2namedList(verhaak.df)
#'
#' gsc.df <- geneSets[["Richards_NatureCancer_2021_sc"]]
#' gsc.list <- wideDF2namedList(gsc.df)
#'
#' neftel.df <- geneSets[["GBM_Hs_Neftel2019"]]
#' neftel.list <- wideDF2namedList(neftel.df)
#'
#' verhaak.list <- lapply(verhaak.list, toupper)
#' gsc.list <- lapply(gsc.list, toupper)
#' neftel.list <- lapply(neftel.list, toupper)
#'
#' # classify cells based on provided genesets
#' v.auc <- runAUC(object = so.query, genelist = verhaak.list, n.workers = 12)
#' v.auc$plot.auc
#' v.auc$plot.nm
#' v.auc$plot.max.score
#'
#' g.auc <- runAUC(object = so.query, genelist = gsc.list, n.workers = 12)
#' g.auc$plot.auc
#' g.auc$plot.nm
#' g.auc$plot.max.score
#'
#' n.auc <- runAUC(object = so.query, genelist = neftel.list, n.workers = 12)
#' n.auc$plot.auc
#' n.auc$plot.nm
#' n.auc$plot.max.score
#'
runAUC <- function(object, genelist, assay = DefaultAssay(object), n.workers = 1, n.repeats = 5, n.iterations = 1000, posterior.p = 0.9, mixture.analysis = T, size = autoPointSize(ncol(object)), ...){

  suppressMessages({
    suppressWarnings({



      require(AUCell)
      so.e.mat <- object@assays[[assay]]@data

      n.auc.cores <- n.workers
      cells_rankings <- AUCell_buildRankings(so.e.mat, plotStats = F, nCores = n.auc.cores, verbose = F)
      cells_AUC <- AUCell_calcAUC(geneSets = genelist, rankings = cells_rankings,
                                  aucMaxRank=nrow(cells_rankings)*0.05,verbose = F, nCores = n.auc.cores)

      auc.val <- getAUC(cells_AUC)

      class.prediction <- list()
      for (i in 1:nrow(auc.val)){

        which.ind <- i
        which.class <- rownames(auc.val)[i]


        best.mix <- NULL
        mixmdl <- NULL

        if (mixture.analysis){

        for (j in 1:n.repeats){
          mixmdl = tryCatch(mixtools::normalmixEM(auc.val[which.ind, ],
                                                  k = 2, maxit = n.iterations, maxrestarts = 10, verb = FALSE), error = function(e) {
                                                    # print("EM algorithm did not converge")
                                                    mixmdl = NULL
                                                  })

          if (!is.null(mixmdl)){
            if (!is.null(best.mix)){
              if (mixmdl$loglik < best.mix$loglik){
                best.mix <- mixmdl
              }
            } else {
              best.mix <- mixmdl
            }
          }

        }




        cur.auc.vals <- auc.val[which.ind, ]
        if (length(cur.auc.vals) > 5000){
          cur.auc.vals0 <- sample(cur.auc.vals,5000, replace = F)
          shap <- shapiro.test(cur.auc.vals0)
        } else {
          shap <- shapiro.test(cur.auc.vals)
        }

        if (shap[["p.value"]]<0.05){
          is.norm <- F
        } else {
          is.norm <- T
        }

        loglik_2 <- mixmdl[["loglik"]]
        loglik_1 <- sum(dnorm(auc.val[which.ind, ], mean(auc.val[which.ind, ]), sd(auc.val[which.ind, ]), log = TRUE))

        } else {
          is.norm <- T
          loglik_2 <- NULL
          loglik_1 <- NULL
        }




        single.comp <- F
        try({
          if (loglik_1 < loglik_2){
            single.comp <- T
          } else {
            single.comp <- F
          }
        }, silent = T)


        if (single.comp & is.norm){
          class.prediction[[which.class]] <- character()
        } else {
          class.prediction[[which.class]] <- colnames(auc.val)[best.mix[["posterior"]][ ,which.max(best.mix[["mu"]])] > posterior.p]
        }

      }

      cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=F, nCores=n.auc.cores, assign=TRUE,  verbose = TRUE)
      df.auc.umap <- data.frame(
        cells = colnames(object),
        x = object@reductions[["umap"]]@cell.embeddings[ ,1],
        y = object@reductions[["umap"]]@cell.embeddings[ ,2]
      )

      # a.list <- list()
      df.auc.umap$class.auc <- NA
      df.auc.umap$class.nm <- NA

      score.mat.nm <- matrix(data = 0, nrow = nrow(df.auc.umap), ncol = length(class.prediction))
      class.mat.nm <- matrix(data = 0, nrow = nrow(df.auc.umap), ncol = length(class.prediction))
      score.mat.auc <- matrix(data = 0, nrow = nrow(df.auc.umap), ncol = length(class.prediction))
      class.mat.auc <- matrix(data = 0, nrow = nrow(df.auc.umap), ncol = length(class.prediction))

      for (i in 1:length(class.prediction)){
        cur.set <- names(class.prediction)[i]
        cur.assignment <- cells_assignment[[i]][["assignment"]]
        cur.assignment2 <- class.prediction[[i]]

        score.mat.auc[ , i] <- auc.val[cur.set, ]
        class.mat.auc[ , i] <- (df.auc.umap$cells %in% cur.assignment)*1

        score.mat.nm[ , i] <- auc.val[cur.set, ]
        class.mat.nm[ , i] <- (df.auc.umap$cells %in% cur.assignment2)*1
        df.auc.umap[ ,cur.set ] <- auc.val[cur.set, ]
      }

      valid.class.scores.nm <- score.mat.nm*class.mat.nm
      valid.class.scores.auc <- score.mat.auc*class.mat.auc

      max.class.nm <- apply(data.frame((valid.class.scores.nm)), 1, which.max)
      max.class.auc <- apply(data.frame((valid.class.scores.auc)), 1, which.max)

      df.auc.umap$class.nm <- rownames(auc.val)[max.class.nm]
      df.auc.umap$class.nm[apply(data.frame((valid.class.scores.nm)), 1, max) == 0] <- "unresolved"
      df.auc.umap$class.auc <- rownames(auc.val)[max.class.auc]
      df.auc.umap$class.auc[apply(data.frame((valid.class.scores.auc)), 1, max) == 0] <- "unresolved"

      max.class <- apply(data.frame(t(auc.val)), 1, which.max)
      df.auc.umap$class.max.score <- rownames(auc.val)[max.class]

      if ((length(genelist) > 10)){
        colpal <- scales::hue_pal()(length(genelist))
      } else if ((length(genelist) <= 10) & (length(genelist) > 1)){
        colpal <- ggthemes::ptol_pal()(length(genelist))
      } else {
        colpal <- "tomato"
      }
      colpal <- c(colpal, "grey")
      names(colpal) <- c(names(genelist), "unresolved")

      df.auc.umap$class.auc <- factor(df.auc.umap$class.auc, levels = names(colpal))
      df.auc.umap$class.nm <- factor(df.auc.umap$class.nm, levels = names(colpal))
      df.auc.umap$class.max.score <- factor(df.auc.umap$class.max.score, levels = names(colpal))

      plt.umap.list <- list()
      for (i in 1:nrow(auc.val)){

        which.set <- rownames(auc.val)[i]
        df.auc.umap2 <- df.auc.umap
        df.auc.umap2$auc <- auc.val[i, ]

        plt.umap.list[[which.set]] <- df.auc.umap2 %>%
          ggplot(aes(x = x, y = y, color = auc)) +
          geom_point(size = size, ...) +
          labs(x = "UMAP 1", y = "UMAP 2", caption = "AUCell-based scoring", title = which.set, subtitle = "AUC scores") +
          theme_miko(center.title = T, legend = T) +
          scale_color_gradient2(high = "red")
      }

      plt.auc.umap <- df.auc.umap %>%
        dplyr::arrange(-as.numeric(class.auc)) %>%
        ggplot(aes(x = x, y = y, color = class.auc)) +
        geom_point(size = size, ...) +
        labs(x = "UMAP 1", y = "UMAP 2", caption = "Original AUCell-based classification", color = "Class") +
        theme_miko(center.title = T, legend = T) +
        scale_color_manual(values = colpal) +
        guides(fill = F, color = guide_legend(override.aes = list(size = 4)))

      plt.nm.umap <- df.auc.umap %>%
        dplyr::arrange(-as.numeric(class.nm)) %>%
        ggplot(aes(x = x, y = y, color = class.nm)) +
        geom_point(size = size, ...) +
        labs(x = "UMAP 1", y = "UMAP 2", caption = "NM-modified AUCell-based classification", color = "Class") +
        theme_miko(center.title = T, legend = T) +
        scale_color_manual(values = colpal) +
        guides(fill = F, color = guide_legend(override.aes = list(size = 4)))


      plt.max.umap <- df.auc.umap %>%
        ggplot(aes(x = x, y = y, color = class.max.score)) +
        geom_point(size = size, ...) +
        labs(x = "UMAP 1", y = "UMAP 2", caption = "Classification based on max AUCell score", color = "Class") +
        theme_miko(center.title = T, legend = T) +
        scale_color_manual(values = colpal) +
        guides(fill = F, color = guide_legend(override.aes = list(size = 4)))

    })
  })

  closeAllConnections()
  rm(object)
  invisible({gc()})

  return(
    list(
      data = df.auc.umap,
      plot.auc = plt.auc.umap,  # auc-based approach (+rejection)
      plot.nm = plt.nm.umap, # loglikelihood and normality test-based classification (+rejection)
      plot.max.score = plt.max.umap, # max score across genesets (-rejection)
      plot.list = plt.umap.list,
      class.prediction = class.prediction
    )
  )

}


#' Run Modular Scoring
#'
#' Wrapper for Seurat::AddModuleScore(). Calculates module scores for feature expression programs in single cells.
#'
#' @param object Seurat object
#' @param genelist Named list of genesets.
#' @param assay Assay used for expression matrix.
#' @param score.key Expression program prefix. default is "MS".
#' @param size UMAP point size.
#' @param ncore Number of workers for parallelized implementation. Default is 10.
#' @param raster Convert points to raster format, default is FALSE.
#' @param rescale rescale values from 0 to 1. Default is FALSE.
#' @param verbose Print progress. Default is TRUE.
#' @param winsorize.quantiles Rescale values to lie between lower and upper bound quanitle. Default  = c(0,1).
#' @param return.plots Logical to compute and return plots in results list. Default is TRUE.
#' @param search Search for symbol synonyms for features in features that don't match features in object? Searches the HGNC's gene names database; see UpdateSymbolList for more details. Default is FALSE.
#' @param ... additional parameters passed to geom_point(...)
#' @name runMS
#' @seealso \code{\link{AddModuleScore}}
#' @author Nicholas Mikolajewicz
#' @return list of results along with ggplot handles visualizing class predictions overlaid on UMAP.
#' @examples
#'
#' # get genesets
#' verhaak.df <- geneSets[["Verhaak_CancerCell_2010"]]
#' verhaak.list <- wideDF2namedList(verhaak.df)
#'
#' gsc.df <- geneSets[["Richards_NatureCancer_2021_sc"]]
#' gsc.list <- wideDF2namedList(gsc.df)
#'
#' neftel.df <- geneSets[["GBM_Hs_Neftel2019"]]
#' neftel.list <- wideDF2namedList(neftel.df)
#'
#' verhaak.list <- lapply(verhaak.list, toupper)
#' gsc.list <- lapply(gsc.list, toupper)
#' neftel.list <- lapply(neftel.list, toupper)
#'
#' # classify cells based on provided genesets
#' v.auc <- runMS(object = so.query, genelist = verhaak.list)
#' v.auc$plot.max.score
#'
#' g.auc <- runMS(object = so.query, genelist = gsc.list)
#' g.auc$plot.max.score
#'
#' n.auc <- runMS(object = so.query, genelist = neftel.list)
#' n.auc$plot.max.score
#'
runMS <- function(object, genelist, assay = DefaultAssay(object), score.key = "MS", size = autoPointSize(ncol(object)), ncore = 10, raster = F, rescale = F, verbose = T, winsorize.quantiles = c(0,1), return.plots = T, search = F, ...){

  require(scales);

  opt.bsize <- optimalBinSize(object, verbose = verbose)

  if (class(genelist) == "character"){
    genelist <- list(geneset = genelist)
  }

  # start cluster
  if (ncore > 1){
    if (ncore > detectCores()){
      ncore <- detectCores()
    } else {
      ncore <- ncore
    }
  }


  miko_message("Scoring gene modules...", verbose = verbose)
  if (ncore > 1){
    object <-  AddSModuleScore(
      object = object,
      features = genelist,
      pool = NULL,
      nbin = opt.bsize,
      ctrl = 100,
      ncore = ncore,
      k = FALSE,
      assay = assay,
      name = score.key,
      seed = 1,
      search = search
    )
  } else{
    object <-   Seurat::AddModuleScore(
      object = object,
      features = genelist,
      pool = NULL,
      nbin = opt.bsize,
      ctrl = 100,
      k = FALSE,
      assay = assay,
      name = score.key,
      seed = 1,
      search = search
    )
  }




  df.ms <- object@meta.data[ ,grepl(score.key, colnames(object@meta.data) )]
  if (is.numeric(df.ms)){
    df.ms <- as.data.frame(df.ms)
    rownames(df.ms) <- colnames(object)
  }
  colnames(df.ms) <- names(genelist)

  miko_message("Consolidating results...", verbose = verbose)
  which.max.score <- apply(df.ms, 1, which.max)
  class.prediction.max <- colnames(df.ms)[which.max.score]

  if (return.plots){
  df.auc.umap <- data.frame(
    cells = colnames(object),
    x = object@reductions[["umap"]]@cell.embeddings[ ,1],
    y = object@reductions[["umap"]]@cell.embeddings[ ,2]
  )

  df.auc.umap <- bind_cols(df.auc.umap, df.ms)
  df.auc.umap$class.ms <- class.prediction.max

  if ((length(genelist) > 10)){
    colpal <- scales::hue_pal()(length(genelist))
  } else if ((length(genelist) <= 10) & (length(genelist) > 1)){
    colpal <- ggthemes::ptol_pal()(length(genelist))
  } else {
    colpal <- "tomato"
  }
  colpal <- c(colpal, "grey")
  names(colpal) <- c(names(genelist), "unresolved")

  df.auc.umap$class.ms <- factor(df.auc.umap$class.ms, levels = names(colpal))
  miko_message("Generating plots...", verbose = verbose)

  if (raster){
    geom_point_fun <- scattermore::geom_scattermore
    size <- 1
  } else {
    geom_point_fun <- geom_point
  }

  plt.umap.list <- list()
  for (i in 1:ncol(df.ms)){

    which.set <- colnames(df.ms)[i]

    if (rescale){
      df.auc.umap[ ,which.set] <- scMiko::rescaleValues(df.auc.umap[ ,which.set])
    }

    if ((winsorize.quantiles[1] != 0) & (winsorize.quantiles[2] != 1)){
      if (winsorize.quantiles[1] < winsorize.quantiles[2]){
        if ((winsorize.quantiles[1] > 0) & (winsorize.quantiles[1] < 1)){
          q1 <- quantile(df.auc.umap[ ,which.set], winsorize.quantiles[1])
        } else {
          q1 <- min(df.auc.umap[ ,which.set], na.rm = T)
        }
        if ((winsorize.quantiles[2] > 0) & (winsorize.quantiles[2] < 1)){
          q2 <- quantile(df.auc.umap[ ,which.set], winsorize.quantiles[2])
        } else {
          q2 <- max(df.auc.umap[ ,which.set], na.rm = T)
        }

        df.auc.umap[df.auc.umap[ ,which.set] < q1 ,which.set] <- q1
        df.auc.umap[df.auc.umap[ ,which.set] > q2 ,which.set] <- q2
      }
    }


    df.auc.umap2 <- df.auc.umap
    df.auc.umap2$auc <- df.auc.umap2[ ,which.set]

    plt.umap.list[[which.set]] <- df.auc.umap2 %>%
      ggplot(aes(x = x, y = y, color = auc)) +
      geom_point_fun(size = size, ...) +
      labs(x = "UMAP 1", y = "UMAP 2", color = "Score",  title = which.set, subtitle = "Modular scores") +
      theme_miko(center.title = T, legend = T) +
      scale_color_gradient2(low = muted("blue"),
                            mid = "white",
                            high = muted("red"))
  }

  } else {
    plt.umap.list <- list()
    df.auc.umap <- df.ms
    df.auc.umap$class.ms <- class.prediction.max
    df.auc.umap$class.ms <- factor(df.auc.umap$class.ms, levels = c(names(genelist), "unresolved"))
  }


  if (return.plots){
    plt.max.umap <- df.auc.umap %>%
      dplyr::arrange(-as.numeric(class.ms)) %>%
      ggplot(aes(x = x, y = y, color = class.ms)) +
      # geom_point() +
      geom_point_fun(size = size, ...) +
      labs(x = "UMAP 1", y = "UMAP 2", caption = "Classification based on max modular score", color = "Class") +
      theme_miko(center.title = T, legend = T) +
      scale_color_manual(values = colpal) +
      guides(fill = F, color = guide_legend(override.aes = list(size = 4)))
  } else {
    plt.max.umap <- NULL
  }



# try({
#   # closeAllConnections()
#   # rm(object)
#   invisible({gc()})
# }, silent = T)


  return(
    list(
      data = df.auc.umap,
      plot.max.score = plt.max.umap, # max score across genesets (-rejection)
      plot.list = plt.umap.list,
      class.prediction = class.prediction.max
    )
  )

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
miko_integrate <- function(object, split.by = "Barcode", min.cell = 50, k.anchor = 20, k.weight = 35, nfeatures = 3000, split.prenorm = F, assay = "SCT", variable.features.n = 3000, verbose = T, use.existing.sct = F, conserve.memory = F, vars.to.regress = "percent.mt"){

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

    # renormalize
    if (verbose.) miko_message("Running SCTransform...", verbose = verbose.)
    object.list <- pbapply::pblapply(X = object.list, FUN = SCTransform, method = "glmGamPoi", verbose = verbose., assay = assay, conserve.memory = conserve.memory,
                                     vars.to.regress = vars.to.regress, variable.features.rv.th = 1.3, variable.features.n = variable.features.n)

  } else {
    if (verbose.) miko_message("Running SCTransform...", verbose = verbose.)
    object <- SCTransform(object, method = "glmGamPoi", verbose = verbose., assay = assay,
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

  if (verbose.) miko_message("Complete!")

  return(object2)

}


#' Get differentially expressed genes
#'
#' Differential expression analysis wrapper for presto::wilcoxauc(). Option to return differentially-expressed gene list (return.list = T) or statistics only (return.list = F)
#'
#' @param object Seurat object
#' @param assay assay. Default "SCT".
#' @param data data slot. Default "data".
#' @param group_by Name of grouping variable (must be present in object's meta.data)
#' @param auc.thresh AUC threshold. Default = 0.6.
#' @param fdr.thresh FDR threshold. Default = 0.01.
#' @param logFC.thresh logFC threshold. Default = NA.
#' @param pct.dif.thresh Difference in expression percentage. Default = NA.
#' @param pct.in.thresh Expression percentages exceeding this threshold are retained. Default is NA.
#' @param pct.out.thresh Expression percentages below this threshold are retained Default is NA.
#' @param return.list If TRUE, return list of differentially expressed genes. If FALSE, returns table with statistics from differential expression analysis.
#' @param return.all If TRUE, all thresholding filters are ignored, and all results are returned.
#' @param sig.figs If specified and return.list = F, rounds statistics to specified significant figure (recommended: 3). Default is NA.
#' @param verbose Print progress. Default is TRUE.
#' @name getDEG
#' @seealso \code{\link{wilcoxauc}}
#' @author Nicholas Mikolajewicz
#' @return data.frame or list
#' @examples
#'
getDEG <- function(object, assay = DefaultAssay(object), data = "data",
                   group_by = "seurat_clusters", auc.thresh = 0.6, fdr.thresh = 0.01, logFC.thresh = NA, pct.dif.thresh = NA, pct.in.thresh = NA, pct.out.thresh= NA, return.list = T, return.all = F, sig.figs = NA, verbose = T){

  require(presto)
  stopifnot("Seurat" %in% class(object) )

  miko_message("Running differential expression analysis...", verbose = verbose)
  deg.dat <- presto::wilcoxauc(X = object , assay = data, seurat_assay = assay, group_by = group_by)
  deg.dat$pct.dif <- deg.dat$pct_in - deg.dat$pct_out

  deg.dat$sensitivity <- deg.dat$pct_in/((100-deg.dat$pct_in) + (deg.dat$pct_in))
  deg.dat$specificity <- (100-deg.dat$pct_out)/((100-deg.dat$pct_out) + (deg.dat$pct_out))
  deg.dat$PPV <-  deg.dat$pct_in/( deg.dat$pct_in + (deg.dat$pct_out))
  deg.dat$NPV <- (100-deg.dat$pct_out)/((100-deg.dat$pct_out) + (100-deg.dat$pct_in))


  if (!return.all){
    miko_message("Applying thresholds...", verbose = verbose)


    if (!is.na(auc.thresh)){
      deg.dat <- deg.dat %>% dplyr::filter(auc > auc.thresh)
    }
    if (!is.na(fdr.thresh)){
      deg.dat <- deg.dat %>% dplyr::filter(padj < fdr.thresh)
    }

    if (!is.na(logFC.thresh)){
      deg.dat <- deg.dat %>% dplyr::filter(logFC > logFC.thresh)
    }

    if (!is.na(pct.dif.thresh)){
      deg.dat <- deg.dat %>% dplyr::filter(pct.dif > pct.dif.thresh)
    }

    if (!is.na(pct.in.thresh)){
      deg.dat <- deg.dat %>% dplyr::filter(pct_in > pct.in.thresh)
    }

    if (!is.na(pct.out.thresh)){
      deg.dat <- deg.dat %>% dplyr::filter(pct_out < pct.out.thresh)
    }


  }

  if (return.list){
    # miko_message("Preparing list output...", verbose = verbose)
    u.group <- unique(deg.dat$group)

    deg.list <- list()
    for (i in 1:length(u.group)){
      deg.list[[u.group[i]]] <- deg.dat$feature[deg.dat$group %in% u.group[i]]
    }
    # miko_message("Complete!", verbose = verbose)
    return(deg.list)
  } else {
    # miko_message("Complete!", verbose = verbose)

    if (!is.na(sig.figs)){
      try({

        deg.dat[ ,c("avgExpr", "logFC", "statistic", "auc", "pval", "padj", "pct_in", "pct_out", "pct.dif", "sensitivity", "specificity", "PPV", "NPV")] <- signif(deg.dat[ ,c("avgExpr", "logFC", "statistic", "auc", "pval", "padj", "pct_in", "pct_out", "pct.dif", "sensitivity", "specificity", "PPV", "NPV")], sig.figs)
      }, silent = T)
    }
    return(deg.dat)
  }

}




#' Run Robust Prinicipal Component Analysis
#'
#' Run a robust PCA (rPCA) dimensionality reduction on single-cell seurat object.
#'
#' @param object Seurat object
#' @param assay Name of Assay rPCA is being run on
#' @param features Features to compute PCA on. If features=NULL, PCA will be run using scaled features for the Assay. Note that the features must be present in the scaled data. Any requested features that are not scaled or have 0 variance will be dropped, and the PCA will be run using the remaining features.
#' @param npcs Total Number of PCs to compute and store (50 by default)
#' @param maxpcs Max Number of PCs to compute and store (50 by default)
#' @param reductions.key dimensional reduction key, specifies the string before the number for the dimension names. RPC by default
#' @param reduction.name dimensional reduction name, rpca by default
#' @param seed.use Set a random seed. By default, sets the seed to 42. Setting NULL will not set a seed.
#' @param verbose Print progress. Default is TRUE.
#' @param method Robust PCA method. default is "hubert".
#' @param maxdir maximal number of random directions to use for computing the outlyingness of the data points. Default is maxdir=100.
#' @param signflip a logical value indicating wheather to try to solve the sign indeterminancy of the loadings - ad hoc approach setting the maximum element in a singular vector to be positive. Default is signflip = TRUE
#' @param ... additional parameters passed to rPCA methods.
#' @name runRPCA
#' @seealso \code{\link{PcaHubert}}
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#' @examples
#'
runRPCA <- function(object, assay = NULL, features = NULL, npcs = 50, maxpcs = 50,  reduction.key = "RPC_", reduction.name = "rpca", seed.use = 42, verbose = T,
                    method = c( "hubert", "robpca", "fasthcs", "pcal"), maxdir = 100,signflip = T, ...){

  # computation time (hubert):
  # 1000 cells ~ 2 min
  # 17000 cells ~ 10 min

  # computation time (pcal):
  # 1000 cells ~ 1 min
  # 17000 cells ~ 8 min

  # recommended algorithms for scRNAseq: hubert, pcal
  set.seed(seed.use)

  if ("Seurat" %in% class(object)){
    if (is.null(assay)){
      assay <- DefaultAssay(object)
    } else {
      DefaultAssay(object) <- assay
    }
    emat <- object@assays[[assay]]@scale.data
  } else if ("matrix" %in% class(object)){
    emat <- object
    assay <- "temp"
  }

  # emat <- object@assays[[assay]]@scale.data
  if (!is.null(features)){
    emat <- emat[rownames(emat) %in% features, ]
  }

  miko_message("Running robust PCA...", verbose = verbose)

  if (method == "hubert"){
    require(rrcov)
    pca.red <- rrcov::PcaHubert (x = emat, k = npcs, kmax = maxpcs, trace = verbose, maxdir = maxdir, signflip = T, ...)
    res.list <- list(
      embeddings = pca.red@loadings,
      loadings = pca.red@scores,
      stdev = pca.red@eigenvalues
    )

  } else if (method == "robpca"){
    require(rospca)
    pca.red <- rospca::robpca (x = emat, k = npcs, kmax = maxpcs,  ndir = maxdir, ...)
    res.list <- list(
      embeddings = pca.red[["loadings"]],
      loadings = pca.red[["scores"]],
      stdev = pca.red[["eigenvalues"]]
    )

  } else if (method == "fasthcs"){
    require(FastHCS)
    # only works for q <= 25
    if (maxpcs > 25) maxpcs <- 25
    pca.red <- FastHCS::FastHCS (x = emat, q = maxpcs, seed = seed.use)
    res.list <- list(
      embeddings = pca.red[["loadings"]],
      loadings = pca.red[["scores"]],
      stdev = pca.red[["eigenvalues"]]
    )
  } else if (method == "pcal"){
    require(rrcov)
    pca.red <- rrcov::PcaLocantore (x = emat, k = npcs, kmax = maxpcs, trace = verbose, signflip = T)
    res.list <- list(
      embeddings = pca.red@loadings,
      loadings = pca.red@scores,
      stdev = pca.red@eigenvalues
    )
  }

  if ("Seurat" %in% class(object)){
    object[[reduction.name]] <- CreateDimReducObject(embeddings = res.list$embeddings,
                                                     loadings = res.list$loadings,
                                                     stdev = res.list$stdev,
                                                     key = reduction.key,
                                                     assay =assay, ...)

  } else if ("matrix" %in% class(object)){
    object <-  CreateDimReducObject(embeddings = res.list$embeddings,
                                    loadings = res.list$loadings,
                                    stdev = res.list$stdev,
                                    key = reduction.key,
                                    assay =assay, ...)
  }

  miko_message("Complete!", verbose = verbose)
  return(object)

}


#' Evaluate signature coherence.
#'
#' Evaluates signature coherence by computing signature score (using runMS) and then calculating correlations between signature components (i.e., gene expression) and signature score. Genes for which correlations do not exceed coherence threshold are dropped and remaining genes are used to calculate new coherent signature.
#'
#' @param object Seurat object
#' @param ms.result Output from runMS() for provided genelist. If not specified, signatures are computed using runMS.
#' @param genelist Named list of genes used to compute gene signature.
#' @param slot Seurat slot used for gene expression correlations. Default is 'data.'
#' @param assay Seurat assay. Default is DefaultAssay(object).
#' @param coherence.threshold Coherence threshold. Default is 0.1.
#' @param show.grid Logical to show grid on coherence graph.
#' @param coherent.ms Logical to recompute coherent signature following coherence thresholding.
#' @param cor.method correlation metod used to determine signature coherence. Options: spearman, pearson. Default is spearman.
#' @param ... additional parameters passed to scMiko::runMS().
#' @name signatureCoherence
#' @seealso \code{\link{runMS}}
#' @author Nicholas Mikolajewicz
#' @return list of results
#' @examples
#'
signatureCoherence <- function(object = NULL, ms.result = NULL, genelist,
                               slot = "data", assay = DefaultAssay(object), coherence.threshold = 0.1, show.grid = T,
                               coherent.ms = T, cor.method = "spearman", ...){


  if (class(genelist) == "character"){
    genelist <- list(geneset = genelist)
  }


  if (is.null(ms.result)){
    ms.result <- runMS(object = object, genelist = genelist, return.plots = F, ...)
  }

  match.sets <- colnames(ms.result[["data"]])[colnames(ms.result[["data"]]) %in% names(genelist)]


  if (assay == "SCT" & slot == "scale"){
    object <-  GetResidual(object = object, features = unique(unlist(genelist)))
  }
  # , as.dense = T
  expr.mat <- getExpressionMatrix(so = object, only.variable = F, which.data = slot, which.assay = assay)
  if (sum(dim(expr.mat)) == 0) stop("Expression matrix is empty.")
  genelist.coherent <- list()
  cor.list <- list()
  plt.list <- list()

  miko_message("Computed signature score correlations...")
  for (i in 1:length(match.sets)){

    sig.score <- (as.matrix(ms.result$data[ ,match.sets[i]]))
    # expr.score <- t(expr.mat[rownames(expr.mat) %in% genelist[[match.sets[i]]], ])
    expr.score <-  tryCatch({
      expr.score <- t(expr.mat[rownames(expr.mat) %in% genelist[[match.sets[i]]], ])
    }, error = function(e){
      expr.mat <- as.matrix(expr.mat)
      expr.score <- t(expr.mat[rownames(expr.mat) %in% genelist[[match.sets[i]]], ])
      return(expr.score)
    })

    if (sum(dim(expr.score)) == 0) next

    if (cor.method == "pearson"){
      if (all(class(expr.score) == class(sig.score))){
        cor.score <- cor(expr.score, sig.score)
      } else {
        cor.score <- qlcMatrix::corSparse(expr.score, sig.score)
      }
    } else if (cor.method == "spearman"){
      cor.score <- cor(as.matrix(expr.score), sig.score)
    }

    cor.score[is.na(cor.score)] <- 0

    # a <- qlcMatrix::corSparse(expr.score, expr.score)
    # rownames(a) <- colnames(a) <- colnames(expr.score)

    rownames(cor.score) <- colnames(expr.score)
    colnames(cor.score) <- match.sets[i]

    df.cor.score <- as.data.frame(cor.score)
    colnames(df.cor.score) <- "score"
    df.cor.score$gene <- rownames(df.cor.score)


    raw.score <- signif(mean(df.cor.score$score), 3)
    coh.score <- signif(mean(df.cor.score$score[df.cor.score$score > coherence.threshold]), 3)

    raw.genes <- length(df.cor.score$score)
    coherent.genes <- sum(df.cor.score$score > coherence.threshold)

    df.cor.score$is.coherent <- df.cor.score$score > coherence.threshold
    df.cor.score$raw.score <- raw.score
    df.cor.score$coh.score <- coh.score

    plt.cor.score <- df.cor.score %>%
      ggplot(aes(x = score, y = reorder(gene, score), color = is.coherent)) +
      geom_point() +
      geom_segment(aes(x = 0, xend = score, y = gene, yend = gene)) +
      geom_vline(xintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = coherence.threshold, linetype = "dashed", color = "tomato") +
      theme_miko(legend = T, grid = show.grid) +
      scale_color_manual(values = c("TRUE" = "black", "FALSE" = "grey")) +
      labs(x = "Correlation with Signature Score (r)", y = "Signature Genes",
           caption = "black dashed line: r = 0\nred dashed line: coherence threshold",
           title = paste0(match.sets[i], " score coherence"),
           subtitle = paste0("Raw score = ", raw.score, " (", raw.genes, "/", raw.genes, " genes)\nCoherent score = ",
                             coh.score, " (", coherent.genes , "/", raw.genes, " genes)"))

    # print(plt.cor.score)
    plt.list[[match.sets[i]]] <- plt.cor.score
    cor.list[[match.sets[i]]] <- df.cor.score
    if (coherent.genes > 0){
      genelist.coherent[[match.sets[i]]] <- df.cor.score$gene[df.cor.score$is.coherent]
    }

  }


  if (coherent.ms){
    miko_message("Recomputing coherent signatures...")
    if (length(genelist.coherent) > 0){
      ms.coherent.result <-  tryCatch({
        runMS(object = object, genelist = genelist.coherent, ...)
      }, error = function(e){
        return(runMS(object = object, genelist = genelist.coherent,return.plots = F, ...))
      })

    } else {
      ms.coherent.result <- list()
    }

  } else {
    ms.coherent.result <- list()
  }


  return(list(
    genelist = genelist,
    genelist.coherent = genelist.coherent,
    coherence.plots = plt.list,
    ms.coherent.result = ms.coherent.result,
    coherence.data = cor.list
  ))

}


