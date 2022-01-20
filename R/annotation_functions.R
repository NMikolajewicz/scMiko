#' Null distribution for standardized module scores.
#'
#' Downsample number of cells in Seurat object by specified factor
#'
#' @param object Seurat Object
#' @param assay Name of assay to use
#' @param n.replicate Number of replicates
#' @param min.gs.size minimum gene set size
#' @param max.gs.size maximum gene set size
#' @param step.size step size; increment of sequence between `min.gs.size` and `max.gs.size`.
#' @param nworkers number of workers used for parallel implementation. Default is 1.
#' @param subsample.n Numeric [1,ncol(object)]. Number of cells to subsample. If specified, overides subsample.factor.
#' @param nbin Number of bins of aggregate expression levels for all analyzed features. See `AddModuleScore` for details.
#' @param verbose Print progress. Default is TRUE.
#' @name nullScore
#' @author Nicholas Mikolajewicz
#' @concept miko_score
#' @return list of results that are used as input into `mikoScore`.
#' @seealso \code{\link{AddSModuleScore}} for standardized module scoring, \code{\link{mikoScore}} for miko scoring, \code{\link{sigScore}} for derivation of p values for miko scores.
#' @examples
#' maxgs <- max(unlist(lapply(gene.list, length)))
#' if (maxgs > 200) maxgs <- 200
#' stepsize <- round(maxgs / 15)
#' if (stepsize < 1) stepsize <- 1
#' if (stepsize > 20) stepsize <- 10
#'
#' # Generating null distributions for genesets
#' bm.null <- nullScore(object = object, assay = "RNA, n.replicate = 10, nbin = 24,
#'                     min.gs.size = 2, max.gs.size = gs.size + 5, step.size = stepsize,
#'                     nworkers =20, verbose = T, subsample.n = 5000)
#'
#' # variance.raw.plot <- bm.null$variance.raw.plot
#' # variance.mean.plot <- bm.null$variance.mean.plot
#' # plt.null.model <- bm.null$mean.plot
#' # plt.corrected.plot <- bm.null$corrected.plot
#'
#' # Compute miko scores for genesets
#' object <- mikoScore(object = object, geneset = gene.list, nbin = 24,
#'                           nullscore = bm.null, assay = "RNA", nworkers = 20)
#'
#' # Get significant  miko scores for scored genesets
#' score.result <- sigScore(object = object,  geneset = gene.list, reduction = "umap")
#'
nullScore <- function(object, assay = DefaultAssay(object), n.replicate = 25,
                      min.gs.size = 2, max.gs.size = 100, step.size = 2, nworkers = 1, subsample.n = NULL, nbin = 24, verbose = T){


  require(foreach)
  if (!is.numeric(min.gs.size)) min.gs.size <- 2
  if (min.gs.size < 0) min.gs.size <- 2
  if (!is.numeric(max.gs.size)) max.gs.size <- 100
  if (max.gs.size < min.gs.size) max.gs.size <- min.gs.size + 10
  if (!is.numeric(step.size)) step.size <- 2
  if (!is.numeric(nbin)) nbin <- 24
  stopifnot("Seurat" %in% class(object))

  if (!is.null(subsample.n)){
    if (subsample.n < ncol(object)){
      set.seed(1023)
      object <- downsampleSeurat(object = object, subsample.n = subsample.n, verbose = verbose)
    }
  }

  # ncores <- nworkers
  nbin <- optimalBinSize(object = object, nbin = nbin)



  miko_message("Generating random gene sets...", verbose = verbose)
  base.seq <- unique(c(3:10, seq(10, 20, by = 2), seq(min.gs.size, max.gs.size)))
  n.rep <- n.replicate
  rsize <- rep(base.seq, n.rep)


  rgene.list <- list()
  for (j in 1:length(rsize)){
    rgene.list[[j]] <- sample(rownames(object), rsize[j])
  }
  names(rgene.list) <- paste0("gs", seq(1, length(rgene.list)))

  miko_message("Scoring gene modules...", verbose = verbose)

  if (nworkers > length(rgene.list)) nworkers <- length(rgene.list)

  object <-  AddSModuleScore(
    object = object,
    features = rgene.list,
    pool = NULL,
    nbin = nbin,
    ctrl = 100,
    k = FALSE,
    nworkers = nworkers,
    assay = assay,
    name = "MS",
    seed = 1,
    search = F
  )


  # df.ms <- object@meta.data[ ,grepl("MS", colnames(object@meta.data) )]
  df.ms <- object@meta.data[ ,paste0("MS", names(rgene.list) )]
  if (is.numeric(df.ms)){
    df.ms <- as.data.frame(df.ms)
    rownames(df.ms) <- colnames(object)
  }
  colnames(df.ms) <- names(rgene.list)

  scores <- as.matrix(df.ms)

  # }


  # compute null model
  miko_message("Fitting null models...", verbose = verbose)
  suppressMessages({
    suppressWarnings({
      df.var <- data.frame(size = rsize,
                           x.var = apply(scores, 2, var, na.rm = T),
                           x.mean = apply(scores, 2, mean, na.rm = T),
                           x.median =  apply(scores, 2, median, na.rm = T)
      )
    })
  })


  df.var <- df.var %>%
    dplyr::group_by(size) %>%
    dplyr::mutate(y.mean = median((x.var), na.rm = T),
                  y.sd = mad((x.var), na.rm = T) )

  df.var.sum <- df.var %>%
    dplyr::group_by(size) %>%
    dplyr::summarize(x.var = median(x.var),
                     x.mean = median(x.mean))

  miko_message("Generating summary plots...", verbose = verbose)
  variance.model <-  tryCatch({
    glm(formula = x.var ~ size, family = Gamma, data = df.var.sum )
  }, error = function(e){
    variance.model <- glm(formula = x.var ~ size, family = inverse.gaussian, data = df.var.sum )
    return(variance.model)
  })

  # df.var$sd.pred <-  sqrt(variance.model[["fitted.values"]])
  df.var$sd.pred <-sqrt(predict.glm(object = variance.model, newdata = df.var, type="response") )
  df.var.sum$sd.pred <-sqrt(predict.glm(object = variance.model, newdata = df.var.sum, type="response") )

  variance.raw.plot <- df.var %>%
    dplyr::arrange(size) %>%
    ggplot(aes(x = size, y = (x.var))) +
    # scattermore::geom_scattermore() +
    geom_point(pch = 21, color = "white", size = 2, fill = "black") +
    geom_line(aes(y = (sd.pred^2)), color = "red", size = 1) +
    labs(x = "Geneset Size", y = "Score Variance", title = "Null variance model") +
    theme_miko(grid = T)

  variance.mean.plot <- df.var.sum %>%
    dplyr::arrange(size) %>%
    ggplot(aes(x = size, y = (x.var))) +
    geom_point(pch = 21, color = "white", size = 2, fill = "black") +
    # scattermore::geom_scattermore() +
    geom_line(aes(y = (sd.pred^2)), color = "red", size = 1) +
    labs(x = "Geneset Size", y = "Score Variance", title = "Null variance model") +
    theme_miko(grid = T)


  model.offset <- 1 + abs(min(df.var.sum$x.mean))

  mean.model <-  tryCatch({
    library(splines)
    list(fit = glm(formula = (x.mean ) ~ bs(size, degree = 3), family = gaussian, data = df.var.sum ),
         model = "spline")
  }, error = function(e){

    # glm(formula = log10(x.mean + model.offset) ~ log10(size), family = inverse.gaussian, data = df.var.sum )

    fit.list <- list(fit = glm(formula = log10(x.mean + model.offset) ~ log10(size), family = gaussian, data = df.var.sum ),
                     model = "gaussian")
    return(fit.list)
  })

  if (mean.model$model == "gaussian"){
    df.var$mean.pred <- ( 10^predict.glm(object = mean.model$fit, newdata =df.var, type="response")) - model.offset
    df.var.sum$mean.pred <-( 10^predict.glm(object = mean.model$fit, newdata =df.var.sum, type="response")) - model.offset
    linear.model <- list(
      fit = mean.model,
      model.offset = model.offset
    )
  } else if (mean.model$model == "spline") {
    df.var$mean.pred <- predict.glm(object = mean.model$fit, newdata =df.var, type="response")
    linear.model <- list(
      fit = mean.model,
      model.offset = 0
    )
  }


  df.var$pi.pred.high <- df.var$mean.pred +  2*(df.var$sd.pred)
  df.var$pi.pred.lo <-  df.var$mean.pred-2*(df.var$sd.pred)

  df.raw.score <- as.data.frame(t(scores))
  score.name <- colnames(df.raw.score)
  df.raw.score$size = rsize
  df.raw.score.long <- pivot_longer(data = df.raw.score, cols = score.name)

  yscale <- c(min(df.var$pi.pred.lo, na.rm = T), max(df.var$pi.pred.high, na.rm = T))
  if (yscale[1] <=0){
    yscale[1] <- yscale[1]*1.1
  } else {
    yscale[1] <- yscale[1]*0.9
  }

  if (yscale[2] <=0){
    yscale[2] <- yscale[2]*0.9
  } else {
    yscale[2] <- yscale[2]*1.1
  }

  plt.mean <- df.var %>%
    ggplot(aes(x = size, y = (x.mean))) +
    scattermore::geom_scattermore(data = df.raw.score.long, aes(x = size, y = value), color = "grey") +
    geom_point(pch = 21, color = "white", size = 2, fill = "black") +
    geom_line(aes(y = mean.pred),  color = "red", size = 1) +
    geom_line(aes(y = pi.pred.high), color = "red", size = 1)  +
    geom_line(aes(y = pi.pred.lo), color = "red", size = 1)   +
    labs(x = "Geneset Size", y = "Score Mean", title = "Null model with 95% CI",
         subtitle = paste0("Empirical False Positive Rate = ", 100*signif(1- sum(df.var$pi.pred.high > df.var$x.mean) / nrow(df.var), 3), "%")) +
    theme_miko(grid = T) +
    coord_cartesian(ylim = yscale)

  plt.null.correction <- df.var %>%
    ggplot(aes(x = size, y = (x.mean - mean.pred)/sd.pred , fill = (x.mean)/sd.pred > 2)) +
    # scattermore::geom_scattermore() +
    geom_point(pch = 21, color = "white", size = 2) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = 2) +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey")) +
    labs(color = "z > 2", x = "Geneset Size", y = "Corrected Score") + theme_miko(grid = T)

  list(
    variance.model = variance.model,
    variance.raw.plot = variance.raw.plot,
    variance.mean.plot = variance.mean.plot,
    mean.model = mean.model,
    mean.plot = plt.mean,
    corrected.plot = plt.null.correction,
    null.scores = scores,
    geneset  = rgene.list
  )


}





#' Calculate standardized module scores for feature expression programs in single cells.
#'
#' Calculates the standardized expression level of each program on single-cell level, by subtracting by aggregated expression of control feature sets (like in `AddModuleScore`) and then scaling the difference by the pooled variance of the gene set of interest and control features. Like `AddModuleScore`, all analyzed features are binned based on averaged expression and control features are randomly selected from each bin.
#'
#' @param object Seurat Object
#' @param features A list of vectors of features for expression programs; each entry should be a vector of feature names
#' @param pool List of features to check expression levels against, defaults to `rownames(x = object)`
#' @param nbin Number of bins of aggregate expression levels for all analyzed features
#' @param ctrl Number of control features selected from the same bin per analyzed feature
#' @param nworkers Number of workers for parallel implementation. Default is 1.
#' @param k Use feature clusters returned from DoKMeans
#' @param assay Name of assay to use
#' @param name Name for the expression programs; will append a number to the end for each entry in features (eg. if `features` has three programs, the results will be stored as name1, name2, name3, respectively)
#' @param seed Set a random seed. If NULL, seed is not set.
#' @param search Search for symbol synonyms for features in features that don't match features in object? Searches the HGNC's gene names database; see UpdateSymbolList for more details
#' @param ... Extra parameters passed to UpdateSymbolList
#' @name AddSModuleScore
#' @concept miko_score
#' @author Nicholas Mikolajewicz
#' @return Returns seurat object with standardized module scores added to object meta data; each module is stored as name# for each module program present in `features`
#' @seealso \code{\link{AddModuleScore}} for original module scoring function implemented in `Seurat`.
AddSModuleScore <- function (object, features, pool = NULL, nbin = NULL, ctrl = 100, nworkers = 1,
                             k = FALSE, assay = NULL, name = "Cluster", seed = 1, search = FALSE, ...){

  require(parallel)
  ncore <- nworkers
  if (is.null(names(features))) names(features) <- paste0("gs", 1:length(features))
  names(features) <- make.names(names(features))

  if (is.null(nbin)) nbin <- optimalBinSize(object)
  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }
  assay.old <- DefaultAssay(object = object)
  assay <- assay %||% assay.old
  DefaultAssay(object = object) <- assay
  assay.data <- GetAssayData(object = object)
  features.old <- features
  if (k) {
    .NotYetUsed(arg = "k")
    features <- list()
    for (i in as.numeric(x = names(x = table(object@kmeans.obj[[1]]$cluster)))) {
      features[[i]] <- names(x = which(x = object@kmeans.obj[[1]]$cluster ==
                                         i))
    }
    cluster.length <- length(x = features)
  } else {
    if (is.null(x = features)) {
      stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
      missing.features <- setdiff(x = x, y = rownames(x = object))
      if (length(x = missing.features) > 0) {

        if (search) {
          tryCatch(expr = {
            updated.features <- UpdateSymbolList(symbols = missing.features,
                                                 ...)
            names(x = updated.features) <- missing.features
            for (miss in names(x = updated.features)) {
              index <- which(x == miss)
              x[index] <- updated.features[miss]
            }
          }, error = function(...) {
            warning("Could not reach HGNC's gene names database",
                    call. = FALSE, immediate. = TRUE)
          })
          missing.features <- setdiff(x = x, y = rownames(x = object))

        }
      }
      return(intersect(x = x, y = rownames(x = object)))
    })
    cluster.length <- length(x = features)
  }

  LengthCheck <- function(values, cutoff = 0) {
    return(vapply(
      X = values,
      FUN = function(x) {
        return(length(x = x) > cutoff)
      },
      FUN.VALUE = logical(1)
    ))
  }

  if (!all(LengthCheck(values = features))) {
    features <- lapply(X = features.old, FUN = CaseMatch,
                       match = rownames(x = object))
  }
  if (!all(LengthCheck(values = features))) {
    stop(paste("The following feature lists do not have enough features present in the object:",
               paste(names(x = which(x = !LengthCheck(values = features)))),
               "exiting..."))
  }
  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]
  data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                         n = nbin, labels = FALSE, right = FALSE)
  names(x = data.cut) <- names(x = data.avg)
  ctrl.use <- vector(mode = "list", length = cluster.length)

  miko_message("Preparing control genesets...", ...)
  ctrl.use2 <-  pbapply::pblapply(seq_along(features), function(i){

    features.use <- features[[i]]

    for (j in 1:length(x = features.use)) {
      which.ind <- which(x = data.cut == data.cut[features.use[j]])
      if (length(which.ind) < ctrl){
        ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which.ind], size = ctrl, replace = T)))
      } else {
        ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which.ind], size = ctrl, replace = FALSE)))
      }

    }

    return(ctrl.use[[i]])
  })

  ctrl.use <- ctrl.use2; rm(ctrl.use2)

  ctrl.use <- lapply(X = ctrl.use, FUN = unique)

  # start cluster
  if (ncore > detectCores()){
    n.work.score <- detectCores()
  } else {
    n.work.score <- ncore
  }

  cl <- parallel::makeCluster(n.work.score)
  doParallel::registerDoParallel(cl)

  # miko_message("Scoring gene sets...", ...)
  # iterate through each gene signature list
  score.results <- foreach(i = 1:cluster.length, .packages = c("Matrix"))  %dopar% {

    features.use <- ctrl.use[[i]]

    ad <- assay.data[features.use, ]
    cscore <- Matrix::colMeans(x = ad)
    czscore <- sparseMatrixStats::colVars(x = ad)

    features.use <- features[[i]]
    data.use <- assay.data[features.use, , drop = FALSE]
    fscore <- Matrix::colMeans(x = data.use)
    fzscore <- sparseMatrixStats::colVars(x = data.use)

    return(
      list(
        cscore = cscore,
        czscore = czscore,
        fscore = fscore,
        fzscore = fzscore
      )
    )

  }

  # stop workers
  parallel::stopCluster(cl)


  ctrl.scores <- ((do.call(rbind, pbapply::pblapply(score.results, function(x){ x[["cscore"]]} ))))
  features.scores <- ((do.call(rbind, pbapply::pblapply(score.results, function(x){ x[["fscore"]]} ))))
  ctrl.var <- ((do.call(rbind, pbapply::pblapply(score.results, function(x){ x[["czscore"]]} ))))
  features.var <- ((do.call(rbind, pbapply::pblapply(score.results, function(x){ x[["fzscore"]]} ))))
  # features.scores.use <- (features.scores - ctrl.scores)
  features.var[features.var == 0] <- median(ctrl.var[ctrl.var != 0])
  ctrl.var[ctrl.var == 0] <- median(ctrl.var[ctrl.var != 0])

  # raw score
  raw.score <- features.scores - ctrl.scores
  # rownames(x = raw.score) <- paste0(name, 1:cluster.length)
  rownames(x = raw.score) <- paste0("raw_", names(features))
  raw.score <- as.data.frame(x = t(x = raw.score))
  rownames(x = raw.score) <- colnames(x = object)
  object@misc[["raw_score"]] <- raw.score

  # standardized score
  features.scores.use <- (features.scores - ctrl.scores)/sqrt(ctrl.var + features.var) #+ features.var
  # rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
  rownames(x = features.scores.use) <- paste0(name, names(features))
  features.scores.use <- as.data.frame(x = t(x = features.scores.use))
  rownames(x = features.scores.use) <- colnames(x = object)
  object[[colnames(x = features.scores.use)]] <- features.scores.use
  CheckGC()
  DefaultAssay(object = object) <- assay.old
  return(object)
}



#' Calculate Miko module scores for feature expression programs
#'
#' Feature expression programs are scored by first computing standardized module scores with `AddSModuleScore` and then scaling the scores by the null distributions calculated using `nullScore`. The resulting scores are robust to gene set size and can be further used to compute whether feature expression program is significantly up-regulated in single-cell population.
#'
#' @param object Seurat Object
#' @param geneset A list of vectors of features for expression programs; each entry should be a vector of feature names
#' @param nullscore `nullScore` output for provided `object`. Must run `nullScore` prior to running `mikoScore`.
#' @param assay Name of assay to use.
#' @param nworkers Number of workers for parallel implementation. Default is 1.
#' @param nbin Number of bins of aggregate expression levels for all analyzed features. See `AddModuleScore` for details.
#' @param verbose Print progress. Default is TRUE.
#' @name mikoScore
#' @author Nicholas Mikolajewicz
#' @concept miko_score
#' @return list of results.
#' @seealso \code{\link{AddSModuleScore}} for standardized module scoring, \code{\link{nullScore}} for calculating null score distributions \code{\link{sigScore}} for derivation of p values for miko scores.
#' @examples
#'
#' maxgs <- max(unlist(lapply(gene.list, length)))
#' if (maxgs > 200) maxgs <- 200
#' stepsize <- round(maxgs / 15)
#' if (stepsize < 1) stepsize <- 1
#' if (stepsize > 20) stepsize <- 10
#'
#' # Generating null distributions for genesets
#' bm.null <- nullScore(object = object, assay = "RNA, n.replicate = 10, nbin = 24,
#'                     min.gs.size = 2, max.gs.size = gs.size + 5, step.size = stepsize,
#'                     nworkers =20, verbose = T, subsample.n = 5000)
#'
#' # variance.raw.plot <- bm.null$variance.raw.plot
#' # variance.mean.plot <- bm.null$variance.mean.plot
#' # plt.null.model <- bm.null$mean.plot
#' # plt.corrected.plot <- bm.null$corrected.plot
#'
#' # Compute miko scores for genesets
#' object <- mikoScore(object = object, geneset = gene.list, nbin = 24,
#'                           nullscore = bm.null, assay = "RNA", nworkers = 20)
#'
#' # Get significant  miko scores for scored genesets
#' score.result <- sigScore(object = object,  geneset = gene.list, reduction = "umap")
#'
mikoScore <- function(object, geneset, nullscore, assay = DefaultAssay(object), nworkers = 1, nbin = 24 ,  verbose = T){

  # ncores <- nworkers
  nbin <- optimalBinSize(object = object, nbin = nbin)

  miko_message("Running cell scoring...", verbose = verbose)

  nonrelevant.field <- colnames(object@meta.data)[grepl("cell|cluster", colnames(object@meta.data))]

  names(geneset) <- make.names(names(geneset))

  if (nworkers > length(geneset)) nworkers <- length(geneset)


  object <-  AddSModuleScore(
    object = object,
    features = geneset,
    pool = NULL,
    nbin = nbin,
    ctrl = 100,
    k = FALSE,
    nworkers = nworkers,
    assay = assay,
    name = "cell_",
    seed = 1,
    search = F
  )


  # raw scores
  df.raw <- object@misc[["raw_score"]]
  if (is.numeric(df.raw)){
    df.ms <- as.data.frame(df.raw)
    rownames(df.raw) <- colnames(object)
  }
  colnames(df.raw) <- names(geneset)

  # raw.scores <- as.matrix(df.raw)
  object@misc[["raw_score"]] <-df.raw

  # cell scores
  df.ms <- object@meta.data[ ,grepl("cell", colnames(object@meta.data) )]
  av.field <- colnames(df.ms)
  nonrelevant.field2 <- nonrelevant.field[nonrelevant.field %in% av.field]
  try({df.ms <- df.ms %>% dplyr::select(-nonrelevant.field2)}, silent = T)

  if (is.numeric(df.ms)){
    df.ms <- as.data.frame(df.ms)
    rownames(df.ms) <- colnames(object)
  }
  colnames(df.ms) <- gsub("cell_", "",  colnames(df.ms))
  gs.scores <- as.matrix(df.ms)


  # clusters scores
  miko_message("Running cluster scoring...", verbose = verbose)
  df.gs.size <- data.frame(size = unlist(lapply(geneset, length)))


  if (nullscore$mean.model$model == "spline"){
    var.pred <- predict.glm(object = nullscore$variance.model, newdata = df.gs.size, type="response")
    mean.pred <- predict.glm(object = nullscore$mean.model$fit, newdata = df.gs.size, type="response")
    gs.score.correct <- (gs.scores - mean.pred)/sqrt(var.pred)
  } else if (nullscore$mean.model$model == "gaussian"){
    var.pred <- predict.glm(object = nullscore$variance.model, newdata = df.gs.size, type="response")
    mean.pred <- (10^(predict.glm(object = nullscore$mean.model$fit, newdata = df.gs.size, type="response"))) - nullscore$mean.model$model.offset
    gs.score.correct <- (gs.scores - mean.pred)/sqrt(var.pred)
  }
  gs.score.correct <- as.data.frame(gs.score.correct)

  colnames(x = gs.score.correct) <- paste0("cluster_", colnames(x = gs.score.correct))
  object[[colnames(x = gs.score.correct)]] <- gs.score.correct

  miko_message("Complete!", verbose = verbose)

  return( object)

}






#' Benchmark Miko scoring pipeline using cluster-specific markers that had been variably contaminated with random genes.
#'
#' Differentially-expressed gene sets for each cluster are derived and variably contaminated with random genes prior to calculating cluster-specific miko scores for each gene set. The performance of miko scoring is then evaluated on cluster-specific and non-specific gene sets.
#'
#' @param object Seurat Object
#' @param geneset.size size of differentially-expressed gene sets
#' @param group_by Name of grouping variable in `object` meta feature. Default is "seurat_clusters".
#' @param assay Name of assay to use
#' @param deg.logFC.threshold logFC threshold for differential-expression analysis. Default is 0.5.
#' @param deg.fdr.threshold FDR threshold for differential-expression analysis. Default is 0.05.
#' @param miko.fdr.threshold FDR threshold for Miko scores. Default is 0.05.
#' @param nworkers Number of workers for parallel implementation. Default is 1.
#' @param verbose Print progress. Default is TRUE.
#' @name benchmarkScores
#' @author Nicholas Mikolajewicz
#' @concept miko_score
#' @return list of benchmark results.
#' @seealso \code{\link{AddSModuleScore}} for standardized module scoring, \code{\link{wilcoxauc}} for differential expression analysis
benchmarkScores <- function(object, geneset.size = 15, group_by = "seurat_clusters", assay = DefaultAssay(object),deg.logFC.threshold = 0.5, deg.fdr.threshold = 0.05, miko.fdr.threshold = 0.05, verbose = T, nworkers = 1){

  miko_message("Performing differential expression analysis...", verbose = verbose)
  require(presto)
  df.prest <- presto::wilcoxauc(X = object, group_by = group_by, assay = "data", seurat_assay = assay)

  gs.size <- geneset.size
  df.presto.sig <- df.prest %>% dplyr::filter(logFC > deg.logFC.threshold, padj < deg.fdr.threshold)

  u.group <- as.numeric(unique(df.presto.sig$group))
  u.group <- u.group[order(u.group)]

  miko_message("Generating benchmark genesets...", verbose = verbose)
  contam.vector <-unique( round(seq(0, gs.size, by = gs.size / 10)))
  all.gene <- rownames(so.query)

  gene.list <- list()
  for (i in 1:length(u.group)){

    df.presto.sig.cur <- df.presto.sig %>% dplyr::filter(group == u.group[i])

    if (nrow(df.presto.sig.cur) < gs.size) next
    nonsig.gene <- all.gene[!c(all.gene %in% df.presto.sig.cur$feature)]
    for (j in 1:length(contam.vector)){

      n.sig <- gs.size - contam.vector[j]
      n.rand <- contam.vector[j]


      gene.list[[paste0("c", u.group[i], "_", contam.vector[j], "contam")]] <- c(sample(df.presto.sig.cur$feature, n.sig),
                                                                                 sample(nonsig.gene, n.rand))

    }
  }


  maxgs <- gs.size
  if (maxgs > 200) maxgs <- 200
  stepsize <- round(maxgs / 15)
  if (stepsize < 1) stepsize <- 1
  if (stepsize > 20) stepsize <- 10
  miko_message("Generating null distributions for benchmarking genesets...", verbose = verbose)
  bm.null <- nullScore(object = object, assay = assay, n.replicate = 10, nbin = 24,
                       min.gs.size = 2, max.gs.size = gs.size + 5, step.size = stepsize,
                       nworkers =nworkers, verbose = T, subsample.n = 5000)

  miko_message("Computing miko scores for benchmarking genesets...", verbose = verbose)
  object <- mikoScore(object = object, geneset = gene.list, nbin = 24,
                      nullscore = bm.null, assay = assay, nworkers = nworkers)

  miko_message("Getting significant  miko scores for benchmarking genesets...", verbose = verbose)
  score.result <- sigScore(object = object,  geneset = gene.list, reduction = "umap")

  miko_message("Computing coherence scores for benchmarking genesets...", verbose = verbose)
  raw.mat <- object@misc[["raw_score"]]
  colnames(raw.mat) <- gsub("raw_", "", colnames(raw.mat))

  df.cscore <- coherentFraction(object = object, score.matrix =raw.mat, nworkers = nworkers,
                              genelist = gene.list, assay = assay, slot = "data", subsample.cluster.n = 500)

  df.score <- score.result$cluster_stats
  df.score$gs <- df.score$name
  df.score$cluster <- df.score$cluster
  df.merge <- merge(df.cscore, df.score, by = c("gs", "cluster"))


  df.val <- df.merge
  df.val$cell.type <- df.val$gs
  df.val$target.cluster <- gsub("c", "", stringr::str_extract(df.val$cell.type, "c[0-9]*"))
  df.val$contam <- as.numeric(gsub("contam", "", gsub("_", "", stringr::str_remove(df.val$cell.type, "c[0-9]*"))))/gs.size

  df.val$is.target <- gsub("c", "", df.val$cluster) == df.val$target.cluster

  # overall summary plots #####################################################
  df.val.sum <- df.val %>%
    dplyr::group_by(is.target, contam) %>%
    dplyr::summarize(
      sig.enrich = sum(fdr < miko.fdr.threshold & miko_score > 0)/length(fdr),
      sig.enrich.80frac = sum( (fdr < miko.fdr.threshold) & (miko_score > 0) & (coherence_fraction >= 0.8)) / length(fdr),
      sig.enrich.90frac = sum( (fdr < miko.fdr.threshold) & (miko_score > 0) & (coherence_fraction >= 0.9)) / length(fdr),
      sig.enrich.100frac = sum( (fdr < miko.fdr.threshold) & (miko_score > 0) & (coherence_fraction >= 1)) / length(fdr),
      .groups = 'drop'
    )
  df.val.sum.long <- pivot_longer(df.val.sum, cols = c("sig.enrich", "sig.enrich.80frac", "sig.enrich.90frac",
                                                       "sig.enrich.100frac"))


  df.val.sum.long$is.target2 <- "cluster non-specific"
  df.val.sum.long$is.target2[df.val.sum.long$is.target] <- "cluster-specific"
  plt.ce.performance <- df.val.sum.long %>%
    ggplot(aes(x = contam*100, y = value, color = name)) +
    geom_point() +
    coord_cartesian(ylim = c(0,1)) +
    geom_line() +
    facet_wrap(~is.target2) +
    theme_miko(legend = T, grid = T, color.palette = "ptol") +
    labs(y = "Fraction of genesets that are significant", x = "Geneset contamination (%)", color = "scoring")



  return(
    list(
      benchmark_data = df.val,
      benchmark_plot = plt.ce.performance,
      benchmark_statistics=  df.val.sum.long,
      benchmark_genesets = gene.list
    )

  )

}




#' Examine cluster-specific scoring performance for cluster-specific gene sets that have been variably contaminated with random genes.
#'
#' Examine cluster-specific scoring performance for cluster-specific gene sets that have been variably contaminated with random genes.
#'
#' @param object Seurat Object
#' @param benchmark_data `benchmark_data` returned in list from `benchmarkScores` function. Must run `benchmarkScores` prior to `benchmarkCluster`.
#' @param benchmark_genesets `benchmark_genesets` returned in list from `benchmarkScores` function. Must run `benchmarkScores` prior to `benchmarkCluster`.
#' @param which.cluster Name of group within grouping variable used to run `benchmarkScores`. E.g., '3' if `benchmarkScores(... , group_by = "seurat_clusters")`
#' @param miko.fdr.threshold FDR threshold for Miko scores. Default is 0.05.
#' @param coherence.fraction.threshold Coherence fraction threshold. Default is 0.8.
#' @name benchmarkCluster
#' @author Nicholas Mikolajewicz
#' @return list of benchmark results for specified cluster.
#' @seealso \code{\link{benchmarkScores}} for miko scoring benchmarking
#' @examples
benchmarkCluster <- function(object, benchmark_data, benchmark_genesets, which_cluster, miko.fdr.threshold = 0.05, coherence.fraction.threshold = 0){

  benchmark_data.cluster <- benchmark_data %>% dplyr::filter(cluster == which_cluster, is.target)
  val.gs <- benchmark_genesets[benchmark_data.cluster$cell.type]



  df.dot.all <- NULL
  for (i in 1:length(val.gs)){

    df.dot <- DotPlot(object = object, features = val.gs[[i]])[["data"]]
    df.dot <- df.dot %>% dplyr::filter(id == which_cluster)
    df.dot$gs = names(val.gs)[i]

    df.dot.all <- bind_rows(df.dot.all, df.dot)

  }

  df.dot.all$target.cluster <- gsub("c", "", stringr::str_extract(df.dot.all$gs, "c[0-9]*"))
  df.dot.all$contam <- as.numeric(gsub("contam", "", gsub("_", "", stringr::str_remove(df.dot.all$gs, "c[0-9]*"))))/unique(unlist(lapply(val.gs, length)))

  df.dot.all <- df.dot.all %>%
    dplyr::group_by(gs) %>%
    dplyr::arrange(-avg.exp.scaled) %>%
    dplyr::mutate(ind = seq(1, length(avg.exp.scaled)))

  plt.dot <- df.dot.all %>%
    ggplot(aes(x = ind, y = as.factor(round(contam*100)), color = avg.exp.scaled, size = pct.exp)) +
    geom_point() +
    scale_color_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
    labs(x = "Gene Index", y = "Gene Set Contamination (%)", title = "Gene Set Expression Profiles", subtitle = "Same cellular population, varying gene set contamination", size = "Percent Expr.", color = "Scaled Expr.") +
    theme_miko(legend = T, grid = T) +
    theme(legend.position = "bottom")


  plt.enrich <- benchmark_data.cluster %>%
    ggplot(aes(x = miko_score, y = as.factor(round(contam*100)), fill = fdr < miko.fdr.threshold & miko_score > 0))+
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey")) +
    theme_miko(legend = T, grid = T) +
    geom_bar(stat = "identity") +
    labs(x=  "Miko Score", y = "Gene Set Contamination (%)", title = "Miko Score", fill = paste0("FDR<", miko.fdr.threshold)) +
    theme(legend.position = "bottom")


  plt.coh <- benchmark_data.cluster %>%
    ggplot(aes(x = coherence_fraction, y = as.factor(round(contam*100)), fill =coherence_fraction  > coherence.fraction.threshold))+
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey")) +
    theme_miko(legend = T, grid = T) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = coherence.fraction.threshold, linetype = "dashed") +
    labs(x=  "Coherence Fraction", y = "Gene Set Contamination (%)", title = "Coherence Fraction", fill = paste0("coh.frac>", coherence.fraction.threshold)) +
    theme(legend.position = "bottom")


  plt.cluster.rep <- cowplot::plot_grid(plt.dot, plt.enrich, plt.coh, nrow = 1, align = "h")


  return(
    list(benchmark_cluster_plot = plt.cluster.rep,
         benchmark_data_cluster = benchmark_data.cluster)
  )

}




#' Calculate Miko score significance
#'
#' Miko scores are transformed to p values and corrected for multiple-comparisons using Benjamini Hochberg method.
#'
#' @param object Seurat Object. Must contain miko scores (i.e., must run `mikoScore` prior to `sigScore`)
#' @param geneset gene set list used for miko scoring (i.e., those provided as input into `mikoScore`).
#' @param reduction Dimensionality reduction to add to cell-level statistics data.frame. Intended for downstream use (e.g., viewing miko scores in UMAP space).
#' @name sigScore
#' @author Nicholas Mikolajewicz
#' @return list of data frames containing cell- and cluster-level statistics. Also return plot comparing raw (unstandardized) module score with Miko scores.
#' @seealso \code{\link{AddSModuleScore}} for standardized module scoring, \code{\link{mikoScore}} for miko scoring
#' @examples
#'
#' maxgs <- max(unlist(lapply(gene.list, length)))
#' if (maxgs > 200) maxgs <- 200
#' stepsize <- round(maxgs / 15)
#' if (stepsize < 1) stepsize <- 1
#' if (stepsize > 20) stepsize <- 10
#'
#' # Generating null distributions for genesets
#' bm.null <- nullScore(object = object, assay = "RNA, n.replicate = 10, nbin = 24,
#'                     min.gs.size = 2, max.gs.size = gs.size + 5, step.size = stepsize,
#'                     nworkers =20, verbose = T, subsample.n = 5000)
#'
#' # variance.raw.plot <- bm.null$variance.raw.plot
#' # variance.mean.plot <- bm.null$variance.mean.plot
#' # plt.null.model <- bm.null$mean.plot
#' # plt.corrected.plot <- bm.null$corrected.plot
#'
#' # Compute miko scores for genesets
#' object <- mikoScore(object = object, geneset = gene.list, nbin = 24,
#'                           nullscore = bm.null, assay = "RNA", nworkers = 20)
#'
#' # Get significant  miko scores for scored genesets
#' score.result <- sigScore(object = object,  geneset = gene.list, reduction = "umap")
#'
sigScore <- function(object, geneset, reduction = "umap"){

  names(geneset) <- make.names(names(geneset))
  df.cs <- object@meta.data[ ,paste0("cluster_", names(geneset))]
  df.rs <- object@misc[["raw_score"]]

  set.names <- colnames(df.cs)
  df.cs$cluster <- object@meta.data$seurat_clusters
  try({
    df.cs$reduction.x <- object@reductions[[reduction]]@cell.embeddings[,1]
    df.cs$reduction.y <- object@reductions[[reduction]]@cell.embeddings[,2]
  }, silent = T)

  df.cs.sum <- df.cs %>%  pivot_longer(cols = set.names)
  df.cs.sum$name <- gsub("cluster_", "", df.cs.sum$name)

  set.names <- colnames(df.rs)
  df.rs$cluster <- object@meta.data$seurat_clusters
  try({
    df.rs$reduction.x <- object@reductions[[reduction]]@cell.embeddings[,1]
    df.rs$reduction.y <- object@reductions[[reduction]]@cell.embeddings[,2]
  }, silent = T)
  df.rs.sum <- df.rs %>%
    pivot_longer(cols = set.names)
  df.cs.sum$value.rs <- df.rs.sum$value

  df.cs.sum <- as.data.frame(df.cs.sum)

  df.cs.sum2 <- df.cs.sum %>%
    dplyr::group_by( name, cluster) %>%
    dplyr::summarise(
      raw_score = mean(value.rs, na.rm = T),
      miko_score = mean(value, na.rm = T),
      .groups = 'drop'
    )

  df.cs.sum2$p <- 2*pnorm(q=abs(df.cs.sum2$miko_score), lower.tail=FALSE)

  df.cs.sum2 <- df.cs.sum2 %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(
      fdr = p.adjust(p, method = "BH")
    )


  plt.score.comp <- df.cs.sum2 %>%
    ggplot(aes(x = raw_score, y = miko_score, fill = fdr < 0.05 & miko_score > 0)) +
    scale_fill_manual(values = c("TRUE" = "tomato", "FALSE" = "grey")) +
    geom_point(pch = 21, color = "white", size = 2) +
    # scattermore::geom_scattermore() +
    theme_miko(legend = T, grid=  T) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(title = "Raw vs. Miko Enrichment Scores", x = "Raw Score", y = "Miko Score",
         subtitle = paste0(signif(sum(df.cs.sum2$fdr < 0.05 & df.cs.sum2$miko_score > 0)/nrow(df.cs.sum2), 3)*100, "% Significance Rate (FDR < 0.05)"),
         fill = "FDR<0.05")


  return(
    list(
      cluster_stats = df.cs.sum2,
      cell_stats = df.cs.sum,
      score_plot = plt.score.comp
    )

  )


}


#' Calculate coherent fraction for feature expression program.
#'
#' Calculate fraction of genes that are correlated with feature expression program. Performed on a per-cluster basis ("seurat_clusters" in `object` meta data).
#'
#' @param object Seurat Object.
#' @param score.matrix matrix of feature expression program. Can be computed using `AddModuleScores`, `AddSModuleScores`, `mikoScore`, among others.
#' @param geneset gene set list used for obtaining `score.matrix` (e.g., gene set provided as input into `mikoScore` or `AddSModuleScores`).
#' @param method Correlation method: "pearson" or "spearman". Default is "pearson".
#' @param assay Name of assay to use.
#' @param slot Use expression data from this slot in `object`.
#' @param subsample.cluster.n number of cells to subsample within each cluster ("seurat_clusters" in `object` meta data). Default is 500.
#' @param nworkers Number of workers for parallel implementation. Default is 1.
#' @param verbose Print progress. Default is TRUE.
#' @name coherentFraction
#' @author Nicholas Mikolajewicz
#' @return data.frame containing cluster-level coherent fractions.
#' @seealso \code{\link{AddSModuleScore}} for standardized module scoring, \code{\link{mikoScore}} for miko scoring
#' @examples
#'
#' so.query <-  AddSModuleScore(object = so.query, features = gene.list)
#'
#' raw.mat <- so.query@misc[["raw_score"]]
#' colnames(raw.mat) <- gsub("raw_", "", colnames(raw.mat))
#
#' df.cscore <- coherentFraction(object = so.query, score.matrix =raw.mat, nworkers = 20,
#'                               genelist = gene.list, assay = DefaultAssay(so.query), slot = "data", subsample.cluster.n = 500)
#'
#'
coherentFraction <- function(object, score.matrix, genelist, method = c("pearson", "spearman"), assay = DefaultAssay(object), slot = "data", subsample.cluster.n = 500, nworkers = 1, verbose = T){

  set.seed(1023)

  stopifnot("'object' is not a Seurat object" = ("Seurat" %in% class(object)))
  names(genelist) <- make.names(names(genelist))
  colnames(score.matrix) <- make.names(colnames(score.matrix))
  all.genes <- unique(unlist(genelist))

  expr.mat <- getExpressionMatrix(so = object, only.variable = F, which.data = slot, which.assay = assay)
  expr.mat <- expr.mat[rownames(expr.mat) %in% all.genes, ]
  expr.mat <- tryCatch({
    t(expr.mat)
  }, error = function(e){
    return(t(as.matrix(expr.mat)))
  })

  if ("seurat_clusters" %in% colnames(object@meta.data)){
    df.meta <- data.frame(cluster = object@meta.data$seurat_clusters)
  } else {
    object <- setResolution(object = object, resolution = 0.8)
    df.meta <- data.frame(cluster = object@meta.data$seurat_clusters)
  }

  uclust <- as.numeric(as.character(unique(df.meta$cluster))); uclust <- uclust[order(uclust)]

  # prerank
  if (method == "spearman"){
    miko_message("Preranking expression and score matrices...", verbose = verbose)
    expr.mat.rank <- pbapply::pbapply(expr.mat, 2, rank)
    rs.res.rank <- pbapply::pbapply(score.matrix, 2, rank) #,  ties.method = "random", na.last = NA
  } else if (method == "pearson"){
    expr.mat.rank <- as.matrix(expr.mat)
    rs.res.rank <- as.matrix(score.matrix)
  }


  df.cscore <- NULL
  miko_message("Computing cluster-level correlations between gene expression and scores...", verbose = verbose)

  # nworkers <- nworkers
  if (length(uclust) < nworkers) nworkers <- length(uclust)
  cl <- parallel::makeCluster(nworkers)
  doParallel::registerDoParallel(cl)


  # iterate through each gene signature list
  coh.results <- foreach(i = 1:length(uclust))  %dopar% { #

    # for (i in 1:length(uclust)){
    # print(i)

    which.cell <- which(df.meta$cluster == uclust[i])
    if (subsample.cluster.n < length(which.cell)){
      which.cell <- sample(which.cell, subsample.cluster.n)
    }

    if (length(which.cell) > 1){

      expr.mat.cur <- expr.mat.rank[which.cell, ]
      rs.res.cur <- rs.res.rank[which.cell, ]

      if (all(class(expr.mat.cur) == class(rs.res.cur))){
        cor.score <- cor(expr.mat.cur, rs.res.cur)
      }
      rownames(cor.score) <- colnames(expr.mat.cur)
      colnames(cor.score) <- colnames(rs.res.cur)
      cor.score[is.na(cor.score)] <- 0

      cor.gs <- lapply(seq_along(genelist), function(j){
        cor.score[rownames(cor.score) %in% genelist[[j]], names(genelist)[j]]
      })


      # suppressWarnings({
      #   suppressMessages({
      # cwct <- lapply(cor.gs, function(x){  wilcox.test(x = x, mu = 0, alternative = "greater")[["p.value"]]})

      cfrac <- lapply(cor.gs, function(x){ sum(x> 0, na.rm = T)/ length(x)})

      cn <- lapply(cor.gs, length)

      cmean <- lapply(cor.gs, median, na.rm = T)
      cvar <- lapply(cor.gs, var, na.rm = T)

      df.res <- data.frame(cluster = uclust[i],
                           gs = names(genelist),
                           coherence_fraction = unlist(cfrac),
                           r_mean = unlist(cmean),
                           r_variance = unlist(cvar),
                           # coh.sd = unlist(csd),
                           gene.n = unlist(cn),
                           cell.n = length(which.cell)
      )




    } else {
      df.res <- data.frame(cluster = uclust[i],
                           gs = names(genelist),
                           coherence_fraction = 0,
                           r_mean = 0,
                           r_variance = 0,
                           gene.n = unlist( lapply(genelist, length)),
                           cell.n = length(which.cell)
      )
    }

    return(df.res)

  }

  parallel::stopCluster(cl)

  df.cscore <- bind_rows(coh.results)

  return(df.cscore)

}



#' Word cloud visualization of cell-type annotations from the Miko scoring pipeline
#'
#' Given outputs from the Miko scoring pipeline, the top cell-type annotations fare visualized using word clouds.
#'
#' @param object Seurat Object.
#' @param object.group name of object meta data field specifying cluster membership. Default is "seurat_clusters".
#' @param score vector of Miko scores
#' @param score.group vector of group memberships
#' @param score.cell.type vector of cell-type names/labels.
#' @param score.p vector of p values.
#' @param score.fdr vector of fdr values. Optional.
#' @param score.coherence.fraction vector of coherence fractions. See coherentFraction(...) for details.
#' @param score.frequent.flier vector of logicals specifying whether score belongs to frequent flier.
#' @param fdr.correction Specify whether p-value should be corrected using Benjamini & Hochberg method. Default is T.
#' @param p.threshold p value threshold. Default is 0.05.
#' @param coherence.threshold Numerical [0,1] specifying minimal coherence required to qualify for visualization. Default is 0.8.
#' @param show.n.terms Maximal number of cell-type terms shown in word cloud. Default is 15.
#' @param verbose Logical, specify whether process is printed. Default is T.
#' @name annotationCloud
#' @author Nicholas Mikolajewicz
#' @return list of ggplot handles
#' @seealso \code{\link{mikoScore}} for miko scoring, \code{\link{coherentFraction}} for coherence scoring
#' @examples
#'
#' df.score_summary <- data.frame(cluster = df.merge$cluster,
#'                               cell.type = df.merge$gs,
#'                               miko_score = signif(df.merge$miko_score, 3) ,
#'                               p =  signif(df.merge$p),
#'                               fdr =  signif(df.merge$fdr),
#'                               coherence_fraction =  signif(df.merge$coherence_fraction))
#'
#'
#' plt.cloud <- annotationCloud(object = so.query_scored,
#'                              object.group = "seurat_clusters",
#'                              score = df.score_summary$miko_score,
#'                              score.group = df.score_summary$cluster,
#'                              score.cell.type = df.score_summary$cell.type,
#'                              score.p = df.score_summary$p,
#'                              score.fdr = df.score_summary$fdr,
#'                              score.coherence.fraction = df.score_summary$coherence_fraction,
#'                              score.frequenct.flier = NULL,
#'                              fdr.correction = T,
#'                              p.threshold = 0.05,
#'                              coherence.threshold = 0.9,
#'                              show.n.terms = 15,
#'                              verbose = T)
#'
annotationCloud <- function(object,
                            object.group = "seurat_clusters",
                            score,
                            score.group,
                            score.cell.type,
                            score.p,
                            score.fdr = NULL,
                            score.coherence.fraction = NULL,
                            score.frequenct.flier = NULL,
                            fdr.correction = T,
                            p.threshold = 0.05,
                            coherence.threshold = 0.8,
                            show.n.terms = 15,
                            verbose = T){


  require(ggwordcloud)

  miko_message("Generating annotation wordclouds...", verbose = verbose)

  if (is.null(score.coherence.fraction)) coherence.threshold <- 0
  if (is.null(score.frequenct.flier)) score.frequenct.flier <- F

  vst.merge.cloud <- data.frame(
    miko_score = score,
    cluster = score.group,
    cell.type = score.cell.type,
    p = score.p,
    fdr = score.fdr,
    coherence_fraction = score.coherence.fraction,
    frequent_flier = score.frequenct.flier
  )

  u.cl <- unique(vst.merge.cloud$cluster)
  plt.ww.list <- list()

  show.w <- show.n.terms


  plt.cluster.umap <- highlightUMAP(object = object, group = object.group,
                                    reduction = "umap", highlight.color = "tomato")
  names(plt.cluster.umap) <- as.character(gsub("group_", "",  names(plt.cluster.umap)))

  u.cl <- object@meta.data[[object.group]]
  #
  # # get unique clusters
  if (object.group == "seurat_clusters"){
    u.cl <- unique(as.numeric(as.character((u.cl))))
    u.cl <- u.cl[order(u.cl)]
  } else {
    u.cl <- unique((as.character((u.cl))))
  }



  for (i in 1:length(u.cl)){

    vst.merge.subset <- vst.merge.cloud[vst.merge.cloud$cluster %in% u.cl[i], ]
    vst.merge.subset$coh.score <- vst.merge.subset$coherence_fraction
    vst.merge.subset$coh.score[vst.merge.subset$coh.score < coherence.threshold] <- 0

    vst.merge.subset <- vst.merge.subset %>% dplyr::filter(miko_score > 0,
                                                           coherence_fraction >= coherence.threshold,
                                                           !frequent_flier)

    if (nrow(vst.merge.subset) >0){


      n.sig.score <- nrow(vst.merge.subset %>% dplyr::filter(fdr < p.threshold, miko_score > 0))
      n.coh.score <- nrow(vst.merge.subset %>% dplyr::filter(fdr < p.threshold,miko_score > 0,
                                                             coherence_fraction >= coherence.threshold)) #, coh.score

      if (fdr.correction){
        enrich.label <- paste0(n.sig.score, "/", length(marker.list), " (", signif(n.sig.score/  ulength(vst.merge.cloud$cell.type), 2)*100, "%) gene sets are enriched (FDR < 5e-2*, 1e-5**, 1e-8***)\nand exceed ", 100*coherence.threshold, "% coherence")
      } else {
        enrich.label <- paste0(n.sig.score, "/", length(marker.list), " (", signif(n.sig.score/  ulength(vst.merge.cloud$cell.type), 2)*100, "%) gene sets are enriched (p < 5e-2*, 1e-5**, 1e-8***)\nand exceed ", 100*coherence.threshold, "% coherence")
      }


      vst.merge.subset$sig.stringent <- vst.merge.subset$fdr < p.threshold &
        vst.merge.subset$coherence_fraction > coherence.threshold
      df.f1 <- vst.merge.subset %>% dplyr::filter(sig.stringent) %>%
        dplyr::filter( miko_score > 0, coherence_fraction >= coherence.threshold) %>% #| (coh.fdr < 0.05)
        dplyr::top_n(show.w, miko_score)
      if (nrow(df.f1) < show.w){
        ndif <- show.w - nrow(df.f1)
        df.f1 <- bind_rows(df.f1, vst.merge.subset %>%
                             dplyr::filter(!sig.stringent) %>%
                             dplyr::filter(miko_score > 0) %>% #| (coh.fdr < 0.05) #fdr < 0.05,
                             dplyr::top_n(ndif, miko_score))
      }
      if (nrow(df.f1) == 0) next

      df.f1$cell.type <- gsub("_|-", " ", df.f1$cell.type)
      df.f1$cell.type <- stringr::str_wrap(df.f1$cell.type, 40)

      maxlim <- max(df.f1$miko_score)
      if (maxlim > 30){
        maxlim <- 30
      }

      minlim.coh <-  min(df.f1$coh.score)
      maxlim.coh <- max(df.f1$coh.score)

      df.f1$miko_score[(df.f1$miko_score) > maxlim] <- maxlim
      df.f1$coh.score[(df.f1$coh.score) > maxlim.coh] <- maxlim.coh

      df.f1$cell.type2 <- df.f1$cell.type
      df.f1$cell.type2[df.f1$fdr < p.threshold] <- paste0( df.f1$cell.type[df.f1$fdr < p.threshold], "*")
      df.f1$cell.type2[df.f1$fdr < 1e-5] <- paste0( df.f1$cell.type[df.f1$fdr < 1e-5], "**")
      df.f1$cell.type2[df.f1$fdr < 1e-8] <- paste0( df.f1$cell.type[df.f1$fdr < 1e-8], "***")

      df.f1$cell.type3 <- df.f1$cell.type2

      limseq <- unique(round( seq(-log10(p.threshold), maxlim, by = signif((maxlim-(-log10(p.threshold)))/4, 1))))

      if (fdr.correction){
        col.label <- "-log10(FDR)"
      } else {
        col.label <- "-log10(p)"
      }

      w1 <- df.f1 %>%
        ggplot(aes(label = cell.type3, color = -log10(fdr), size = abs(miko_score))) +
        geom_text_wordcloud(scale_size_area = 40, rm_outside = TRUE, eccentricity = 1, show.legend = T) +
        theme_minimal() +
        labs(title  = "Enriched", subtitle = enrich.label,  size = "Score", color = col.label) +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)) +
        theme(legend.position = "bottom") +
        scale_color_gradient2(low = "grey", mid = "grey", high = scales::muted("red"), midpoint = -log10(p.threshold)) +
        scale_size(range = c(1, 5)) +
        theme( #
          legend.title = element_text(color = "black", size = 10),
          legend.text = element_text(color = "black", size = 8) ,
          legend.box.background = element_rect(colour = "black")
        )

      plt.ww.list[[as.character(u.cl[i])]] <- cowplot::plot_grid(plt.cluster.umap[[as.character(u.cl[i])]],w1, ncol = 2)

    }
  }
  return(plt.ww.list)
}


