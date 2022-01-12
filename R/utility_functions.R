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
#' @param object Seurat object
#' @param gNames.list named character vector mapping ensembl (names) to gene symbols (entries)
#' @param convert.RNA logical to convert ensembl in RNA assay to gene symbol
#' @name ens2sym.so
#' @author Nicholas Mikolajewicz
#' @return Seurat object with genes represented as symbols
#'
ens2sym.so <- function(object, gNames.list, convert.RNA = TRUE){


  # var features
  try({
    if ("SCT" %in% names(object@assays)){
      warning("Converting SYMBOL to ENSEMBLE in SCT assay...\n")
      so_temp <- object@assays[["SCT"]]@var.features
      if ( checkGeneRep (gNames.list, so_temp) == "ensembl") object@assays[["SCT"]]@var.features <- as.vector((gNames.list[so_temp]))
    }
  }, silent = T)


  # scale data
  try({
    if ("SCT" %in% names(object@assays)){
      so_temp <- object@assays[["SCT"]]@scale.data
      if (checkGeneRep (gNames.list, row.names(so_temp)) == "ensembl") {
        row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
        object@assays[["SCT"]]@scale.data <- so_temp
      }
    }
  }, silent = T)

  # data
  try({
    if ("SCT" %in% names(object@assays)){
      so_temp <- object@assays[["SCT"]]@data
      if ( checkGeneRep (gNames.list, row.names(so_temp)) == "ensembl") {
        row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
        object@assays[["SCT"]]@data <- so_temp
      }
    }
  }, silent = T)

  # metadata
  try({
    if ("SCT" %in% names(object@assays)){
      so_m <- object@assays[["SCT"]]@meta.features
      so_m$ENSEMBLE <- rownames(so_m)
      if (checkGeneRep (gNames.list,  so_m$ENSEMBLE) == "ensembl") {
        so_m$SYMBOL <- as.vector((gNames.list[so_m$ENSEMBLE]))
        rownames(so_m) <-  make.names(so_m$SYMBOL, unique = T)
        object@assays[["SCT"]]@meta.features <- so_m
      }
    }
  }, silent = T)

  # dimnames
  try({
    if ("SCT" %in% names(object@assays)){
      so_temp <- object@assays[["SCT"]]@counts@Dimnames[[1]]
      if (checkGeneRep (gNames.list, so_temp) == "ensembl") {
        object@assays[["SCT"]]@counts@Dimnames[[1]] <- as.vector((gNames.list[so_temp]))
      }
    }
  }, silent = T)

  # pca feature loading
  try({
    if ("pca" %in% names(object@reductions)){
      so_temp <-  object@reductions[["pca"]]@feature.loadings
      if (checkGeneRep (gNames.list, row.names(so_temp)) == "ensembl") {
        row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
        object@reductions[["pca"]]@feature.loadings <- so_temp
      }
    }
  }, silent = T)

  # vst
  try({
    if ("SCT" %in% names(object@assays)){
      if (checkGeneRep (gNames.list,  names(object@assays[["SCT"]]@misc[["vst.out"]][["genes_log_gmean_step1"]])) == "ensembl") {
        names(object@assays[["SCT"]]@misc[["vst.out"]][["genes_log_gmean_step1"]]) <- as.vector((gNames.list[names(object@assays[["SCT"]]@misc[["vst.out"]][["genes_log_gmean_step1"]])]))
      }
    }
  }, silent = T)

  try({
    if ("SCT" %in% names(object@assays)){
      if (!is.null(object@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]])){
        if (checkGeneRep (gNames.list,   object@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]]@Dimnames[[1]]) == "ensembl") {
          object@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]]@Dimnames[[1]] <- as.vector((gNames.list[ object@assays[["SCT"]]@misc[["vst.out"]][["umi_corrected"]]@Dimnames[[1]]]))
        }
      }
    }
  }, silent = T)

  if ("SCT" %in% names(object@assays)){
    try({
      if (checkGeneRep (gNames.list,  rownames(object@assays[["SCT"]]@misc[["vst.out"]][["gene_attr"]])) == "ensembl") {
        rownames(object@assays[["SCT"]]@misc[["vst.out"]][["gene_attr"]]) <- as.vector((gNames.list[rownames(object@assays[["SCT"]]@misc[["vst.out"]][["gene_attr"]])]))
      }
    }, silent = T)
  }

  if ("SCT" %in% names(object@assays)){
    try({
      so.temp <- rownames(object@assays[["SCT"]]@misc[["vst.out"]][["model_pars"]])
      if (checkGeneRep(gNames.list, so.temp) == "ensembl") {
        rownames(object@assays[["SCT"]]@misc[["vst.out"]][["model_pars"]]) <- as.vector((gNames.list[so.temp]))
      }
    }, silent = T)
  }

  if ("SCT" %in% names(object@assays)){
    try({
      so.temp <- rownames(object@assays[["SCT"]]@misc[["vst.out"]][["model_pars_fit"]])
      if (checkGeneRep (gNames.list,  so.temp) == "ensembl") {
        rownames(object@assays[["SCT"]]@misc[["vst.out"]][["model_pars_fit"]]) <- as.vector((gNames.list[so.temp]))
      }
    }, silent = T)
  }

  # RNA ASSAY
  if (convert.RNA == TRUE){
    warning("Converting SYMBOL to ENSEMBLE in RNA assay...\n")
    # var features
    try({
      so_temp <- object@assays[["RNA"]]@var.features
      if (length(so_temp) > 0){
        if ( checkGeneRep (gNames.list, so_temp) == "ensembl") {
          object@assays[["RNA"]]@var.features <- as.vector((gNames.list[so_temp]))
        }
      }
    }, silent = T)


    # data
    try({
      so_temp <- object@assays[["RNA"]]@data
      if (checkGeneRep (gNames.list, row.names(so_temp)) == "ensembl") {
        row.names(so_temp) <-   as.vector((gNames.list[ row.names(so_temp)]))
        object@assays[["RNA"]]@data <- so_temp
      }
    }, silent = T)

    # counts
    try({
      so_temp <- object@assays[["RNA"]]@counts
      if ( checkGeneRep (gNames.list, row.names(so_temp)) == "ensembl") {
        row.names(so_temp) <-   as.vector((gNames.list[ row.names(so_temp)]))
        object@assays[["RNA"]]@counts <- so_temp
      }
    }, silent = T)

  }

  if ("integrated" %in% names(object@assays)){
    warning("Converting SYMBOL to ENSEMBLE in Integrated assay...\n")

    # var features
    try({
      so_ens <- object@assays[["integrated"]]@var.features
      if (checkGeneRep (gNames.list, so_ens) == "ensembl") {
        object@assays[["integrated"]]@var.features <-as.vector((gNames.list[so_ens]))
      }
    }, silent = T)

    # scale data
    try({
      so_sd <- object@assays[["integrated"]]@scale.data
      if ( checkGeneRep (gNames.list,  row.names(so_sd)) == "ensembl") {
        row.names(so_sd) <-  as.vector((gNames.list[row.names(so_sd)]))
        object@assays[["integrated"]]@scale.data <- so_sd
      }
    }, silent = T)

    # data
    try({
      so_d <- object@assays[["integrated"]]@data
      if (checkGeneRep (gNames.list,  row.names(so_d)) == "ensembl") {
        row.names(so_d) <-  as.vector((gNames.list[row.names(so_d)]))
        object@assays[["integrated"]]@data <- so_d
      }
    }, silent = T)

  }

  # ensure dim names are correctly specified
  object <- updateDimNames(object)

  return(object)
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
#' @param object Seurat Object
#' @name fixBarcodeLabel
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#'
fixBarcodeLabel <- function (object){

  # merge CellType and Barcode, if necessary
  meta.data.names <- names(object@meta.data)

  if (("CellType" %in% meta.data.names) & ("Barcode" %in% meta.data.names)){
    if ("integrated" %in% names(object@assays)){
      barcode <- object@meta.data[["Barcode"]]
      celltype <- object@meta.data[["CellType"]]
      barcode[is.na(barcode)] <- celltype[is.na(barcode)]
    } else {
      barcode <- object@meta.data[["CellType"]]
    }
  } else if (!("CellType" %in% meta.data.names) & ("Barcode" %in% meta.data.names)) {
    barcode <- object@meta.data[["Barcode"]]
  } else if (("CellType" %in% meta.data.names) & !("Barcode" %in% meta.data.names)) {
    barcode <- object@meta.data[["CellType"]]

  } else {
    barcode <- object@meta.data[["orig.ident"]]
    # stop("Problem with CellType/Barcode metadata detected. Troubleshooting required")
    }

  object@meta.data[["Barcode"]] <- barcode

  return(object)
}

#' Set cluster resolution
#'
#' Set 'Seurat_Clusters' metadata entry to specified cluster resolution [0, inf]. Wrapper for Seurat's FindClusters() function.
#'
#' @param object Seurat Object
#' @param resolution Numeric [0, inf] specifying cluster resolution. Values [0.1,1] typically perform well.
#' @param assay Seurat assay to check for existing clustering at specified resolution.
#' @param use.existing.clusters Logical flag specifying whether to use existing clustering solution if it already exists for specified resolution.
#' @param ... additional arguments passed to Seurat::FindClusters(...)
#' @name setResolution
#' @seealso \code{\link{FindClusters}}
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#'
setResolution <- function(object, resolution, assay = DefaultAssay(object), use.existing.clusters = T, ...){

  if (!("Seurat" %in% class(object))) stop("input is not a seurat object")
  if (!("numeric" %in% class(resolution))) stop("resolution is not a numeric")

  target.entry <- paste0(assay, "_snn_res.", resolution)

  if ((target.entry %in% names(object@meta.data)) & use.existing.clusters){
    object@meta.data[["seurat_clusters"]] <- object@meta.data[[target.entry]]
    Idents(object) <- object@meta.data[["seurat_clusters"]]
  } else {

    my.assay <- DefaultAssay(object)
    if ("integrated" %in% names(object@assays)){
      DefaultAssay(object) <- "integrated"
      object <- tryCatch({
        miko_message("Computing clusters using integrated assay...")
        object <- FindClusters(object = object, resolution = resolution, verbose = F,...)
      }, error = function(e){
        print(e)
        miko_message("Finding neighbors...\n")


        pca.prop <- propVarPCA(object)
        target.pc <- max(pca.prop$pc.id[pca.prop$pc.cum_sum<0.9])+1
        object <- FindNeighbors(object, verbose = F, reduction = "pca", dims = 1:target.pc)
        miko_message("Computing clusters using integrated assay...")
        object <- FindClusters(object = object, resolution = resolution, verbose = F, ...)
        DefaultAssay(object) <- my.assay
        return(object)
      })
      DefaultAssay(object) <- my.assay
    } else {
      object <- tryCatch({
        miko_message(paste0("Computing clusters using ", assay, " assay..."))
        object <- FindClusters(object = object, resolution = resolution, verbose = F,...)
      }, error = function(e){
        print(e)



        if (!("pca" %in% names(object@reductions))){
          miko_message("PCA is missing, running necessary preprocessing...")
          if (length(VariableFeatures(object) == 0)){
            object <- FindVariableFeatures(object = object)
          }

          object <- tryCatch({
            RunPCA(object = object, verbose = F)
          }, error = function(e){
            object <- scNormScale(
              so = object,
              method = "NFS",
              vars2regress = NULL,
              enable.parallelization = F,
              n.workers = 1,
              max.memory = (20480 * 1024^2),
              variable.features.n = NULL,
              variable.features.rv.th = 1.3,
              return.only.var.genes = F,
              mean.cutoff = c(0.1, 8),
              dispersion.cutoff = c(1, Inf),
              conserve.memory = T,
              assay = DefaultAssay(object)
            )
            object <- RunPCA(object = object, verbose = F)
            return(object)
          })

        }

        pca.prop <- propVarPCA(object)
        target.pc <- max(pca.prop$pc.id[pca.prop$pc.cum_sum<0.9])+1

        miko_message("Finding neighbors...")
        object <- FindNeighbors(object, verbose = F, reduction = "pca", dims = 1:target.pc)
        miko_message(paste0("Computing clusters using ", assay, " assay..."))
        object <- FindClusters(object = object, resolution = resolution, verbose = F, ...)
        return(object)
      })

    }

  }

  return(object)
}



#' Prepare gene to ensemble conversion vector
#'
#' Prepare gene to ensemble conversion vector
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
prepGeneList <- function (so, global.enviroment, species = "Hs") {


  if (!exists("gNames.list")){
    gNames.list <- NULL

    try({

    if (("gNames.list_master" %in% global.enviroment)) {
      for (i in 1:length(gNames.list_master)) {
        gNames.list <- c(gNames.list, gNames.list_master[[i]])
      }
      gNames.df <- data.frame(n = gsub("\\..*", "", as.vector(names(gNames.list))),
                              g = as.vector(gNames.list))
      gNames.df <- unique(gNames.df)
      gNames.list <- as.vector(gNames.df$g)
      names(gNames.list) <- as.vector(gNames.df$n)
    } else {
      try({av.meta <- so@assays[["RNA"]]@meta.features}, silent = T)
      try({av.meta2 <- so.query@assays[["SCT"]]@meta.features}, silent = T)
      if (exists("av.meta") && all(c("SYMBOL", "ENSEMBL") %in% colnames(av.meta)) && (sum(is.na(av.meta$ENSEMBL)) == 0)) {
        gNames.list <- as.vector(av.meta$SYMBOL)
        names(gNames.list) <- as.vector(av.meta$ENSEMBL)
      } else if (exists("av.meta2") && all(c("SYMBOL", "ENSEMBL") %in% colnames(av.meta2)) && (sum(is.na(av.meta2$ENSEMBL)) == 0)){
        gNames.list <- as.vector(av.meta2$SYMBOL)
        names(gNames.list) <- as.vector(av.meta2$ENSEMBL)
      }
    }

    if (!(exists("gNames.list"))){

      my.rep <- (rownames(so))


      ens.sum <-  sum(grepl("ENS", my.rep))
      ens.mus.sum <-  sum(grepl("ENSMUS", my.rep))
      hi.cap.sum <-  sum(my.rep == toupper(my.rep))
      lo.cap.sum <-  sum(my.rep == firstup(my.rep))

      df.rep <-as.data.frame(t(data.frame(
        ens.sum = ens.sum,
        ens.mus.sum = ens.mus.sum,
        hi.cap.sum = hi.cap.sum,
        lo.cap.sum = lo.cap.sum
      )))


      which.rep <- rownames(df.rep)[which.max(df.rep[, 1])]

      if ((which.rep %in% c("ens.sum", "ens.mus.sum")) & (ens.sum > ens.mus.sum)) {
        db <- org.Hs.eg.db::org.Hs.eg.db
        rep <- "ensemble"
      } else if ((which.rep %in% c("ens.sum", "ens.mus.sum")) & (ens.sum <= ens.mus.sum)) {
        db <- org.Mm.eg.db::org.Mm.eg.db
        rep <- "ensemble"
      } else if (which.rep == "hi.cap.sum") {
        rep <- "symbol"
        db <- org.Hs.eg.db::org.Hs.eg.db
      } else if (which.rep == "lo.cap.sum") {
        db <- org.Mm.eg.db::org.Mm.eg.db
        rep <- "symbol"
      }

      if (rep == "symbol"){
        my.gene <- AnnotationDbi::select(db,
                                         keys = my.rep,
                                         columns = c("ENSEMBL", "SYMBOL"),
                                         keytype = "SYMBOL",
                                         multiVals = first)
      } else if (rep == "ensemble"){
        my.gene <- AnnotationDbi::select(db,
                                         keys = my.rep,
                                         columns = c("ENSEMBL", "SYMBOL"),
                                         keytype = "ENSEMBL",
                                         multiVals = first)
      }

      my.gene <- my.gene[complete.cases(my.gene), ]
      gNames.list <- as.vector(my.gene$SYMBOL)
      names(gNames.list) <- as.vector(my.gene$ENSEMBL)

    }

    }, silent = T)

  }

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



#' Subset Seurat Object
#'
#' Subset Seurat object according to specified meta data field. Only specified meta data entries are retained, while remaining of data is omitted.
#'
#'
#' @param object Seurat object
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
subsetSeurat <- function (object, subset.df){

  # assertion
  if (class(object) != "Seurat") stop("'object' must be Seurat Object")

  # check if subset input is validd
  if (is.na(unique(subset.df$field))){
    subset.flag <- FALSE
    warning("Invalid subsetting field provided...\n")
  } else if ( unique(subset.df$field) %in% names(object@meta.data)) {
    subset.flag <- TRUE
  } else {
    subset.flag <- FALSE
    warning("Subsetting failed...\n")
  }


  # subset data
  if (subset.flag){
    require("Seurat", quietly = T)
    cur.field <- as.vector(unique(subset.df$field))
    if (cur.field == "seurat_clusters"){
      my.cells <- colnames(object)[(as.character(object@meta.data[[cur.field]]) %in% as.vector(subset.df$subgroups))]
    } else {
      pattern <- paste( as.vector(subset.df$subgroups), collapse="|")
      pattern <- gsub(" ", "", pattern)
      match.ind.1 <- grepl(pattern, as.character(object@meta.data[[cur.field]]))
      match.ind.2 <- as.character(object@meta.data[[cur.field]]) %in% as.vector(subset.df$subgroups)
      if (sum(match.ind.1) != sum(match.ind.2)) warning("exact vs. partial matching results in inconsistent number of matches...\n")
      my.cells <- grep(pattern, as.character(object@meta.data[[cur.field]]))
    }

    object <- subset(x = object, cells = my.cells)
    try({object <- UpdateSeuratObject(object)}, silent = T)
  }

  return(object)
}


# collections (subcollections):
# H: Hallmark
# C1: positional genesets
# C2: curated genesets (CGP, CP)
# C3: regulatory target genesets (MIR, TFT)
# C4: computational genesets (CGN, CM)
# C5: ontology genesets (GO, GO:BP, GO:CC, GO:MF, HPO)
# C6: oncogenic signature genesets
# C7: immunologic signatures
# C8: cell type signature genesets

#' Returns list of annotations for given Entrez gene IDs
#'
#' Returns list of Reactome or GO annotations for given Entrez gene IDs
#'
#' @param query.genes Entrez IDs of query genes
#' @param db Database to retrieve annotations from. One of:
#' \itemize{
#' \item "Bader" - Default. List of pathway annotations curated by Bader lab (http://baderlab.org/GeneSets)
#' \item "Reactome"
#' \item "GO" - Requires additional specification of ontology.
#' \item "msigdb" - Requires additional specification of msigdb.collection, and optionally msigdb.subcollection.
#' }
#' @param ontology GO ontologies to retrieve if GO db is selected. One of:
#' \itemize{
#' \item "BP" - Default. Biological processes
#' \item "MF" - Molecular functions
#' \item "CC" - Cellular components
#' }
#' @param species "Mm" or "Hs". Default is "Hs".
#' @param msigdb.collection Geneset collection (only if msigdb database used). See msigdbr() for additional details. One of:
#' \itemize{
#' \item "H" - Hallmark (Default)
#' \item "C1" - positional genesets
#' \item "C2" - curated genesets (CGP, CP)
#' \item "C3" - regulatory target genesets (MIR, TFT)
#' \item "C4" - computational genesets (CGN, CM)
#' \item "C5" - ontology genesets (GO, GO:BP, GO:CC, GO:MF, HPO)
#' \item "C6" - oncogenic signature genesets
#' \item "C7" - immunologic signatures
#' \item "C8" - cell type signature genesets
#' }
#' @param msigdb.subcollection Subcollection corresponding to specified msigdb collection. Possible subcollections for each collection are indicated in parantheses above.
#' @name getAnnotationPathways
#' @return Named list of vectors with gene sets (Entrez format).
#'
getAnnotationPathways <- function(query.genes, db = c("Bader"), ontology = c("BP"), species = c("Hs"),
                                  msigdb.collection = "H", msigdb.subcollection = NULL){

  if (db == "GO"){

    which.ontology <- ontology

    if (species == "Hs"){
      require("org.Hs.eg.db", quietly = T)
      go.e2g <- org.Hs.egGO
      go.g2e <- as.list(org.Hs.egGO2EG)
    } else if (species == "Mm"){
      require("org.Mm.eg.db", quietly = T)
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

  } else if (db == "msigdb"){

    if (is.null(msigdb.collection)) stop("'msigdb.collection' must be specified if using msigdb pathways\n")
    if (!(require(msigdbr))) {
      BiocManager::install("msigdbr");
      library(msigdbr, quietly = T)
    }

    if (species == "Hs"){
      msd.species <- "Homo sapiens"
    } else if (species == "Mm"){
      msd.species <- "Mus musculus"
    }

    df.msd <- msigdbr(species = msd.species, category = msigdb.collection)

    if (!is.null(msigdb.subcollection)) df.msd <-df.msd[grepl(msigdb.subcollection, df.msd$gs_subcat), ]
    if (nrow(df.msd) == 0) df.msd <- msigdbr(species = msd.species, category = msigdb.collection)

    df.msd.subset <- df.msd[df.msd$entrez_gene %in% query.genes, ]

    u.paths <- unique(df.msd$gs_name)

    pathways <- list()
    for (i in 1:length(u.paths)){
      pathways[[u.paths[i]]] <- df.msd.subset$entrez_gene[df.msd.subset$gs_name %in% u.paths[i]]
    }
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
#' @param which.assay Seurat assay to get data from. Default is DefaultAssay(so).
#' @param which.data Specify which data to use (refers to slots in Seurat object assay). One of:
#' \itemize{
#' \item "scale" - Default
#' \item "data"
#' }
#' @param use.additional.genes Character vector of additional genes to include (in addition to varibale, if variable flag is specificed). Default is NA.
#' @param as.dense Logical to convert sparse to dense matrix. Only applies if which.data is 'data'. Default is FALSE.
#' @name getExpressionMatrix
#' @author Nicholas Mikolajewicz
#' @return gene x cell expression matrix
#'
getExpressionMatrix <- function(so, only.variable = F, which.assay = NULL, which.data = "scale", use.additional.genes = NA, as.dense = F){

  # specify assay
  if (is.null(which.assay)) which.assay <- DefaultAssay(so)

  # get complete matrix
  if (which.data == "scale"){
    exp.mat.complete <- so@assays[[which.assay]]@scale.data
  } else if (which.data == "data"){
    if (as.dense){
      exp.mat.complete <- as.matrix(so@assays[[which.assay]]@data)
    } else {
      exp.mat.complete <- (so@assays[[which.assay]]@data)
    }

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
    input.dir <- "C:/Users/Owner/Dropbox/PDF Projects - JM/Data/scRNA-seq/01_sci-RNA-seq3_Hong_Kevin_Jason/NM_HH/Data/Reference_Datasets/"
    output.dir <- "C:/Users/Owner/Dropbox/PDF Projects - JM/R Packages/scMiko/data/"
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
#' @param object Seurat Object
#' @param subsample.factor Numeric [0,1]. Factor to downsample data by.
#' @param subsample.n Numeric [1,ncol(object)]. Number of cells to subsample. If specified, overides subsample.factor.
#' @param sample.group Character. Meta data grouping variable in which min.group.size will be enforced.
#' @param min.group.size Numeric [1,ncol(object)]. Minimum number of cells to downsample to within sample.group. If there are insufficient cells to achieve the target min.group.size, only the available cells are retained.
#' @param seed seed for sampling. Default is 1023.
#' @param verbose Print progress. Default is TRUE.
#' @name downsampleSeurat
#' @author Nicholas Mikolajewicz
#' @return Seurat Object
#'
downsampleSeurat <- function(object, subsample.factor = 1, subsample.n = NULL,  sample.group = NULL, min.group.size = 500, seed = 1023, verbose = T){


  if (is.null(seed)){
    set.seed(1023)
  } else {
    set.seed(seed)
  }


  use.factor <- (subsample.factor < 1) & (subsample.factor > 0 )
  use.n <-( !is.null(subsample.n)) & (is.numeric(subsample.n))

  n.subset <- NULL
  if (!(use.factor | use.n)){
    return(object)
  } else if (use.n & use.factor){
    if (subsample.n > ncol(object)) {
      subsample.n <- ncol(object)
    }
    n.subset <- round(subsample.n)
  } else if (use.n){
    if (subsample.n > ncol(object)) {
      subsample.n <- ncol(object)
    }
    n.subset <- round(subsample.n)
  } else if (use.factor){
    n.subset <- round(subsample.factor *ncol(object))
  } else {
    return(object)
  }

  if (!is.null(sample.group) && (sample.group %in% colnames(object@meta.data))){

    cell.ind <- sample(x = seq(1, ncol(object)), size = n.subset, replace = FALSE, prob = NULL)
    df.tally.orig <- as.data.frame(table(object@meta.data[ ,sample.group]))
    df.tally.sampled <- as.data.frame(table(object@meta.data[cell.ind ,sample.group]))
    colnames(df.tally.orig) <-  c("group", "n_original")
    colnames(df.tally.sampled) <- c("group", "n_sampled")


    df.tally <- merge(df.tally.orig, df.tally.sampled, by = "group")

    df.tally$dif <- df.tally$n_original - df.tally$n_sampled
    df.tally$flip <- (df.tally$n_original < min.group.size) != (df.tally$n_sampled < min.group.size)
    df.tally$below.target <- df.tally$n_sampled < min.group.size
    df.tally$nextra <- df.tally$n_sampled - min.group.size

    n.spare <- sum(df.tally$nextra[df.tally$nextra > 0])
    n.missing <- abs(sum(df.tally$nextra[df.tally$nextra < 0]))

    n.av <- min(c(n.spare, n.missing))

    df.meta <- data.frame(ind = seq(1, ncol(object)), group = object@meta.data[ ,sample.group])
    for (i in 1:nrow(df.tally)){
      df.tally.sampled <- as.data.frame(table(object@meta.data[cell.ind ,sample.group]))
      # n.spare <- sum(df.tally$nextra[df.tally$nextra > 0])
      n.spare <- df.tally$nextra[i]
      n.av <- min(c(n.spare, n.missing))
      if (n.av > 0){
        df.meta$spare <- (df.meta$ind %in% cell.ind) & (df.meta$group %in% df.tally$group[i] )
        df.meta$missing <- (!(df.meta$ind %in% cell.ind)) & (df.meta$group %in% df.tally.sampled$Var1[df.tally.sampled$Freq < min.group.size] )

        if (sum(df.meta$spare) == 0) next
        if (sum(df.meta$missing) == 0) next

        omit.ind <- sample(x = df.meta$ind[df.meta$spare], size = n.av, replace = FALSE, prob = NULL)
        n.av <- min(c(n.av, length(df.meta$ind[df.meta$missing] )))
        add.ind <- sample(x = df.meta$ind[df.meta$missing], size = n.av, replace = FALSE, prob = NULL)

        cell.ind <- unique(c(cell.ind[!(cell.ind %in% omit.ind)], add.ind))

        n.missing <- n.missing - length(add.ind)

        df.tally.sampled <- as.data.frame(table(object@meta.data[cell.ind ,sample.group]))
      }
    }


  } else {
    cell.ind <- sample(x = seq(1, ncol(object)), size = n.subset, replace = FALSE, prob = NULL)
  }


  if (!is.null(n.subset)){
    object <-  tryCatch({
      miko_message(paste0("Sampling ", n.subset, "/", ncol(object), " (", signif(100*n.subset/ncol(object), 4), "%) cells"), verbose = verbose)

      object <- subset(object , cells = cell.ind)
    }, error = function(e){
      warning("Failed to downsample seurat object. Troubleshooting required.")
      print(e)
      return(object)
    })

    return(object)
  } else {
    return(object)
  }

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
aggGroupExpression <-  function(so, which.data = "data", which.assay = DefaultAssay(so), which.center = "mean", which.group = "seurat_clusters", which.features = NULL, do.parallel = F, verbose = T){
  return(avgGroupExpression(so, which.data, which.assay, which.center, which.group, which.features, do.parallel, verbose))
}

#' Get summary of group expression in Seurat object
#'
#' Get summary group expression in Seurat object. Can include mean, median, fraction (of expressing cells), sd, or cv.
#'
#' @param so Seurat Object
#' @param which.data Character specifying which data slot. Default is "data".
#' @param which.assay Character specifying which assay to use.
#' @param which.center Character specifying which summary measure to use. Must be one of "mean", "median", "fraction", "sum", "sd", or "cv". If unspecified, default is "mean".
#' @param which.group Character specifying group field in Seurat metadata. Default is "seurat_clusters".
#' @param which.features Character specifying which genes to include. All if unspecified.
#' @param do.parallel Logical specifying whether to perform computations in parallel. Default is F. Uses future.apply package.
#' @param verbose Logical
#' @name avgGroupExpression
#' @author Nicholas Mikolajewicz
#' @return data.frame (gene rows, group columns)
#'
avgGroupExpression <-  function(so, which.data = "data", which.assay = DefaultAssay(so), which.center = "mean",
                                which.group = "seurat_clusters", which.features = NULL, do.parallel = F, verbose = T){
  # which.center options: "mean", "fraction", "median", "sum", "sd", "cv"

  # inititate parallel processes
  if (do.parallel){
    library(future.apply, quietly = T)
    plan(multisession) ## Run in parallel on local computer
  }

  # entire matrix
  exp.mat.complete <- getExpressionMatrix(so, which.data = which.data, which.assay = which.assay)

  # subset matrix
  if (!is.null(which.features)){
    which.features <- which.features[which.features %in% rownames(exp.mat.complete)]
    if (length(which.features) == 0) stop("'which.features' are not available in Seurat object")
    exp.mat.complete <- exp.mat.complete[rownames(exp.mat.complete) %in% which.features, ]
  }

  if ((is.null(dim(exp.mat.complete))) && (length(which.features) == 1)){
    exp.mat.complete <- t(as.matrix(exp.mat.complete))
    rownames(exp.mat.complete) <- which.features
  }

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
  rm(so); invisible({gc()})

  # compute measure of centrality
  avg.mat <- matrix(nrow = length(gene.list), ncol = length(u.clusters))
  if (verbose) miko_message(paste0("Computing ", which.center, "...") , verbose = verbose)
  for (i in 1:length(u.clusters)){

    current.mat <- exp.mat.complete[ ,cluster.membership %in% u.clusters[i]]
    # if (is.numeric(current.mat)){
    #   if (which.center %in% c("mean", "median")){
    #     avg.mat[,i] <- current.mat
    #   } else {
    #     avg.mat[,i] <- NA
    #   }
    # } else {

    if (which.center == "mean"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(current.mat, 1, function(x) log(mean(expm1(x), na.rm = T)+1))
      } else {
        avg.mat[,i] <- apply(current.mat, 1, function(x) log(mean(expm1(x), na.rm = T)+1))
      }
    } else if (which.center == "median"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(current.mat, 1, function(x) log(median(expm1(x), na.rm = T)+1))
      } else {
        avg.mat[,i] <- apply(current.mat, 1, function(x) log(median(expm1(x), na.rm = T)+1))
      }
    } else if (which.center == "sum"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(current.mat, 1, function(x) sum(x, na.rm = T))
      } else {
        avg.mat[,i] <- apply(current.mat, 1, function(x) sum(x, na.rm = T))
      }
    } else if (which.center == "sd"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(current.mat, 1, function(x) log(sd(expm1(x), na.rm = T)+1))
      } else {
        avg.mat[,i] <- apply(current.mat, 1, function(x) log(sd(expm1(x), na.rm = T)+1))
      }
    } else if (which.center == "cv"){
      if (do.parallel){
        sd.cur <- future_apply(current.mat, 1, function(x) log(sd(expm1(x), na.rm = T)+1))
        av.cur <- future_apply(current.mat, 1, function(x) log(mean(expm1(x), na.rm = T)+1))
      } else {
        sd.cur <- apply(current.mat, 1, function(x) log(sd(expm1(x), na.rm = T)+1))
        av.cur <- apply(current.mat, 1, function(x) log(mean(expm1(x), na.rm = T)+1))
      }
      avg.mat[,i] <- sd.cur / abs(av.cur)
    } else if (which.center == "fraction"){
      e.subset <- exp.mat.complete[ ,cluster.membership %in% u.clusters[i]]
      if (is.numeric(e.subset)){
        avg.mat[,i] <- mean(e.subset>0)
      } else {
        if (do.parallel){
          avg.mat[,i] <- future_apply(e.subset, 1, function(x) sum(x>0)/length(x))
        } else {
          avg.mat[,i] <- apply(e.subset, 1, function(x) sum(x>0)/length(x))
        }
      }
    } else {
      stop("which.center must be specified as 'mean', 'median', 'fraction', 'sd', or 'cv'")
    }
    # }
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
#' @param verbose Print progress. Default is TRUE
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
runWGCNA <- function(e.mat, s.mat = NULL, cor.metric = "rho_p", soft.power = 2, use.TOM = T, network.type = "signed", TOM.type = "unsigned", rescale.adjacency = F, verbose = T, ...){

  # similarity matrix - using proportionality metric for scRNAseq data.
  if (is.null(s.mat)){
    miko_message("Computing similarity matrix...", verbose = verbose)
    s.mat <-  dismay::dismay(e.mat, metric = cor.metric)
  }

  # adjacency matrix
  miko_message("Computing adjacency matrix...", verbose = verbose)
  a.mat <-  sim2adj(s.mat, soft.power, network.type)

  # rescale value if needed
  if (rescale.adjacency) a.mat <- recaleValues(a.mat, new.min = 0, new.max = 1)

  # compute topological overlap matix (TOM)
  if (use.TOM){
    miko_message("Computing topological overlap matix...", verbose = verbose)
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
  miko_message("Computing dissimilarity matix...", verbose = verbose)
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
  library(flashClust, quietly = T)

  tree <- flashClust(as.dist(d.mat), method = method, ...)

  return(tree)
}


#' Cut clustered tree at varying heights and overlay with dendrogram to find optimal parameter set.
#'
#' WGCNA::cutreeDynamic is run for varying levels of deepSplit parameter (0:4), and cluster membership is assigned for each parameter set. Default cut method is 'hybrid'.
#'
#' @param tree h.clust object generated by dist2hclust.
#' @param d.mat distance matrix used to generate tree.
#' @param genes vector of gene names corresponding rows/col of distance matrix (d.mat). If specified, additional "genes" column is provided in output.
#' @param cut.height Maximum joining heights that will be considered. Default is 0.998.
#' @param pam.respects.dendro Logical, only used for method "hybrid". If TRUE, the PAM stage will respect the dendrogram in the sense that objects and small clusters will only be assigned to clusters that belong to the same branch that the objects or small clusters being assigned belong to. Default is FALSE.
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
optimalDS <- function(tree, d.mat, genes = NULL, cut.height = 0.998, pam.respects.dendro = FALSE,  ...){

  require(dynamicTreeCut, quietly = T)
  require(WGCNA, quietly = T)

  mColorh <- NULL
  for (ds in 0:4){
    cut.tree <- cutreeDynamic(dendro = tree,distM= d.mat, cutHeight = cut.height, deepSplit=ds, pamRespectsDendro = pam.respects.dendro, ...)
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

  library(WGCNA, quietly = T)
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

  is.success <- F

  try({

    # set lower bound to zero
    old.min <- min(values, na.rm = T)
    if (old.min < 0) {
      values <- values + abs(old.min)
    } else if (old.min > 0) {
      values <- values - abs(old.min)
    }

    # set upper bound to one
    old.max <- (max(values, na.rm = T))
    values <- values/old.max

    new.range <- new.max - new.min
    if (new.range != 0){
      values <- values * new.range
      values <- values + new.min
    }

    if(min(values) != new.min) warning(paste0("minimum of rescaled values is not ", new.min, "\n"))
    if(max(values) != new.max) warning(paste0("maximum of rescaled values is not ", new.max, "\n"))

    is.success <- T
  }, silent = T)

  if (!is.success) values <- new.max * values / max(values, na.rm = T)

  return(values)

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
#' @param r2.target Retrieves soft power that corresponds to network topology corresponding to r2 >= r2.target. Default is 0.9.
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

getSoftThreshold2 <- function(s.mat, power =c(seq(0.5,5, by = 0.5), seq(6,10)), network.type = "signed", nBreaks = 20, removeFirst = T, rescale.adjacency = F, n.cores = 4, r2.target = 0.9){


  plt.sf.list <- list()
  powers <- power
  r2.sf <- c()

  # start cluster
  cl <- parallel::makeCluster(n.cores)
  doParallel::registerDoParallel(cl)

  sf.list <- list()

  sf.list <- (foreach(i = 1:length(powers), .packages = c("scMiko", "dplyr", "ggplot2", "ggpmisc"))) %dopar% {

    power.cur <- powers[i]
    a.cur <-  sim2adj(s.mat, soft.power = power.cur, network.type = network.type)

    if (rescale.adjacency)  a.cur <- recaleValues(a.cur, new.min = 0, new.max = 1)

    # get connectivity for specified power
    net.connectivity.df <- getConnectivity(a.cur, rownames(a.cur), flag.top.n = 20)

    rm(a.cur); invisible({gc()});

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

    r2.sf.current <- summary(lm1)[["r.squared"]] * sign(lm1[["coefficients"]][["log.dk"]])

    df.sf <- data.frame(x = log.dk, y = log.p.dk)

    # store node linkage distribution plot
    plt.sf.current <- df.sf %>%
      ggplot(aes(x=x, y=y)) +
      geom_smooth(method = "lm", color = "tomato", fill = "tomato", formula = 'y ~ x') +
      geom_point(size = 3) +
      xlab("N Links (Log)") +
      ylab("N Nodes (Log)") +
      labs(title = "Node Linkages", subtitle = paste0("Soft Power = ", power.cur)) +
      theme_classic() +
      stat_fit_glance(method = "lm",
                      label.y = "bottom",
                      method.args = list(formula = y ~ x),
                      mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',
                                                    stat(r.squared), stat(p.value))),
                      parse = TRUE)



    sf.list <- list(
      r2 = r2.sf.current,
      df = df.sf,
      plt = plt.sf.current
    )

    return(sf.list)

  }

  parallel::stopCluster(cl)

  # unpack results
  r2.sf <- c(); plt.sf.list <- list();
  for (i in 1:length(sf.list)){
    power.cur <- powers[i]
    r2.sf[i] <- sf.list[[i]][["r2"]]
    plt.sf.list[[as.character(power.cur)]] <- sf.list[[i]][["plt"]]
  }

  # store powers and r2
  df.r2.sf <- data.frame(sf = powers, r2 = r2.sf)

  # get best power estimate
  r2.target <- -1*abs(r2.target) # ensure correct sign
  if (sum(df.r2.sf$r2 < r2.target) > 0){
    best.power <- df.r2.sf$sf[which(df.r2.sf$r2 < r2.target)[1]]
    r2.opt <- df.r2.sf$r2[which(df.r2.sf$r2 < r2.target)[1]]
  } else {
    best.power <-  df.r2.sf$sf[which.min(df.r2.sf$r2)]
    r2.opt <- df.r2.sf$r2[which.min(df.r2.sf$r2)]
  }


  # optimizatio plot
  plt.opt.sf <- df.r2.sf %>%
    ggplot(aes(x = sf, y = r2)) +
    geom_hline(yintercept = r2.opt, color = "tomato") +
    geom_vline(xintercept = best.power, color = "tomato") +
    geom_smooth(method = "loess", color = "black", fill = "grey", formula = 'y ~ x') +
    geom_point(size = 3) +
    xlab("Soft Power") +
    ylab("R2 (Scale Free Topology)") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(title = "Soft Power Optimization",
         subtitle = paste0("Optimal power = ", best.power, ", r2 = ", signif(r2.opt, 3)))

  # store results
  output <- list(
    powerEstimate = best.power,
    r2Estimate =r2.opt,
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
#' Remove duplicate genes from Seurat Object.
#'
#' @param object Seurat Object
#' @param retain.graphs Logical specifying whether to retain graphs in subset Seurat object because Seurat::subset function inherently removed the graphs. Default is TRUE.
#' @name rmDuplicateGenes
#' @return Seurat object
#' @author Nicholas Mikolajewicz
#' @examples
#'
#' so.query <- rmDuplicateGenes(so.query, retain.graphs = T)
#'
rmDuplicateGenes <- function(object, retain.graphs = T){

  # get gene names
  all.genes <- rownames(object)

  # store graphs
  # subset function has unexpected behavior of removing graphs when subseting
  if (retain.graphs){
    graph.holder <- object@graphs
  }

  # if duplicates exist, subset seurat object
  if (sum(duplicated(all.genes)) > 0){
    which.unique <- all.genes[!duplicated(all.genes)]
    object <- subset(object, features = which.unique)
  }

  # retrieve graphs
  if (retain.graphs){
    object@graphs <- graph.holder
    rm(graph.holder)
  }

  return(object)

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
#' @param y Optional second gene set. If provided, resulting matrix is gene.sets x y.
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
jaccardSimilarityMatrix <- function(gene.sets, y = NULL, assert.unique = T){

  if (is.null(y)){

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


  } else {
    n.sets.x <- length(gene.sets)
    n.sets.y <- length(y)
    j.mat <- matrix(nrow = n.sets.x, ncol = n.sets.y)
    for (i in 1:n.sets.x){
      for (j in 1:n.sets.y){
        i.name <- names(gene.sets)[i]
        j.name <- names(y)[j]
        j.mat[i, j] <- scMiko::getJaccard(gene.sets[[i.name]], y[[j.name]], assert.unique = assert.unique)
      }
    }

    rownames(j.mat) <- names(gene.sets)
    colnames(j.mat) <- names(y)
  }

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
    dplyr::group_by(cluster) %>%
    dplyr::summarize(x.center = get.center(x, which.center),
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


  require("parsnip", quietly = T)
  require("yardstick", quietly = T)

  # get data for highly variable genes (hvg)
  cur.data <- Seurat::GetAssayData(so, slot = slot, assay = assay)
  # cur.data <-  as.matrix(as.data.frame(cur.data))
  match.ind <- which(rownames(cur.data) %in% hvg)
  # dat_use <- as.data.frame(t(Seurat::GetAssayData(so, slot = slot, assay = assay)[match.ind,]))
  is.success <- F
  try({
    dat_use <- as.data.frame(t(cur.data[match.ind,]))
    is.success <- T
  }, silent = T)

  if (!is.success){
    dat_use <- as.data.frame(t(as.matrix(cur.data[match.ind,])))
  }

  rm(cur.data)

  # merge expression data and pseudotime
  dat_use_df <- cbind(pseudotimes, dat_use)
  colnames(dat_use_df)[1] <- "pseudotime"
  dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

  # Define training, testing and validation sets
  colnames(dat_use_df) <- make.names(colnames(dat_use_df) , unique = TRUE, allow_ = TRUE)
  dat_split <- rsample::initial_split(dat_use_df)
  dat_train <- rsample::training(dat_split)
  dat_val <- rsample::testing(dat_split)

  # Train Model
  model <- parsnip::rand_forest(mtry =mtry, trees = trees, min_n = min_n, mode = mode) %>%
    set_engine("ranger", importance =importance, num.threads = num.threads) %>%
    fit(pseudotime ~ ., data = dat_train)

  # Evaluate Model
  val_results <- dat_val %>%
    dplyr::mutate(estimate = predict(model, .[,-1]) %>%
                    pull()) %>%
    dplyr::select(truth = pseudotime, estimate)
  model.metrics <- yardstick::metrics(data = val_results, truth, estimate)

  # store results
  df.metrics <- data.frame(lineage = lineage.name,
                           rmse = signif(model.metrics[[".estimate"]][1],3),
                           rsq = signif(model.metrics[[".estimate"]][2],3),
                           mae = signif(model.metrics[[".estimate"]][3],3))

  rm(model.metrics);
  invisible({gc()})

  output <- list(
    model = model,
    prediction = val_results,
    performance = df.metrics
  )

  return(output)


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





#' Variance explained by each principal component.
#'
#' For given Seurat Object, retrieve principal components and compute proportion of explained variance.
#'
#' @param so Seurat Object
#' @param reduction.name Name of reduction to use. Default is 'pca'.
#' @return data.frame summarize proportion of variance explained by each principal component
#' @author Nicholas Mikolajewicz
#' @name propVarPCA
propVarPCA <- function(so, reduction.name = "pca"){

  # get pca reduction
  if (reduction.name %in% names(so@reductions)){
    pc.std <- so@reductions[[reduction.name]]@stdev
  } else {
    if (reduction.name == "pca"){

      if (length(VariableFeatures(so) == 0)){
        so <- FindVariableFeatures(object = so)
      }

      so <- tryCatch({
        RunPCA(object = so, verbose = F)
      }, error = function(e){
        so <- scNormScale(
          so = so,
          method = "NFS",
          vars2regress = NULL,
          enable.parallelization = F,
          n.workers = 1,
          max.memory = (20480 * 1024^2),
          variable.features.n = NULL,
          variable.features.rv.th = 1.3,
          return.only.var.genes = F,
          mean.cutoff = c(0.1, 8),
          dispersion.cutoff = c(1, Inf),
          conserve.memory = T,
          assay = DefaultAssay(so)
        )
        so <- RunPCA(object = so, verbose = F)
        return(so)
      })

      pc.std <- so@reductions[[reduction.name]]@stdev

    } else {
      stop("No valid reduction available.")
    }
  }


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
    # original
    # c.dim <- length(so@assays[[all.assays[i]]]@counts@Dimnames)
    # d.dim <- length(so@assays[[all.assays[i]]]@data@Dimnames)

    # alternative
    c.dim <- 2-sum(unlist(lapply(so@assays[[all.assays[i]]]@counts@Dimnames, is.null)))
    d.dim <- 2-sum(unlist(lapply(so@assays[[all.assays[i]]]@counts@Dimnames, is.null)))

    if ((c.dim == 1) & (d.dim == 2)){
      if ((so@assays[[all.assays[i]]]@counts@Dim[1] == so@assays[[all.assays[i]]]@data@Dim[1]) &
          (so@assays[[all.assays[i]]]@counts@Dim[2] == so@assays[[all.assays[i]]]@data@Dim[2])){
        so@assays[[all.assays[i]]]@counts@Dimnames <- so@assays[[all.assays[i]]]@data@Dimnames
      }
    } else if ((c.dim == 2) & (d.dim == 1)){
      if ((so@assays[[all.assays[i]]]@counts@Dim[1] == so@assays[[all.assays[i]]]@data@Dim[1]) &
          (so@assays[[all.assays[i]]]@counts@Dim[2] == so@assays[[all.assays[i]]]@data@Dim[2])){
        so@assays[[all.assays[i]]]@data@Dimnames <- so@assays[[all.assays[i]]]@counts@Dimnames
      }
    } else if ((c.dim == 1) & (d.dim == 1)){
      meta.features <- so@assays[[all.assays[i]]]@meta.features
      if ("SYMBOL" %in% colnames(meta.features)){
        if (dim(so@assays[[all.assays[i]]]@counts)[1] == nrow(meta.features)){
          so@assays[[all.assays[i]]]@counts@Dimnames[[1]] <- as.character(meta.features$SYMBOL)
        }
        if (dim( so@assays[[all.assays[i]]]@data)[1] == nrow(meta.features)){
          so@assays[[all.assays[i]]]@data@Dimnames[[1]] <- as.character(meta.features$SYMBOL)
        }
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
  try({
    so <- updateDimNames(so)
  }, silent = T)


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
  # subset(x = so, cells = WhichCells(so, cells = include.which.all))
  # object <- subset(x = so, cells = WhichCells(so, cells = include.which.all))
  # SubsetData(so, cells = WhichCells(so, cells = include.which.all))
  return(so[ ,colnames(so) %in% include.which.all])
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
scMikoUpdate <- function(token = "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", ...){

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
#' @param input.data Input dataset (that was preprocessed in current module). A character.
#' @param input.subset data subset. A character specifying whether a certain subset of data was analyzed. A character.
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
updateCentralLog <- function(Module, input.data = NA, input.subset = NA, clog.file = "moduleLog.csv", log.path = if(exists("data.path")) data.path else "", user.id = if(exists("user")) user else "guest", run.notes = NA, pdf.flag = save.pdf){

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
    Identifier = paste0("R", current.run,"_", Module, "_", user),
    Module = Module,
    User = user.id,
    Date = format(Sys.time(), '%d %B, %Y'),
    Input = paste(input.data, collapse = ", "),
    Subset = paste(input.subset, collapse = ", "),
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
  warning("Central log update successful!\n")

  return(current.clog$Identifier)
}


#' Save figure as pdf
#'
#' Save figure as pdf
#'
#' @param file.name path/filename of output pdf. A character.
#' @param plot.handle plot handle that will be saved as pdf. E.g., GGplot handle.
#' @param fig.width Numeric for plot width. Default is 5.
#' @param fig.height Numeric for plot height Default is 5.
#' @param save.flag Logical to save figure. Default is TRUE.
#' @author Nicholas Mikolajewicz
#' @name savePDF
#' @examples
#'
#' # save module figure
#' savePDF(file.name = "M01_QC_violin.pdf", plot.handle = plt.QC_violin, fig.width = 5, fig.height = 5, save.flag = save.pdf)
#'
savePDF <- function(file.name, plot.handle, fig.width = 5, fig.height = 5, save.flag = T){
  if (save.flag){
    if (!grepl(".pdf", file.name)) file.name <- paste0(file.name, ".pdf")
    pdf(file=file.name, width = fig.width, height = fig.height)
    print(plot.handle)
    garbage <- dev.off()
  }
}

#' Save figure as html
#'
#' Save figure as html
#'
#' @param file.name path/filename of output html A character.
#' @param plot.handle plot handle that will be saved as html E.g., GGplot handle.
#' @param fig.width Numeric for plot width. Default is 5.
#' @param fig.height Numeric for plot height Default is 5.
#' @param save.flag Logical to save figure. Default is TRUE.
#' @author Nicholas Mikolajewicz
#' @name saveHTML
#' @examples
#'
#' # save module figure
#' saveHTML(file.name = "M01_QC_violin.pdf", plot.handle = plt.QC_violin, fig.width = 5, fig.height = 5, save.flag = save.pdf)
#'
saveHTML <- function(file.name, plot.handle, fig.width = 5, fig.height = 5, save.flag = T){

  if (save.flag){

    if (!grepl(".html", file.name)) file.name <- paste0(file.name, ".html")

    library(R3port, quietly = T);
    library(rmarkdown, quietly = T)

    try({
      R3port::html_plot(
        plot = plot.handle,
        out = file.name,
        title = "",
        titlepr = "",
        footnote = "",
        pwidth = fig.width,
        pheight = fig.height,
        res = NULL,
        fontsize = 12,
        units = "in",
        rawout = NULL,
        cleancur = T,
        show = F
      )
    }, silent = T)

    suppressMessages({
      suppressWarnings({
        self_contained_html <- function (input, output){
          require(xfun, quietly = T)
          input <- normalizePath(input)
          if (!file.exists(output))
            file.create(output)
          output <- normalizePath(output)
          template <- tempfile(fileext = ".html")
          on.exit(unlink(template), add = TRUE)
          write_utf8("$body$", template)
          from <- if (pandoc_available("1.17"))
            "markdown_strict"
          else "markdown"
          pandoc_convert(input = input, from = from, output = output,
                         options = c("--self-contained", "--template", template))
          invisible(output)
        }

        self_contained_html(input = file.name, output = file.name)

      })
    })

    # clean header
    library(xml2, quietly = T)
    h <- as_list(read_html(file.name))
    if (h[["html"]][["body"]][["p"]][[1]] == "<!DOCTYPE html>  "){
      h[["html"]][["body"]][["p"]][[1]] <- NULL
      catch_out <- write_html(as_xml_document(h$html), file.name,
                              options=c("format","no_declaration"))
    }
  }


}



#' Convert sparse matrix to dense matrix
#'
#' Convert sparse matrix to dense matrix. Developed to handle large datasets, by constructing dense matrix block-by-block.
#'
#' @param mat.sparse Sparse matrix.
#' @param block.size If large dataset, construct dense matrix column-wise block-by-block, where block size is specified by numeric block.size parameter. Default is 10000.
#' @param transpose Logical specifying whether to transpose data.
#' @param verbose.error Logical specifying whether to print error message in case of large dataset which requires block-by-block construction.
#' @author Nicholas Mikolajewicz
#' @name sparse2dense
#' @return dense matrix
#' @seealso \code{\link{sparse2df}}
#' @examples
#'
#' # Get sparse matrix
#' exp.mat <- so@assays[[current.assay]]@data
#'
#' # convert to dense matrix
#' mat.dense <- sparse2dense(mat.sparse, transpose = T)
#'
sparse2dense <- function(mat.sparse, block.size = 10000, transpose = F, verbose.error = F){

  mat.dense <- tryCatch({
    mat.dense <- as.matrix(mat.sparse)
    if (transpose) mat.dense <- t(mat.dense)
  }, error = function(e){
    if (verbose.error) print(e)
    mat.dense <- NULL
    n.blocks <- round(ncol(mat.sparse) / block.size)+1
    for (i in 1:(n.blocks)){
      start.range <- ((i-1)*block.size)+1
      end.range <- i*block.size
      if (start.range > ncol(mat.sparse)) next
      if (end.range > ncol(mat.sparse)) end.range <- ncol(mat.sparse)
      if (transpose){
        mat.dense <- rbind(mat.dense, t(as.matrix(mat.sparse[, start.range:end.range])))
      } else {
        mat.dense <- cbind(mat.dense, as.matrix(mat.sparse[, start.range:end.range]))
      }
    }
    if (transpose){
      colnames(mat.dense) <- rownames(mat.sparse)
      rownames(mat.dense) <- colnames(mat.sparse)
    } else {
      colnames(mat.dense) <- colnames(mat.sparse)
      rownames(mat.dense) <- rownames(mat.sparse)
    }

    return(mat.dense)
  })

  if (transpose){
    colnames(mat.dense) <- rownames(mat.sparse)
    rownames(mat.dense) <- colnames(mat.sparse)
  } else {
    colnames(mat.dense) <- colnames(mat.sparse)
    rownames(mat.dense) <- rownames(mat.sparse)
  }
  rm(mat.sparse); invisible({gc()})

  return(mat.dense)
}


#' Convert sparse matrix to data.frame
#'
#' Convert sparse matrix to data.frame. Developed to handle large datasets, by constructing data.frame block-by-block.
#'
#' @param mat.sparse Sparse matrix.
#' @param block.size If large dataset, construct data.frame column-wise block-by-block, where block size is specified by numeric block.size parameter. Default is 10000.
#' @param transpose Logical specifying whether to transpose data.
#' @param verbose.error Logical specifying whether to print error message in case of large dataset which requires block-by-block construction.
#' @author Nicholas Mikolajewicz
#' @name sparse2df
#' @return data frame
#' @seealso \code{\link{sparse2dense}}
#' @examples
#'
#' # Get sparse matrix
#' exp.mat <- so@assays[[current.assay]]@data
#'
#' # convert to data.frame
#' df <- sparse2df(mat.sparse, transpose = T)
#'
sparse2df <- function(mat.sparse, block.size = 10000, transpose = F, verbose.error = F){

  df <- tryCatch({
    df <- as.data.frame(mat.sparse)
    if (transpose) {
      df <- as.data.frame(t(df))
    }
  }, error = function(e){
    if (verbose.error) print(e)
    df <- NULL
    n.blocks <- round(ncol(mat.sparse) / block.size)+1
    for (i in 1:(n.blocks)){
      start.range <- ((i-1)*block.size)+1
      end.range <- i*block.size
      if (start.range > ncol(mat.sparse)) next
      if (end.range > ncol(mat.sparse)) end.range <- ncol(mat.sparse)
      if (transpose){
        df <- bind_rows(df,as.data.frame(t(as.matrix(mat.sparse[, start.range:end.range]))))
      } else {
        df <- bind_cols(df,as.data.frame(mat.sparse[, start.range:end.range]))
      }
    }
    if (transpose){
      colnames(df) <- rownames(mat.sparse)
      rownames(df) <- colnames(mat.sparse)
    } else {
      colnames(df) <- colnames(mat.sparse)
      rownames(df) <- rownames(mat.sparse)
    }
    return(df)
  })

  if (transpose){
    colnames(df) <- rownames(mat.sparse)
    rownames(df) <- colnames(mat.sparse)
  } else {
    colnames(df) <- colnames(mat.sparse)
    rownames(df) <- rownames(mat.sparse)
  }
  rm(mat.sparse); invisible({gc()})

  return(df)

}


#' Sort factor levels in numerical order
#'
#' Sort factor levels in numerical order
#'
#' @param factor factor with numerical entries.
#' @author Nicholas Mikolajewicz
#' @name orderedFactor
#' @return ordered factor
#'
orderedFactor <- function(f){
  f <- factor(f, levels = unique(f)[order(as.numeric(unique(as.character(f))))])
  return(f)
}



#' Get local density (z) of bivariate relationship (x,y)
#'
#' Get local density (z) of bivariate relationship (x,y)
#'
#' @param x numeric (equal length as y)
#' @param y numeric
#' @param ... additional parameters passed to MASS::kde2d
#' @author Nicholas Mikolajewicz
#' @name getDensity
#' @return density values
#' @author Kamil Slowikowski (https://slowkow.com/notes/ggplot2-color-by-density/)
#' @examples
#' # get data and compute densities
#' df.meta <- so@meta.data # so is seurat object
#' df.meta$density1 <- getDensity(df.meta$nCount_RNA, df.meta$percent.mt, n = 100)
#'
#' # generate scatter plot with overlayed density values
#' plt.handle1 <- df.meta %>% ggplot(aes(x = nCount_RNA, y = percent.mt, color = density1)) + geom_point() + theme_miko(legend = T) +
#'    xlab("UMI/cell") + ylab("Mitochondrial Content (%)") +
#'    labs(title = paste0("r = ", rho1p, "; rho = ", rho1s)) + scale_color_viridis("Density") +
#'    theme(legend.position="bottom", legend.key.width=unit(legend.width,"cm"))
#'
getDensity <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}



#' Parallelized correlation
#'
#' Parallelized correlation
#'
#' @param max matrix
#' @param method Correlation method, a character. Options: 'pearson' or 'spearman'.
#' @param do.par Logical indicating whether to use parallel implementation. Default is T.
#' @param n.workers number of workers to use for parallel computation.
#' @name parCor
#' @value correlation matrix
#' @examples
#' @author Nicholas Mikolajewicz
#'
parCor <- function(mat, method = "spearman", do.par = T, n.workers = 4){

  # check if eligible method
  stopifnot(method %in% c("spearman", "pearson"))

  # get necessary packages
  require("doParallel", quietly = T)
  require("data.table", quietly = T)
  require("WGCNA", quietly = T)

  # pre-rank matrix (if spearman)
  if (method == "spearman"){
    warning("pre-ranking matrix...\n")
    mat <- apply(mat, 2, function(x) rank(x))
  }

  # parallelize
  if (do.par){

    # specify n.workers and initiate clusters
    warning("initiating clusters...\n")
    cl <- parallel::makeCluster(n.workers)
    doParallel::registerDoParallel(cl)

    # chunk sizes
    chunk.ranges <- seq(1, ncol(mat), by = round(ncol(mat) / n.workers) +1)
    if (max(chunk.ranges) < ncol(mat)) chunk.ranges <- c(chunk.ranges, ncol(mat))

    # compute correlations
    warning("computing correlations...\n")
    res <- foreach(i = 1:(length(chunk.ranges)-1),
                   .combine = rbind,
                   .multicombine = TRUE,
                   .inorder = FALSE,
                   .packages = c('data.table', 'doParallel', 'WGCNA')) %dopar% {
                     start.ind <- chunk.ranges[i]
                     if (i == (length(chunk.ranges)-1)){
                       end.ind <- chunk.ranges[i+1]
                     } else {
                       end.ind <- chunk.ranges[i+1] - 1
                     }

                     WGCNA::cor(x = mat[,start.ind:end.ind], y = mat, method = 'pearson')
                   }

    stopCluster(cl)

  } else {

    # fast correlation implementation (no parallelized)
    fast_cor <- function(m) {
      m <- t(m)
      m <- m - rowMeans(m)           # center
      m <- m / sqrt(rowSums(m^2))    # scale
      tcrossprod(m)                  # cross-product
    }

    # compute correlation
    warning("computing correlations...\n")
    res <- fast_cor(mat)

  }

  return(res)

}


#' Run gene-set enrichment analysis (GSEA)
#'
#' Run gene-set enrichment analysis (GSEA)
#'
#' @param gene Character vector of gene names
#' @param value Numeric vector of values used to rank genes (e.g., logFC, expression, etc.). Must be same length as gene argument.
#' @param species Species. One of "Mm" (mouse) or "Hs" (human)
#' @param my.entrez Data.frame with SYMBOL and ENTREZID columns. Used to relate gene symbols to Entrez IDs. Retrieved internally unless explicitly provided.
#' @param my.pathway Named list of pathways, with each entry containing vector of Entrez IDs. Retrieved internally unless explicitly provided.
#' @param min.size Minimum gene set size. Default is 3.
#' @param max.size Minimum gene set size. Default is 300.
#' @param do.plot Logical to return dotplot visualizing top GSEA ranked pathways. Default is T.
#' @param plot.top.n Numeric specifying how many top pathways to visualize. Default is 10.
#' @param path.name.wrap.width Numeric specifying width of pathway names to wrap around. Argument passed to stringr::str_wrap(..., width = path.name.wrap.width)
#' @return list of enrichment results
#' @author Nicholas Mikolajewicz
#'
runGSEA <- function(gene, value, species, db = "GO", my.entrez = NULL, my.pathway = NULL, min.size = 3,
                    max.size = 300, do.plot = T, plot.top.n = 10, path.name.wrap.width = 40){


  suppressMessages({
    suppressWarnings({


      if (is.null(my.entrez)){
        # get entrez to gene symbol mapping
        my.symbol <-gene
        my.entrez <- sym2entrez(my.symbol, my.species = species )
        my.entrez <- my.entrez[complete.cases(my.entrez), ]
      }

      if (is.null(my.pathway)){
        # get pathways
        my.pathway <- getAnnotationPathways(query.genes = my.entrez$ENTREZID, db = db, ontology = "BP", species = species)
      }

      # prep genelist
      gene.list <- value
      names(gene.list) <- gene
      match.ind <- match(names(gene.list), my.entrez$SYMBOL)
      names(gene.list) <- as.character(my.entrez$ENTREZID[match.ind])
      gene.list <- gene.list[!is.na(names(gene.list))]
      gene.list = sort(gene.list, decreasing = TRUE)

      # clean list
      df.ent <- data.frame(names = names(gene.list), values = as.vector(gene.list))
      df.ent <- df.ent[complete.cases(df.ent), ]
      df.ent <- df.ent[!is.infinite(df.ent$values), ]
      gene.list.clean <- df.ent$values
      names(gene.list.clean) <- df.ent$names

      # pathway gsea enrichment

      # fgseaMultilevel {fgsea}
      gse.pathway <- fgsea::fgsea(my.pathway, gene.list.clean, nperm=1000, minSize = min.size, maxSize=max.size)

      # make human readable
      gse.pathway <- gse.pathway
      gse.pathway$set <- lapply(gse.pathway$leadingEdge,
                                mapvalues,from = my.entrez$ENTREZID, to = my.entrez$SYMBOL)
      gse.pathway$set <- lapply(gse.pathway$set, paste,collapse = ", ")
      gse.pathway <- gse.pathway %>% dplyr::select(-c("leadingEdge"))

      if (do.plot){
        # get top GSEA
        gsea.top <- (gse.pathway %>% dplyr::filter(NES > 0) %>% dplyr::arrange(log1p(pval)))[1:plot.top.n, ]
        gsea.bottom <- (gse.pathway %>% dplyr::filter(NES < 0) %>% dplyr::arrange(log1p(pval)))[1:plot.top.n, ]
        gse.pathway.top <- bind_rows(gsea.top, gsea.bottom)

        # plot
        gse.pathway.top$path.trun <- stringr::str_wrap(gse.pathway.top$pathway, path.name.wrap.width)
        plt.gsea <- gse.pathway.top %>%
          ggplot(aes(x = NES, y = reorder(path.trun, NES), fill = -log1p(pval), size = -log1p(pval))) +
          geom_segment(aes(x = 0, xend = NES, y = reorder(path.trun, NES), yend = reorder(path.trun, NES)), size = 0.05) +
          geom_point(pch=21) +
          theme_miko(legend = T) +
          xlab("NES") + ylab("") +
          labs(title = "GSEA", fill = "-log(p)", size =  "-log(p)") +
          geom_vline(xintercept = 0, linetype = "dashed") +
          viridis::scale_fill_viridis(option ="B")

        return(list(
          gse.pathway = gse.pathway,
          gse.pathway.top = gse.pathway.top,
          plt.gsea = plt.gsea
        ))
      } else {
        return(gse.pathway)
      }

    })
  })

}


#' Run hypergeometric gene enrichment analysis.
#'
#' Run hypergeometric gene enrichment analysis.
#'
#' @param gene.list Named list of genesets to run enrichment on (symbol format).
#' @param gene.universe Gene universe (symbol format)
#' @param species Species. One of "Mm" (mouse) or "Hs" (human)
#' @param pathway.db Name of databse to get annotation genelists from. Default is "Bader". See scMiko::getAnnotationPathways() for options.
#' @param n.workers Number of workers for parallelization. Default is 16.
#' @param my.pathway Named list of pathways, with each entry containing vector of Entrez IDs. Retrieved internally unless explicitly provided.
#' @param my.pathway.representation If my.pathway is provided, species which format, options = "SYMBOL" or "ENTREZ". ENTREZ is Default.
#' @param min.size min geneset size. Default is 2.
#' @param max.size max geneset size. Default is 300.
#' @param e2s entrez to symbol mapping (computationally demanding). Default False.
#' @param go.ontology BP, MF or CC. Ignored if pathway.db != "GO".
#' @param verbose Print progress. Default is TRUE.
#' @return enrichment results
#' @author Nicholas Mikolajewicz
#'
runHG <- function(gene.list, gene.universe,species, pathway.db = "Bader", n.workers = 16, my.pathway = NULL, my.pathway.representation = "ENTREZ", min.size = 2, max.size = 300, e2s = F, go.ontology = "BP", verbose = T){


  if (is.null(names(gene.list)) | ("" %in% names(gene.list)))  stop("One or more entries in the provided gene list is unnamed.")

  suppressMessages({
    suppressWarnings({

      require("foreach", quietly = T); require("parallel", quietly = T); require("fgsea", quietly = T); require("plyr", quietly = T); require("dplyr", quietly = T)

      my.symbol <- gene.universe
      my.entrez <- sym2entrez(my.symbol, my.species = species )
      my.entrez <- my.entrez[complete.cases(my.entrez), ]

    })})

  miko_message("Preparing pathway genesets...", verbose = verbose)
  if (is.null(my.pathway)){
    pathways <- getAnnotationPathways(query.genes = my.entrez$ENTREZID, db = pathway.db, ontology = go.ontology, species = species)
  } else {
    pathways <- my.pathway
  }

  pathway.size <- unlist(lapply(pathways, length))
  pathways <- pathways[(pathway.size >= min.size) & (pathway.size <= max.size)]
  names(pathways) <- make.unique(names(pathways))

  g2e.list <- my.entrez$ENTREZID
  names(g2e.list) <- my.entrez$SYMBOL
  e2g.list <- names(g2e.list)
  names(e2g.list) <- g2e.list

  gene.universe <- unique(gene.universe)
  gene.universe.original <- gene.universe
  gene.universe <- unique(g2e.list[gene.universe])

  n.av.core <- parallel::detectCores()
  if (n.workers > n.av.core) n.workers <- n.av.core
  if (n.workers > length(gene.list)) n.workers <- length(gene.list)
  cl <- parallel::makeCluster(n.workers)
  doParallel::registerDoParallel(cl)

  miko_message("Running hypergeometric enrichment...", verbose = verbose)
  res.h.list <- list()
  res.h.list <- foreach(i = 1:length(gene.list), .packages = c("dplyr", "fgsea", "plyr"))  %dopar% {

    if (my.pathway.representation == "ENTREZ"){
      if (!is.null(g2e.list)){
        current.genes <- unique(g2e.list[gene.list[[i]]])
      } else {
        current.genes <- unique(gene.list[[i]])
      }
    } else {
      gene.universe <- gene.universe.original
      current.genes <- unique(gene.list[[i]])
    }

    res.hyper <-  fgsea::fora(pathways = pathways, genes = current.genes, universe = gene.universe, minSize = 1, maxSize = Inf)

    if (e2s){
      suppressMessages({
        suppressWarnings({
          res.hyper$overlapGenes <- lapply(res.hyper$overlapGenes,
                                           mapvalues,from = my.entrez$ENTREZID, to = my.entrez$SYMBOL)
          res.hyper$overlapGenes <- lapply(res.hyper$overlapGenes, paste,collapse = ", ")
        })})
    }

    return(res.hyper)

  }

  parallel::stopCluster(cl)
  names(res.h.list) <- names(gene.list)

  miko_message("Complete!", verbose = verbose)

  return(res.h.list)

}



#' Convert gene ensembl to symbol
#'
#' Gene ensemble id is converted to symbol using org.Hs.eg.db or org.Mm.eg.db annotation databases.
#'
#' @param my.ensembl Character. Vector of ENSEMBL ids.
#' @param my.species Character. Species, one of "Mm" or "Hs".
#' @name sym2ens
#' @author Nicholas Mikolajewicz
#' @return data.frame mapping gene Ensemble to Symbol
#'
ensembl2sym <- function(my.ensembl, my.species){

  my.ensembl <- as.vector(my.ensembl)
  if (my.species == "Hs"){
    db <- org.Hs.eg.db::org.Hs.eg.db
  } else if (my.species == "Mm"){
    db <- org.Mm.eg.db::org.Mm.eg.db
  }
  my.symbol <- AnnotationDbi::select(db,
                                     keys = my.ensembl,
                                     columns = c("ENSEMBL", "SYMBOL"),
                                     keytype = "ENSEMBL")
  return(my.symbol)
}





#' Clean clusters
#'
#' Filter cells that exceed N median absolute deviations from cluster center in UMAP coordinate space. Assumes that Seurat::FindClusters() has already been run.
#'
#' @param object Seurat object.
#' @param mad.threshold Filter threshold. Represents minimal number of median absolute deviations from median required to omit cell from object. Default is 3.
#' @param return.plots Logical indicating whether to return plots. If true, returns list containing object along with ggplot handles. If false, only filtered object is returned.
#' @param verbose Logical indicating whether to print progress messages.
#' @name cleanCluster
#' @author Nicholas Mikolajewicz
#' @return list (return.plots == T) or seurat object (return.plots == F)
#'
cleanCluster <- function(object, mad.threshold = 3, return.plots = F, verbose = T){

  # assertion
  if (!("seurat_clusters" %in% colnames(object@meta.data))) stop("Cluster data not found. Run FindCluster() prior to Seurat::cleanClusters()")

  # get unique clusters
  u.cluster <- unique(as.character(object@meta.data[["seurat_clusters"]]))
  u.cluster <- u.cluster[order(as.numeric(u.cluster))]

  # get umap data
  df.umap <- data.frame(object@reductions[["umap"]]@cell.embeddings)
  colnames(df.umap) <- c("x", "y")
  df.umap$cluster <- object@meta.data[["seurat_clusters"]]
  df.umap$cells <- rownames(df.umap)

  # clean clusters...
  miko_message("Cleaning clusters...", verbose = verbose)
  if (return.plots) plot.list <- list()
  remove.which <- c()
  for (i in 1:length(u.cluster)){

    cluster.id <- u.cluster[i]
    df.umap.current <- df.umap[df.umap$cluster %in% cluster.id, ]

    # cluster centers
    df.center <- data.frame(
      median.x = median(df.umap.current$x, na.rm = T),
      median.y = median(df.umap.current$y, na.rm = T)
    )

    # calculate distances from center
    df.umap.current$dist <- sqrt(((df.umap.current$x - df.center$median.x)^2) + ((df.umap.current$y - df.center$median.y)^2))

    # filter cells
    median.dist = median(df.umap.current$dist, na.rm = T)
    mad.dist = mad(df.umap.current$dist, na.rm = T)
    df.umap.current$do.color <- "black"
    df.umap.current$do.color[df.umap.current$dist >  ( median.dist + (mad.dist * mad.threshold))] <- "tomato"
    current.cells.omitted <- df.umap.current$cells[df.umap.current$do.color == "tomato"]
    remove.which <- c(remove.which, current.cells.omitted)

    # get filtering statistics
    n.omitted <- sum(df.umap.current$do.color == "tomato")
    n.retained <- sum(df.umap.current$do.color == "black")
    p.omitted <- round(100*n.omitted/(n.omitted + n.retained), 3)

    if (return.plots){
      plt.umap <-  ggplot() +
        geom_point(data =df.umap.current, aes(x = x, y = y), color = df.umap.current$do.color) +
        geom_point(data = df.center, aes(x = median.x, y = median.y), color ="tomato", size = 4, shape = 4) +
        xlab("UMAP 1") + ylab("UMAP 2") +
        theme_miko() +
        labs(caption = "x = center; black = retain; red = omit", title = "",
             subtitle = paste0(p.omitted, "% cells omitted"))

      plt.hist <- df.umap.current %>%
        ggplot(aes(x = (dist))) +
        geom_histogram(bins = 30) +
        theme_miko() +
        geom_vline(xintercept = median.dist, color  = "tomato") +
        geom_vline(xintercept = median.dist + (mad.dist * mad.threshold), linetype = "dashed", color  = "tomato") +
        xlab("Distance from center (euclidean)") + ylab("Count") +
        labs(title = paste0("Cluster", cluster.id), subtitle = "Distance from center distribution", caption = c("solid line: median, dashed line = filter threshold"))

      plot.list[[paste0("c", cluster.id)]] <- cowplot::plot_grid(plt.hist, plt.umap)
    }
  }

  # final filter tally
  n.cells.all <- ncol(object)
  p.omitted.all <- round(100*length(remove.which) / n.cells.all, 3)

  # get before plot
  if (return.plots){
    original.plot <- cluster.UMAP(object) +
      labs(title = "Pre-filtering", subtitle = "") +
      theme_miko()
  }

  # subset object to return
  keep.which <- colnames(object)[!(colnames(object) %in% remove.which)]
  object <- object[ ,(colnames(object) %in% keep.which)]

  if (return.plots){
    cleaned.plot <- cluster.UMAP(object) +
      labs(title = "Post-filtering", subtitle = paste0(p.omitted.all, "% cells omitted")) +
      theme_miko()
    plot.all <- cowplot::plot_grid(original.plot, cleaned.plot, ncol = 2)
  }

  # return cleaned object
  if (return.plots){
    return(list(
      object = object,
      cluster.plots = plot.list,
      pre_post.plot = plot.all
    ))
  } else {
    return(object)
  }

}

#' Specify inputs for variance decomposition analysis
#'
#' Step 2 of variance decomposition analysis (see examples). Given Seurat object and vd_Formula output, input list for variance decomposition are generated.
#'
#' @param object Seurat object.
#' @param vd_model.list Output from vd_Formula.
#' @param Features Features to include in analysis. If specified, pct.min and variable features are ignored.
#' @param pct.min Minimal expressing fraction for genes to be included in analysis. Default is 0.
#' @param variable.features Logical specifying whether to use variable features only. If true, looks for variable features within provided Seurat object.
#' @param subsample.factor Numeric [0,1] specfying how to subsample (i.e., downsample) data. Default is 1 (no subsampling)
#' @name vd_Inputs
#' @concept vd
#' @author Nicholas Mikolajewicz
#' @return list of inputs for vd_Run() function.
#' @seealso \code{\link{vd_Run}}
#' @examples
#'
#'parameter.list <- list(
#'  covariates = c( "cluster", "percent.mt", "batch", "cycle", "seq.coverage"),
#'  interactions = c("batch:cluster")
#')
#'
#'# step 1: model formulation
#'vd_model.list <- vd_Formula(object = so.query,
#'                            covariates = parameter.list$covariates,
#'                            interactions = parameter.list$interactions)
#'
#'# step 2: prep model inputs
#'vd_inputs.list <- vd_Inputs(object = so.query, vd_model.list = vd_model.list, features = NULL,
#'                            pct.min =  0.9, variable.features = F, subsample.factor = 1)
#'
#'# step 3: run variance decomposition
#'vd_results.list <- vd_Run(vd_inputs.list, n.workers = 20)
#'
#'# step 4 (optional): visualize UMAP distribution of covariates
#'plt.umap.list <- vd_UMAP(object = so.query, vd_model.list = vd_model.list)
#'
#'# step 5 (optional): visualize decomposition
#'res.var2 <- vd_results.list$varPart.format1
#'plt.var <- plotVarPart( res.var2 ) +
#'  theme_miko() +
#'  labs(title = "Variance Decomposition", subtitle = "Linear Mixed-Effects Model")  +
#'  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
vd_Inputs <- function(object, vd_model.list, features = NULL, pct.min =  0, variable.features = F, subsample.factor = 1){

  require(lme4, quietly = T);

  # if features is specified, pct.min and variable features are ignored.

  # prep expression matrix
  e.mat <- object@assays[[DefaultAssay(object)]]@data

  if (is.null(features)){
    p.mat <- e.mat > 0
    p.exp <- BiocGenerics::rowSums(p.mat)/ncol(p.mat)
    which.gene <- which(p.exp > pct.min)
    e.mat.sub <- as.matrix(e.mat[which.gene, ])
  } else {
    if (sum(features %in% rownames(e.mat)) > 0){
      e.mat.sub <- as.matrix(e.mat[rownames(e.mat) %in% features, ])
    } else{
      stop("Specified features are not found in object. Cannot perform variance decomposition.")
    }
  }


  # use variable genes
  if ((length(object@assays[[DefaultAssay(object)]]@var.features) > 0) & (variable.features)& (is.null(features))){
    var.genes <- object@assays[[DefaultAssay(object)]]@var.features
    which.genes.match <- rownames(e.mat.sub) %in% var.genes
    if (sum(which.genes.match) > 0){
      e.mat.sub <- e.mat.sub[which.genes.match, ]
    }
  }

  # subsample matrix
  if ((subsample.factor < 1) & (subsample.factor > 0)){
    subsample.ind <- sample(seq(1, ncol(e.mat.sub)), round(subsample.factor * ncol(e.mat.sub)))
    e.mat.sub <- e.mat.sub[ ,subsample.ind]
    df.meta.sub <- df.meta.sub[subsample.ind ,]
  }

  fitVarPartModel2 <- function( exprObj, formula, data){

    # exprObj = e.mat.sub[ ,subsample.ind]
    # formula = form
    # data = df.meta.sub[ subsample.ind,]
    REML=FALSE
    useWeights=TRUE
    weightsMatrix=NULL
    showWarnings=TRUE
    fxn=identity
    colinearityCutoff=.999
    control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )

    # nested helper function
    exprIter = function( exprObj, weights, useWeights = TRUE, scale=TRUE, iterCount = "icount"){

      n_features = nrow(exprObj)

      if( iterCount == 'icount2'){
        xit <- icount2( n_features )
      }else{
        xit <- icount( n_features )
      }

      nextEl <- function() {
        j <- nextElem(xit)

        if( is.null(j) || j > n_features){
          res = NULL
        }else{
          if( useWeights && !is.null(weights) ){
            # scale weights to have mean of 1, otherwise it affects the residual variance too much
            if(scale){
              w = weights[j,] /  mean(weights[j,])
            }else{
              w = weights[j,]
            }
          }else{
            w = NULL
          }

          res = list(E = exprObj[j,], weights = w, n_iter = j, max_iter = n_features)
        }
        res
      }
      it <- list(nextElem = nextEl)
      class(it) <- c("abstractiter", "iter")
      it
    }

    formula = stats::as.formula( formula )

    # only retain columns used in the formula
    data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]

    # check dimensions of reponse and covariates
    if( ncol(exprObj) != nrow(data) ){
      stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
    }

    # check if all genes have variance
    if( ! is(exprObj, "sparseMatrix")){
      # check if values are NA
      countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
      if( countNA > 0 ){
        stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
      }

      rv = apply( exprObj, 1, var)
    }else{
      rv = c()
      for( i in seq_len(nrow(exprObj)) ){
        rv[i] = var( exprObj[i,])
      }
    }
    if( any( rv == 0) ){
      idx = which(rv == 0)
      stop(paste("Response variable", idx[1], 'has a variance of 0'))
    }

    # if weightsMatrix is not specified, set useWeights to FALSE
    if( useWeights && is.null(weightsMatrix) ){
      # warning("useWeights was ignored: no weightsMatrix was specified")
      useWeights = FALSE
    }

    # if useWeights, and (weights and expression are the same size)
    if( useWeights && !identical( dim(exprObj), dim(weightsMatrix)) ){
      stop( "exprObj and weightsMatrix must be the same dimensions" )
    }

    # If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
    if( ! identical(colnames(exprObj), rownames(data)) ){
      warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
    }

    form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

    responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
    possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,control=control, ... ), error = function(e) e)

    # detect error when variable in formula does not exist
    if( inherits(possibleError, "error") ){
      err = grep("object '.*' not found", possibleError$message)
      if( length(err) > 0 ){
        stop("Variable in formula is not found: ", gsub("object '(.*)' not found", "\\1", possibleError$message) )
      }
    }

    # fit first model to initialize other model fits - this make the other models converge faster
    responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
    # fitInit <- lmer( eval(parse(text=form)), data=data, REML=REML, control=control )

    fitInit <- tryCatch({
      lmer( eval(parse(text=form)), data=data, REML=REML, control=control )
    }, error = function(e){
      return(lm( eval(parse(text=form)), data=data ))
    }
    )

    # specify gene explicitly in data
    data2 = data.frame(data, expr=responsePlaceholder$E, check.names=FALSE)
    form = paste( "expr", paste(as.character( formula), collapse=''))

    input.data <- list(
      E = exprObj,
      weights = matrix(data = 1, nrow = nrow(exprObj), ncol = ncol(exprObj))
    )

    output.list <- tryCatch({
      list(
        input.data = input.data,
        data = data2,
        form =  form,
        REML =  REML,
        theta =  fitInit@theta,
        fxn = fxn,
        control = control,
        na.action=stats::na.exclude
      )}, error = function(e){
        return( list(
          input.data = input.data,
          data = data2,
          form =  form,
          REML =  REML,
          theta = NULL,
          fxn = fxn,
          control = control,
          na.action=stats::na.exclude
        ))
      })

    return(output.list)
  }

  # get model parameter list
  par.list <- fitVarPartModel2( exprObj = e.mat.sub, formula = vd_model.list$formula, data = vd_model.list$data)

  return(par.list)

}


#' Perform Variance Decomposition Analysis
#'
#' Step 3 of variance decomposition analysis (see examples).
#'
#' @param vd_inputs.list Output from scMiko::vd_Input() function.
#' @param n.workers Number of workers to use (for parallel implementation; uses foreach package)
#' @name vd_Run
#' @author Nicholas Mikolajewicz
#' @concept vd
#' @return List of results summarizing variance explained by each model covariate.
#' @seealso \code{\link{vd_Inputs}}
#' @examples
#'
#'parameter.list <- list(
#'  covariates = c( "cluster", "percent.mt", "batch", "cycle", "seq.coverage"),
#'  interactions = c("batch:cluster")
#')
#'
#'# step 1: model formulation
#'vd_model.list <- vd_Formula(object = so.query,
#'                            covariates = parameter.list$covariates,
#'                            interactions = parameter.list$interactions)
#'
#'# step 2: prep model inputs
#'vd_inputs.list <- vd_Inputs(object = so.query, vd_model.list = vd_model.list, features = NULL,
#'                            pct.min =  0.9, variable.features = F, subsample.factor = 1)
#'
#'# step 3: run variance decomposition
#'vd_results.list <- vd_Run(vd_inputs.list, n.workers = 20)
#'
#'# step 4 (optional): visualize UMAP distribution of covariates
#'plt.umap.list <- vd_UMAP(object = so.query, vd_model.list = vd_model.list)
#'
#'# step 5 (optional): visualize decomposition
#'res.var2 <- vd_results.list$varPart.format1
#'plt.var <- plotVarPart( res.var2 ) +
#'  theme_miko() +
#'  labs(title = "Variance Decomposition", subtitle = "Linear Mixed-Effects Model")  +
#'  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
vd_Run <- function(vd_inputs.list, n.workers = 20){

  require(foreach, quietly = T);
  require(lme4, quietly = T);

  # initiate clusters
  cl <- parallel::makeCluster(n.workers)
  doParallel::registerDoParallel(cl)

  # determine chunk size
  chunk.size <- ceiling(nrow(vd_inputs.list$input.data$E) / n.workers)
  fxn <- vd_inputs.list$fxn
  chunk.start <- seq(1, n.workers*chunk.size, by = chunk.size)

  res <- foreach(j = 1:(length(chunk.start)), .packages = c("lme4"))  %dopar% {

    if (j < length(chunk.start)){
      seq.range <- (chunk.start[j]:(chunk.start[j+1]-1))
    } else if (j == length(chunk.start)){
      seq.range <- (chunk.start[j]:(chunk.start[j] + chunk.size))
    }

    if (min(seq.range) > nrow(vd_inputs.list$input.data$weights)) return(list(res.cur = list(), which.gene = c()))

    if (max(seq.range) > nrow(vd_inputs.list$input.data$weights)) seq.range <- (chunk.start[j]: nrow(vd_inputs.list$input.data$weights))
    E.cur <- vd_inputs.list$input.data$E[seq.range,]
    W.cur <- vd_inputs.list$input.data$weights[seq.range,]
    data3 <- vd_inputs.list$data
    res.cur <- list()
    which.gene <- c()

    if (class(E.cur) == "numeric"){
      E.cur <- t(as.matrix(E.cur))
      rownames(E.cur) <- rownames(vd_inputs.list$input.data$E)[seq.range]
      W.cur <- t(as.matrix(W.cur))
      rownames(W.cur) <- rownames(vd_inputs.list$input.data$weights)[seq.range]
    }


    for (k in 1:nrow(E.cur) ){

      try({
        data3$expr <- E.cur[k,]
        # res.cur[[rownames(E.cur)[k]]] <- fxn( lmer( eval(parse(text=vd_inputs.list$form)), data=data3, REML=F,
        #                                             weights=W.cur[k,],
        #                                             control=vd_inputs.list$control,na.action= vd_inputs.list$na.action,
        #                                             start = list( theta = vd_inputs.list$theta)))

        res.cur[[rownames(E.cur)[k]]] <-   tryCatch({
          fxn( lmer( eval(parse(text=vd_inputs.list$form)), data=data3, REML=F,
                     weights=W.cur[k,],
                     control=vd_inputs.list$control,na.action= vd_inputs.list$na.action,
                     start = list( theta = vd_inputs.list$theta)))
        }, error = function(e){
          fxn( lm( eval(parse(text=vd_inputs.list$form)), data=data3,
                   weights=W.cur[k,],
          ))
        })

        which.gene <- c(which.gene, rownames(E.cur)[k])

      }, silent = T)

    }

    return(list(res.cur, which.gene))

  }

  parallel::stopCluster(cl)

  # unpack results
  res.var <-list()
  for (i in 1:length(res)){
    res.var <- c(res.var, res[[i]][[1]])
  }

  # helper function to extract results
  extractVarPart2 <- function( modelList, showWarnings=TRUE,... ){

    require("variancePartition", quietly = T)

    # get results from first model to enumerate all variables present
    # singleResult = calcVarPart( modelList[[1]], showWarnings=showWarnings,... )

    # for each model fit, get R^2 values
    entry <- 1
    varPart <- lapply( modelList, function( entry )
      calcVarPart( entry, showWarnings=showWarnings,... )
    )

    varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
    colnames(varPartMat) <- names(varPart[[1]])
    rownames(varPartMat) <- names(modelList)

    modelType = ifelse(class(modelList[[1]])[1] == "lm", "anova", "linear mixed model")

    new("varPartResults", varPartMat, type=modelType, method="Variance explained (%)")

  }

  # extract results
  res.var2 <- extractVarPart2(res.var)

  # convert to data.frame
  df.var <- as.data.frame(res.var2)
  df.var$gene <- rownames(df.var)

  return(list(
    varPart.format1 = res.var2,
    varPart.format2 = df.var)
  )

}


#' Specify model formula for variance decomposition.
#'
#' Step 1 of variance decomposition analysis (see examples).
#'
#' @param object Seurat Object
#' @param covariates vector of model covariates
#' @param interactions vector of interaction terms.
#' @name vd_Formula
#' @concept vd
#' @author Nicholas Mikolajewicz
#' @return List containing model formula and data.frame of covariates.
#' @seealso \code{\link{vd_Run}}
#' @examples
#'
#'parameter.list <- list(
#'  covariates = c( "cluster", "percent.mt", "batch", "cycle", "seq.coverage"),
#'  interactions = c("batch:cluster")
#')
#'
#'# step 1: model formulation
#'vd_model.list <- vd_Formula(object = so.query,
#'                            covariates = parameter.list$covariates,
#'                            interactions = parameter.list$interactions)
#'
#'# step 2: prep model inputs
#'vd_inputs.list <- vd_Inputs(object = so.query, vd_model.list = vd_model.list, features = NULL,
#'                            pct.min =  0.9, variable.features = F, subsample.factor = 1)
#'
#'# step 3: run variance decomposition
#'vd_results.list <- vd_Run(vd_inputs.list, n.workers = 20)
#'
#'# step 4 (optional): visualize UMAP distribution of covariates
#'plt.umap.list <- vd_UMAP(object = so.query, vd_model.list = vd_model.list)
#'
#'# step 5 (optional): visualize decomposition
#'res.var2 <- vd_results.list$varPart.format1
#'plt.var <- plotVarPart( res.var2 ) +
#'  theme_miko() +
#'  labs(title = "Variance Decomposition", subtitle = "Linear Mixed-Effects Model")  +
#'  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
vd_Formula <- function(object, covariates = NULL, interactions  = NULL){

  if(class(object) != "Seurat") stop("Seurat object must be provided as input.")
  df.meta <- object@meta.data

  # get model covariates
  if (!is.null(covariates)){
    covariate.names <- covariates
  } else {
    stop("Covariates specification is missing. ")
  }

  # get model interactions
  if (!is.null(interactions)){
    interaction.names <- interactions
  } else {
    interaction.names <- c()
  }


  # check which.available
  covariates.av <- covariate.names[covariate.names %in% colnames(df.meta)]

  # check if interactions are available
  df.interaction.pairs <- NULL
  if (length(interaction.names) > 0){
    for (i in 1:length(interaction.names)){
      cov.int <- c(unlist(strsplit(interaction.names[i], ":")))
      if(all(cov.int %in% covariates.av)){
        df.interaction.pairs <- bind_rows(
          data.frame(v1 = cov.int[1], v2 = cov.int[2])
        )
      }
    }
  }

  if (!is.null(df.interaction.pairs)){
    do.interaction <- T
  } else {
    do.interaction <- F
  }


  append.term <- function(func, term){
    if (length(func) == 0){
      func <- paste0("~ ",term)
    } else {
      func <- paste0(func, " + ", term)
    }
    return(func)
  }

  # variable types
  df.var.type <- NULL
  form <- c()
  for (i in 1:length(covariates.av)){

    var.current <-  df.meta[ ,covariates.av[i]]
    if ((is.numeric(var.current)) | (is.integer(var.current))){
      df.var.type <- bind_rows(df.var.type,
                               data.frame(
                                 variable = covariates.av[i],
                                 type = "continuous"
                               ))

      form <- append.term(form, covariates.av[i])
    } else if ((is.character(var.current)) | (is.factor(var.current))){
      df.var.type <- bind_rows(df.var.type,
                               data.frame(
                                 variable = covariates.av[i],
                                 type = "categorical"
                               ))
      form <- append.term(form, paste0("(1|", covariates.av[i], ")"))
    }


  }

  # append interaction terms
  if (!is.null(df.interaction.pairs)){
    for (i in 1:nrow(df.interaction.pairs)){

      term.current <-  paste0(df.interaction.pairs[i,1], ":", df.interaction.pairs[i,2])
      form <- append.term(form,  paste0("(1|", term.current, ")"))

    }
  }

  # specify model formula
  form2 <- as.formula(form)

  # get covariate data
  df.meta <- object@meta.data
  df.meta.sub <- df.meta[ ,covariates.av]

  if (length(covariates.av) == 1){
    if (!(("data.frame") %in% class(df.meta.sub))){
      df.meta.sub <- as.data.frame(df.meta.sub)
      rownames(df.meta.sub) <- rownames(df.meta)
      colnames(df.meta.sub) <- covariates.av
    }
  }


  return(
    list(
      formula = form2,
      data = df.meta.sub,
      variable.type = df.var.type
    )
  )

}


#' Generate UMAPs with each variance decomposition covariate overlaid.
#'
#' Generate UMAPs with each variance decomposition covariate overlaid.
#'
#' @param object Seurat Object
#' @param vd_model.list Output from vd_Formula().
#' @name vd_UMAP
#' @concept vd
#' @author Nicholas Mikolajewicz
#' @return List ggplot handles, one for each model covariate.
#' @seealso \code{\link{vd_Run}}
#' @examples
#'
#'parameter.list <- list(
#'  covariates = c( "cluster", "percent.mt", "batch", "cycle", "seq.coverage"),
#'  interactions = c("batch:cluster")
#')
#'
#'# step 1: model formulation
#'vd_model.list <- vd_Formula(object = so.query,
#'                            covariates = parameter.list$covariates,
#'                            interactions = parameter.list$interactions)
#'
#'# step 2: prep model inputs
#'vd_inputs.list <- vd_Inputs(object = so.query, vd_model.list = vd_model.list, features = NULL,
#'                            pct.min =  0.9, variable.features = F, subsample.factor = 1)
#'
#'# step 3: run variance decomposition
#'vd_results.list <- vd_Run(vd_inputs.list, n.workers = 20)
#'
#'# step 4 (optional): visualize UMAP distribution of covariates
#'plt.umap.list <- vd_UMAP(object = so.query, vd_model.list = vd_model.list)
#'
#'# step 5 (optional): visualize decomposition
#'res.var2 <- vd_results.list$varPart.format1
#'plt.var <- plotVarPart( res.var2 ) +
#'  theme_miko() +
#'  labs(title = "Variance Decomposition", subtitle = "Linear Mixed-Effects Model")  +
#'  theme(axis.text.x = element_text(angle = 45, hjust = 1))
#'
vd_UMAP <- function(object, vd_model.list){

  df.var.type <- vd_model.list$variable.type
  plt.umap.list <- list()
  for (i in 1:nrow(df.var.type)){


    field.name <-df.var.type$variable[i]
    var.type <- df.var.type$type[i]

    if (var.type == "categorical"){
      plt.umap.list[[field.name]] <-  cluster.UMAP(object, group.by = field.name) + theme_miko(legend = T) +
        labs(title = field.name, subtitle = "UMAP")
    } else if (var.type == "continuous"){
      plt.umap.list[[field.name]] <-  FeaturePlot(object, feature = field.name,
                                                  cols = c("lightgrey", "darkgreen")) +
        theme_miko(legend = T) +
        labs(title = field.name, subtitle = "UMAP") +
        xlab("UMAP 1") + ylab("UMAP 2")
    }

  }

  return(plt.umap.list)

}



#' Assign column entries in data.frame to row names.
#'
#' Assign column entries in data.frame to row names.
#'
#' @param df data.frame
#' @param col Character specifying column name in df.
#' @name col2rowname
#' @author Nicholas Mikolajewicz
#' @return df
#' @examples
#'
#' f.mat <-aggGroupExpression(
#' so = so.query,
#'   which.data = "data",
#'   which.assay = DefaultAssay(so.query),
#'   which.center = "fraction",
#'   which.group = "seurat_clusters",
#'   do.parallel = F
#' )
#'
#' f.mat <- col2rowname(f.mat, col = "genes")
#'
col2rowname <- function(df, col){

  df <- as.data.frame(df)
  if (!(col) %in% colnames(df)) {
    warning("The provided col argument is a not a valid column name. Returning unmodified data.frame.")
    return(df)
  } else {
    df[,col] <- as.character(df[,col])
    rownames(df) <- make.unique(df[,col])
    df <- df %>% dplyr::select(-c(all_of(col)))
    return(df)
  }

}


#' Get UMAP data and plot from Seurat object.
#'
#' Get UMAP data and plot from Seurat object.
#'
#' @param object Seurat Object.
#' @param umap.key Character UMAP key slot in seurat object. Default is "umap"
#' @param node.type "point" or "text"
#' @param meta.features character vector specifying which meta features to retrieve. Default is seurat_clusters
#' @param ... additional parameters to geom_point or geom_text
#' @name getUMAP
#' @author Nicholas Mikolajewicz
#' @return list containing UMAP data.frame and ggplot handle
#' @examples
#'
#' wnnUMAP.list <- getUMAP(so.gene, umap.key = "wnn.umap", node.type = "text")
#'
#' df.wnn.umap <- wnnUMAP.list$df.umap
#' plt.wnn.umap <- wnnUMAP.list$plt.umap
#'
getUMAP <- function(object, umap.key = "umap", node.type = "point", meta.features = "seurat_clusters", size = autoPointSize(ncol(object)), ...){

  # node.type "text" or "point"

  df.umap <- data.frame(object@reductions[[umap.key]]@cell.embeddings)
  colnames(df.umap) <- c("x", "y")
  df.umap$var <- rownames(df.umap)

  for (i in 1:length(meta.features)){
    if (meta.features[i] %in% colnames(object@meta.data)){
      df.umap[ ,meta.features[i]]<- object@meta.data[[meta.features[i]]]
    }
  }

  if ("seurat_clusters" %in% meta.features){
    plt.umap = df.umap %>%
      ggplot(aes(x = x, y = y, label = var, color = seurat_clusters)) +
      theme_miko(legend = T)
  } else {
    plt.umap = df.umap %>%
      ggplot(aes(x = x, y = y, label = var)) +
      theme_miko(legend = F)
  }


  if (node.type == "text"){
    plt.umap <- plt.umap +  geom_text(size = size, ...)
  } else if (node.type == "point"){
    plt.umap <- plt.umap +  geom_point(size = size, ...)
  }

  return(
    list(
      df.umap = df.umap,
      plt.umap = plt.umap
    )
  )
}


#' Generate categorical ColorBrewer palette.
#'
#' Generate categorical ColorBrewer palette. Specify one of 'labels' or 'n' to generate categorical color palette for use with 'scale_fill_manual(values = col.pal)' and 'scale_color_manual(values = col.pal)'
#'
#' @param labels vector of category names that if specified will have color assigned to. Default is NULL.
#' @param n number of colors. Default is NULL. Ignored if 'labels' is specified.
#' @param palette Name of palette. Default: "Spectral"
#' @name categoricalColPal
#' @author Nicholas Mikolajewicz
#' @return Vector of colors. Named if labels provided.
#' @examples
#'
#' # generate color palette
#' col.pal <- categoricalColPal(n = 10)
#'
#' # add to ggplot handle
#' gg.plt <- gg.plt + scale_fill_manual(values = col.pal)
#'
categoricalColPal <- function(labels = NULL, n = NULL, palette = "Spectral"){

  require("RColorBrewer", quietly = T)

  if (is.null(labels) & is.null(n)){
    miko_message("'Error: 'labels' or 'n' were not specified. Color palette was not generated.")
    return(NULL)
  } else if (!is.null(labels)){
    n <- length(labels)
  }

  all.pal <- brewer.pal.info

  if (palette %in% rownames(all.pal)){
    col.pal <- colorRampPalette(RColorBrewer::brewer.pal(all.pal[palette, "maxcolors"], palette))(n)
  } else {
    miko_message("Error: Specified palette must belong to RColorBrewer. See brewer.pal.info for available palettes.")
    return(NULL)
  }

  if (!is.null(labels)){
    names(col.pal) <- labels
  }

  return(col.pal)

}



#' Identify expressed genes
#'
#' Identify expressed genes in Seurat object.
#'
#' @param object Seurat Object
#' @param min.pct minimum expressing fraction. Default: 0.1
#' @param group Character specifying metadata field to group cells by. If not specified, global expression fraction is evaluated. If specified, group-level gene lists are combined used group.boolean.
#' @param group.boolean Boolean used to combine group genelists. One of "OR" or "AND". Default: "OR". Argument is ignored if 'group' is not specified.
#' @name getExpressedGenes
#' @author Nicholas Mikolajewicz
#' @return vector of gene names
#' @examples
#'
#' split.var <- "seurat_clusters"
#' which.genes <- getExpressedGenes(object = so.query, min.pct = 0.1, group = split.var, group.boolean = "OR")
#'
getExpressedGenes <- function(object, min.pct = 0.1, group = NA, group.boolean = "OR"){

  # min.groups: specify min.pct fraction grouping
  # NA: min.pct satisfied global
  # character vector speifying metadata feature to group: min.pct satisfied within object meta.data grouping

  # min.type.operator
  # AND: must satisfy criteria in ALL subgroup
  # OR: must satisfy critreria in ATLEAST ONE group

  emat <- getExpressionMatrix(object, which.data = "data")

  if (is.na(group)){

    # pct.rep <- apply(emat, 1, function(x) mean(x>0))
    pct.rep <-Matrix::rowMeans(emat>0) # faster implementation

    expressed.genes <- names(pct.rep)[pct.rep > min.pct]
  } else if (group %in% colnames(object@meta.data)){
    group.var <- as.character(object@meta.data[ ,group])
    u.gv <- unique(group.var)
    expressed.genes.all <- c()
    for (i in 1:length(u.gv)){
      pct.rep <- NULL
      if (!is.null(dim(emat[ ,group.var %in% u.gv[i]]))){
        pct.rep <- Matrix::rowMeans(emat[ ,group.var %in% u.gv[i]]>0)
        expressed.genes.all <- c(expressed.genes.all, names(pct.rep)[pct.rep > min.pct])
      }
    }

    df.exp.gene <- data.frame(table(expressed.genes.all))

    if (group.boolean == "OR"){
      expressed.genes <- unique(as.character(df.exp.gene$expressed.genes.all))
    } else if (group.boolean == "AND"){
      expressed.genes <- unique(as.character(df.exp.gene$expressed.genes.all[df.exp.gene$Freq == max(df.exp.gene$Freq, na.rm = T)]))
    } else {
      miko_message(paste0("'", group.boolean, "' is not an accepted argument for group.boolean. Must be either 'OR' or 'AND'. 'OR' was used as default argument"))
      expressed.genes <- unique(as.character(df.exp.gene$expressed.genes.all))
    }

  } else {
    miko_message(paste0( "'", group,  "' was not found in seurat object. Failed to identify expressed genes. "))
  }

  return(expressed.genes)
}




#' Subsample cells in seurat object to be balanced (sample-size-wise) across conditions.
#'
#' Subsample cells in seurat object to be balanced (sample-size-wise) across conditions.
#'
#' @param object Seurat Object
#' @param group Character specifying group (i.e., metadata field) to sub-sample within.
#' @param balance.size Numeric. Target sub sample size for each condition. If target exceeds number of available cells, only available cells are included. If unspecified, smallest-sized group is set as target (Default).
#' @name balanceSamples
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#' @examples
#'
#' so.query <- balanceSamples(object = so.query, group = "Barcode")
#'
balanceSamples <- function(object, group, balance.size = NA){

  # balances samples based on min sample present - otherwise if balance siz eis specified, matches to specification.

  stopifnot(class(object) == "Seurat")

  # get sample sizes
  df.meta <- object@meta.data
  stopifnot(group %in% colnames(df.meta))

  u.group <- unique(as.character(df.meta[,group]))

  if (is.na(balance.size)){
    df.tally <- data.frame(table(as.character(df.meta[,group])))
    colnames(df.tally) <- c("group", "n")
    target.n <- min(df.tally$n)
  } else {
    target.n <- balance.size
  }

  # get sampled indices
  all.select.cells <- c()
  for (i in 1:length(u.group)){
    av.cells <- rownames(df.meta)[df.meta[,group] %in% u.group[i]]
    n.cells <- length(av.cells)
    if (n.cells < target.n){
      cur.target <- n.cells
    } else {
      cur.target <- target.n
    }
    select.cells <- sample(x = av.cells, size = cur.target, replace = F)
    all.select.cells <- c(all.select.cells, select.cells)
    miko_message(paste0(signif(100*length(select.cells)/n.cells, 3), "% cells (",length(select.cells) , "/" ,n.cells, ") sampled from ", u.group[i], " group"))

  }

  # subsample
  object <- object[, colnames(object) %in% all.select.cells]

  return(object)

}


#' Project dimensionally-reduced features onto UMAP
#'
#' Project dimensionally-reduced features (e.g., PC, IC) onto UMAP
#'
#' @param object Seurat Object
#' @param reduction Character specifying reduction key. Default: "pca"
#' @param n.components Numeric. Number of top components to project.
#' @param show.n.features Numeric. Show top n features alongside projection. Default is 50.
#' @param pca.min.var.exp Threshold for minimum variance explained by principal component.
#' @param pca.cum.pca.thresh Threshold for maximum cumulative variance explained by principal component.
#' @param umap.reduction Character specifying reduction key for umap coordinates.
#' @param rel.width Relative widths of umap and top feature plots. Default: c(2.5,1)
#' @param ... additional parameters passed to geom_point(). Generally used to specify point size.
#' @name projectReduction
#' @author Nicholas Mikolajewicz
#' @return list of ggplot handles
#' @examples
#'
#'
projectReduction <- function(object, reduction = "pca", n.components = NA, show.n.features = 50, pca.min.var.exp = 0.05, pca.cum.pca.thresh = 0.8,   umap.reduction = "umap", rel.width  = c(2.5, 1), size = autoPointSize(ncol(object)),   ...){

  # so.query <- RunPCA(so.query, verbose = FALSE)

  # specify number of PCA components to use for downstream analysis

  if (grepl("pc", reduction)){
    pca.var.threshold <- pca.cum.pca.thresh
    pca.components <- propVarPCA(object, reduction = reduction)
    red.name <- "PC"
  }  else if (grepl("ic", reduction)){
    red.name <- "IC"
  } else {
    stop("invalid reduction")
  }

  if (is.na(n.components) & grepl("ic", reduction)){
    n.components <- 30
  } else if (!is.na(n.components) & grepl("pc", reduction)){
    # n.components <- 30
    pca.min.var.exp <- 0
    # n.pca <- max(pca.components$pc.id[pca.components$pc.cum_sum<pca.var.threshold])+1
  } else if (is.na(n.components) & grepl("pc", reduction)) {
    n.pca <- max(pca.components$pc.id[pca.components$pc.cum_sum<pca.var.threshold])+1
    n.components <- n.pca
  }




  df.pca.umap <- data.frame(
    x = object@reductions[[umap.reduction]]@cell.embeddings[, 1],
    y = object@reductions[[umap.reduction]]@cell.embeddings[, 2]
  )

  plt.list <- list()
  for (i in 1:n.components){

    if ( grepl("pc", reduction)){
      if (pca.components$pc.prop_var[i] < pca.min.var.exp) next
      var.exp <- paste0("(", signif(pca.components$pc.prop_var[i]* 100, 3), "% variance)")
    } else {
      var.exp <- ""
    }
    df.pca.umap$pc <- object@reductions[[reduction]]@cell.embeddings[ ,i]

    df.factor <- data.frame(
      feature = rownames(object@reductions[[reduction]]@feature.loadings),
      loading = (object@reductions[[reduction]]@feature.loadings[,i])
    )

    df.factor <- (df.factor %>% dplyr::arrange(-abs(loading)))[1: show.n.features, ]
    # dplyr::top_n(show.n.features, abs(loading))

    plt.pc.activity <- df.pca.umap %>%
      ggplot(aes(x = x, y = y, color = pc, size = abs(pc))) +
      geom_point(size = size, ...) +
      scale_color_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
      theme_miko(legend = T) +
      labs(x = "UMAP 1", y = "UMAP 2", title = paste0(red.name, " ", i, " ", var.exp), color = reduction)


    plt.pc.loading<- df.factor %>%
      ggplot(aes(y = reorder(feature, loading), x = loading, fill = loading)) +
      scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
      geom_bar(stat = "identity") +
      labs(x = "Factor Loading", y = "Feature", title = "Factor Loading",  fill = "Loading") +
      theme_miko(legend = T) +
      geom_vline(xintercept = 0, linetype = "dashed") + theme(axis.text.y = element_text(size = 8))

    plt.pc.merge <- cowplot::plot_grid(plt.pc.activity, plt.pc.loading, ncol = 2,rel_widths = rel.width)


    # print(plt.splice.pc.merge)

    plt.list[[paste0(red.name,i)]] <- plt.pc.merge

  }


  return(plt.list)


}



#' Get top loaded features for PCA or ICA dimensional reduction.
#'
#' Get top loaded features for PCA or ICA dimensional reduction.
#'
#' @param feature.loading Feature loading matrix
#' @param sig.threshold Numeric specifying top nth percentile of features to return. Default: 0.99.
#' @param top.is.positive Logical specifying whether to enforce top loaded feature being positive. If negative, loading sign is inverted.
#' @name getReductionGenes
#' @author Nicholas Mikolajewicz
#' @return list of ggplot handles
#' @examples
#'
#' feature.loading <- as.matrix(seurat.object@reductions[[reduction]]@feature.loadings)
#' red.res <- getReductionGenes(feature.loading = feature.loading,  sig.threshold = 0.99)
#'
getReductionGenes <- function(feature.loading, sig.threshold = 0.99, top.is.positive = T){

  if (!("matrix" %in% class(feature.loading))) stop("feature.loading input is not a matrix")

  ica.kme <- t(feature.loading)
  ica.sign <- apply(ica.kme, 1, function(x) sign(x[which.max(abs(x))]))

  if (top.is.positive){

    for (i in 1:nrow(ica.kme)){
      ica.kme[i,] <- ica.kme[i,] * ica.sign[i]
    }
  }

  module.genes.pos <-  apply(ica.kme, 1, function(x) colnames(ica.kme)[((x> quantile( x, sig.threshold)) & (x > 0))])
  module.genes.neg <-  apply(ica.kme, 1, function(x) colnames(ica.kme)[((x < quantile( x, 1-sig.threshold)) & (x < 0))])

  module.genes.pos.list <- list()
  module.genes.neg.list <- list()
  for (i in 1:ncol(module.genes.pos)){
    module.genes.pos.list[[colnames(module.genes.pos)[i]]] <- module.genes.pos[ ,i]
  }
  for (i in 1:ncol(module.genes.neg)){
    module.genes.neg.list[[colnames(module.genes.neg)[i]]] <- module.genes.neg[ ,i]
  }

  return(list(
    feature.lodaing = ica.kme,
    sign.correction = ica.sign,
    module.genes.pos = module.genes.pos.list,
    module.genes.neg = module.genes.neg.list))

}



#' Get unique features from metadata column in seurat object.
#'
#' Get unique features from metadata column in seurat object.
#'
#' @param object Seurat object
#' @param group name of metadata column to get features for.
#' @param is.numeric Logical to assert numerical class and order features.
#' @name uniqueFeatures
#' @author Nicholas Mikolajewicz
#' @return Seurat object
#' @examples
#'
#' so.query <- uniqueFeatures(object = so.query, group = "Barcode, is.numeric = F)
#'
uniqueFeatures <- function(object, group, is.numeric = F){
  # object: Seurat object
  # group: name of metadata column to get features for.

  df.meta <- object@meta.data

  if (!(group %in% colnames(df.meta))){
    miko_message(paste0( "'",group, "' is not a valid metadata column. No unique features returned. "))
    return(NULL)
  } else {
    all.bc <- df.meta[ ,group]
    u.bc <- unique(as.vector(all.bc))
    if (is.numeric){
      u.bc <- as.numeric(u.bc)
      u.bc <- order(u.bc)
    }

  }

  return(u.bc)

}


#' Create pseudo-replicates, stratified by grouping variable.
#'
#' Create pseudo-replicates, stratified by grouping variable.
#'
#' @param object Seurat object
#' @param split.by name of metadata column to create pseudo-replicates for.
#' @param n numerical, number of replicates to create.
#' @name pseudoReplicates
#' @author Nicholas Mikolajewicz
#' @return Returns a Seurat object where latest pseudo replicates results will be stored in object metadata under 'pseudo_replicates'. Note that 'pseudo_replicates' will be overwritten every time pseudoReplicates is run
#' @examples
#'
#' so.query <- uniqueFeatures(object = so.query, group = "Barcode, is.numeric = F)
#'
pseudoReplicates <- function(object, split.by, n = 2){

  # n: number of pseudosamples
  # split.within: name of metadata column to split samples within
  df.meta <- object@meta.data
  u.bc <- uniqueFeatures(object = object, group = split.by, is.numeric = F)
  pseudosample.name <- paste0(split.by, "_", n, "_split")
  df.meta$pseudo_replicates <- NA
  df.meta[ ,pseudosample.name] <- NA
  all.bc <- df.meta[ ,split.by]
  for (i in 1:length(u.bc)){
    which.bc <- which(all.bc %in% u.bc[i])
    n.bc <- sum(all.bc %in% u.bc[i])
    n.per.group <- round(n.bc/n)
    for (j in 1:n){

      if (j < n){
        g1 <- sample(which.bc, n.per.group)
        df.meta$pseudo_replicates[g1] <- paste0(u.bc[i], "_", j)
        which.bc <- which.bc[!(which.bc %in% g1)]
      } else {
        # g1 <- sample(which.bc, n.per.group)
        df.meta$pseudo_replicates[which.bc] <- paste0(u.bc[i], "_", j)
      }
    }

    df.meta[ ,pseudosample.name] <- df.meta$pseudo_replicates
  }
  object@meta.data <- df.meta
  return(object)
}





#' Number of unique values
#'
#' Number of unique values.
#'
#' @param x object
#' @name ulength
#' @author Nicholas Mikolajewicz
#' @return numeric
#' @examples
#'
#' x <- c(1,1,1,1,2,2,2,3,)
#' n.x <- ulength(x)
#'
ulength <- function(x){

  length(unique(x))
}



#' Compute purity of each cell's neighborhood, as defined by KNN graph.
#'
#' Compute purity of each cell's neighborhood, as defined by KNN graph. For each neighborhood, the proportion of cells belonging to the most represented cluster is computed as the purity score. The Purity score is stored in the metadata of the returned seurat object.
#'
#' @param object Seurat object
#' @param graph name of KNN graph to use for analysis. Must be present within provided seurat object. If absent, run FindNeighbors().
#' @param cluster.field name of cluster metadata field to use for cluster membership data. Default is 'seurat_clusters'.
#' @param verbose Print progress. Default is TRUE.
#' @name neighborPurity
#' @author Nicholas Mikolajewicz
#' @return Seurat object with purity score stored in metadata
#' @examples
#'
#' object <- FindNeighbors(object)
#' object <- FindClusters(object, resolution = 1)
#' object <- neighborPurity(object, "RNA_nn")
#'
neighborPurity <- function(object, graph, cluster.field = "seurat_clusters", verbose = T){

  if (graph %in% names(object@graphs)){
    nn.graph <- (object@graphs[[graph]])
    miko_message("Getting nearest neighbors...", verbose = verbose)
    nn.graph.ind <- apply(nn.graph, 1, function(x) which(x>0))
    # nn.graph.ind <- apply((nn.graph > 0), 1, which) # faster
    miko_message("Mapping nearest neighbors to cluster memberships...", verbose = verbose)
    cluster.vec <- as.vector(object@meta.data[ ,cluster.field])

    if ("matrix" %in% class(nn.graph.ind) ){
      nn.graph.clustID <- apply(nn.graph.ind, 2, function(x) cluster.vec[x])
      miko_message("Computing neighborhood purity...", verbose = verbose)
      nn.graph.clust.tally <- apply(nn.graph.clustID, 2, function(x) table(x))
      tally.size <- unlist(lapply(nn.graph.clust.tally, function(x) max(x/sum(x))))
      object@meta.data$purity <-tally.size
    } else if ("list" %in% class(nn.graph.ind)){
      nn.graph.clustID <- lapply(nn.graph.ind,  function(x) cluster.vec[x])
      miko_message("Computing neighborhood purity...", verbose = verbose)
      nn.graph.clust.tally <- lapply(nn.graph.clustID, function(x) table(x))
      tally.size <- unlist(lapply(nn.graph.clust.tally, function(x) max(x/sum(x))))
      object@meta.data$purity <-tally.size

    }
    miko_message("Complete!", verbose = verbose)

  } else {
    miko_message(paste0("'", graph, "' not found. Returning unmodified object."), verbose = verbose)
  }
  return(object)

}


#' Determine species based on gene representation
#'
#' Determine species (Hs, Mm) based on gene (symbol, ensembl) representation
#'
#' @param object Seurat object, gene expression matrix (rownames are genes), or gene character vector.
#' @name detectSpecies
#' @author Nicholas Mikolajewicz
#' @return Species (Hs or Mm)
#' @examples
#'
#' my.species <- detectSpecies(object)
#'
detectSpecies <- function(object){

  if (sum(c("Seurat", "Matrix", "dgCMatrix") %in% class(object)) > 0){
    my.rep <- (rownames(object))
  } else if ( "character" %in% class(object)) {
    my.rep <- object
  } else {
    miko_message("Input does not belong to Seurat, Matrix, dgCMatrix or character class. Unable to detect species.")
    return(NULL)
  }


  ens.sum <-  sum(grepl("ENS", my.rep))
  ens.mus.sum <-  sum(grepl("ENSMUS", my.rep))
  hi.cap.sum <-  sum(my.rep == toupper(my.rep))
  lo.cap.sum <-  sum(my.rep == firstup(my.rep))

  df.rep <-as.data.frame(t(data.frame(
    ens.sum = ens.sum,
    ens.mus.sum = ens.mus.sum,
    hi.cap.sum = hi.cap.sum,
    lo.cap.sum = lo.cap.sum
  )))

  which.rep <- rownames(df.rep)[which.max(df.rep[, 1])]

  if ((which.rep %in% c("ens.sum", "ens.mus.sum")) & (ens.sum > ens.mus.sum)) {
    species <- "Hs"
  } else if ((which.rep %in% c("ens.sum", "ens.mus.sum")) & (ens.sum <= ens.mus.sum)) {
    species <- "Mm"
  } else if (which.rep == "hi.cap.sum") {
    species <- "Hs"
  } else if (which.rep == "lo.cap.sum") {
    species <- "Mm"
  }

  return(species)
}

#' Print message
#'
#' Print message
#'
#' @param ... Character string to print
#' @param time Logical to print time
#' @param verbose Logical to print message
#' @name miko_message
#' @author Nicholas Mikolajewicz
miko_message <- function(..., time = T, verbose = T){
  if (verbose){
    if (time){
      message(Sys.time() , ": ",  ...)
    } else {
      message(...)
    }
  }
}


#' Winsorize values at lower and upper quantiles.
#'
#' Winsorize values at lower and upper quantiles.
#'
#' @param x Numeric vector
#' @param lower.quantile lower quantile [0, 1]. Default = 0.01
#' @param upper.quantile lower quantile [0, 1]. Default = 0.99
#' @name snip
#' @author Nicholas Mikolajewicz
snip <- function(x, lower.quantile = 0.01, upper.quantile = 0.99){
  lq <- quantile(x, lower.quantile, na.rm = T)
  uq <- quantile(x, upper.quantile, na.rm = T)

  x[x < lq] <- lq
  x[x > uq] <- uq

  return(x)
}

#' Returns intersection of all list entries.
#'
#' Returns intersection of all list entries.
#'
#' @param x list
#' @name lintersect
#' @author Nicholas Mikolajewicz
lintersect <- function(x){
  Reduce(intersect,  x)
}



#' Function to draw ggplot heatmaps
#'
#' Wrapper for pheatmap::pheatmap() and ggplotify::as.ggplot.
#'
#' @param mat numeric matrix
#' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. Corresponding values are "row", "column" and "none"
#' @param symmetric_scale Enforce symmetrical color scale. Default is TRUE.
#' @param scale.lim Apply ceiling and floor to all values in matrix. Default is NA.
#' @param color vector of colors used in heatmap
#' @param ... additional parameters passed to pheatmap::pheatmap(...)
#' @name miko_heatmap
#' @author Nicholas Mikolajewicz
miko_heatmap <- function(mat, scale = "none", symmetric_scale = T, scale.lim = NA, color = colorRampPalette(rev(brewer.pal(n = 7, name =  "RdBu")))(100), ...){


  require(RColorBrewer)

  if (symmetric_scale & (scale == "none")){
    if (is.na(scale.lim)){
      scale.lim <- max(abs(mat))
    }
    my.breaks <- seq(-scale.lim, scale.lim, by = (2*scale.lim)/length(color))
    plt.hm <- ggplotify::as.ggplot(pheatmap::pheatmap(mat, breaks = my.breaks, color = color, silent = T, ...) )
  } else {

    # remove zero variance col/row entries
    zero.var.row <-   apply(mat, 1, function(x) var(x) == 0)
    zero.var.col<-   apply(mat, 2, function(x) var(x) == 0)
    mat <- mat[!zero.var.row, !zero.var.col]

    plt.hm <- ggplotify::as.ggplot(pheatmap::pheatmap(mat, color = color, scale = scale, silent = T, ...) )

    if (scale != "none"){
      plt.hm <- plt.hm + labs(caption = paste0(scale, " scaled"))
    }
  }

  plt.hm <- plt.hm +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.subtitle = element_text(hjust = 0.5))

  return(plt.hm)

}



#' Annotate glioblastoma (GBM) subtype  based on Neftel 2019 scoring pipeline.
#'
#' Annotate glioblastoma (GBM) subtype  based on Neftel 2019 scoring pipeline.
#'
#' @param object seurat object
#' @param species Species, either "Mm" or "Hs".
#' @param verbose Print progress. Default is TRUE.
#' @param do.snip truncate scores at 1st and 99th percentiles (protects against outliers). Default is T.
#' @name scoreGBM
#' @author Nicholas Mikolajewicz
scoreGBM <- function(object, species = detectSpecies(object), verbose = T, do.snip = T){

  # get GBM genes
  miko_message("Getting GBM genesets...", verbose = verbose)
  gbm.genes <- data.frame(geneSets[["GBM_Hs_Neftel2019"]])

  # convert to list
  gbm.list <- list()
  for (i in 1:ncol(gbm.genes)){

    # enforce correct species
    if (species == "Hs"){
      gbm.genes[ ,i] <- toupper(gbm.genes[,i])
    } else if (species == "Mm"){
      gbm.genes[ ,i] <- firstup(gbm.genes[,i])
    }

    # get module name
    module.name <- colnames(gbm.genes[i])

    # assign to list
    gbm.list[[module.name]] <- gbm.genes[!is.na(gbm.genes[ ,i]) ,i]
  }

  module.scores <- matrix(ncol = length(gbm.list), nrow = ncol(object))
  miko_message("Scoring cells...", verbose = verbose)
  suppressWarnings({
    suppressMessages({


  for (i in 1:length(gbm.list)){
    current.list <- list(g1 = gbm.list[[i]])


    st <- AddModuleScore(
      object,
      features = current.list,
      pool = NULL,
      nbin = 17,
      ctrl = 100,
      k = FALSE,
      assay = DefaultAssay(object),
      name = "ModuleScore",
      seed = 1,
      search = FALSE,
      verbose = F
    )

    if (do.snip){
      module.scores[ ,i] <- snip(st@meta.data[["ModuleScore1"]])
    }

  }

    })
  })

  colnames(module.scores) <- names(gbm.list)
  rownames(module.scores) <- rownames(object@meta.data)
  rm(st)


  miko_message("Assigning GBM subtypes...", verbose = verbose)
  MES.score <- apply(module.scores[ ,c("MES1", "MES2")], 1, mean)
  NPC.score <- apply(module.scores[ ,c("NPC1", "NPC2")], 1, mean)

  pooled.scores <- cbind(MES.score, NPC.score)
  colnames(pooled.scores) <- c("MES", "NPC")
  module.scores <- cbind(module.scores, pooled.scores)

  # subset module scores
  class.names <- list(
    OPC = "OPC",
    NPC = "NPC",
    AC = "AC",
    MES = "MES"
  )


  ms <- module.scores[ ,unlist(class.names)]
  opc_npc.subset <- c(class.names$OPC, class.names$NPC)
  ac_mes.subset <- c(class.names$AC, class.names$MES)

  # define y-axis
  SC_opc_npc_max <- apply(ms[,opc_npc.subset], 1, function(x) max(x))
  SC_ac_mes_max <- apply(ms[,ac_mes.subset], 1, function(x) max(x))

  # y-axis value
  D.y <- SC_opc_npc_max-SC_ac_mes_max

  # define y axis
  # x-axis sign class
  D.y.sign <- sign(D.y)

  # D > 0 class
  Dpos.x <- log((abs(ms[,class.names$NPC] - ms[,class.names$OPC]) + 1), base = 2)
  Dpos.x.class <- Dpos.x.dif <- ms[,class.names$NPC] - ms[,class.names$OPC]
  Dpos.x.class[Dpos.x.dif > 0] <- 1  #NPC
  Dpos.x.class[Dpos.x.dif < 0] <- -1 #OPC
  Dpos.x <- Dpos.x*Dpos.x.class

  # D < 0 class
  Dneg.x <- log((abs(ms[,class.names$MES] - ms[,class.names$AC]) + 1), base = 2)
  Dneg.x.class <- Dneg.x.dif <- ms[,class.names$MES] - ms[,class.names$AC]
  Dneg.x.class[Dneg.x.dif > 0] <- 1   #MES
  Dneg.x.class[Dneg.x.dif < 0] <- -1  #AC
  Dneg.x <- Dneg.x*Dneg.x.class

  df.state <- data.frame(D.y,
                         Dpos.x,
                         Dneg.x,
                         D.y.sign
  )

  try({df.state$cell <- colnames(object)}, silent = T)
  try({df.state$seurat_clusters <- object@meta.data[["seurat_clusters"]]}, silent = T)
  try({df.state$Barcode <- object@meta.data[["Barcode"]]}, silent = T)
  # deviation.score = 1-object@meta.data[["max.tumor.transfer.score"]])

  df.state$y <- df.state$D.y
  df.state$x <- df.state$Dpos.x
  df.state$x[df.state$y < 0] <- df.state$Dneg.x[df.state$y < 0]

  scale.max <- max(c(abs(df.state$x), abs(df.state$y))) * 1.1


  df.labels <- data.frame(x = c(scale.max*0.9,-scale.max*0.9,-scale.max*0.9,scale.max*0.9),
                          y = c(scale.max,scale.max,-scale.max,-scale.max),
                          label = c("NPC", "OPC", "AC", "MES"))

  miko_message("Generating plots...", verbose = verbose)
  color.pal <- "slategray" # lightgrey
  # plt.metascores
  plt.metascores <- df.state %>%
    # dplyr::arrange(deviation.score) %>%
    # , color = deviation.score , color = "Deviation\nScore"
    ggplot(aes(x, y)) +
    xlab("Relative meta-module score\n[log(|SC1-SC2|+1)]") +
    ylab("Relative meta-module score\n[log(|SC1-SC2|+1)]") +
    labs(title = "GBM Subtypes") +
    theme_miko(legend = T, center.title = T) +
    theme(panel.border = element_rect(colour = color.pal, fill=NA, size=4)) +
    annotate("rect", xmin = -scale.max, xmax = 0, ymin = -scale.max, ymax = scale.max, fill= "white")  +
    annotate("rect", xmin = 0, xmax = scale.max, ymin = 0, ymax = scale.max , fill= "white") +
    annotate("rect", xmin = 0, xmax = scale.max, ymin = -scale.max, ymax = scale.max, fill= "white") +
    annotate("rect", xmin = -scale.max, xmax = 0, ymin = 0, ymax = scale.max, fill= "white") +
    geom_hline(yintercept=0, color = color.pal, size=1.5) +
    geom_vline(xintercept=0, color = color.pal, size=1.5) +
    geom_point(size = 0.9, alpha = 1.0)  +
    geom_label(data = df.labels, aes(x = x, y = y, label = label), fill = color.pal, color="white") +
    viridis::scale_color_viridis()

  # SPECIFY GBM SUBTYPE
  df.state$GBMstate <- NA
  state.threshold <- 0
  df.state$GBMstate[df.state$x > state.threshold & df.state$y > state.threshold] <- "NPC"
  df.state$GBMstate[df.state$x > state.threshold & df.state$y < -state.threshold] <- "MES"
  df.state$GBMstate[df.state$x < -state.threshold & df.state$y > state.threshold] <- "OPC"
  df.state$GBMstate[df.state$x < -state.threshold & df.state$y < -state.threshold] <- "AC"

  df.ms <- as.data.frame(module.scores)
  df.ms$cell <- rownames(df.ms)

  df.ms$x <- object@reductions[["umap"]]@cell.embeddings[ ,1]
  df.ms$y <- object@reductions[["umap"]]@cell.embeddings[ ,2]

  df.ms$state <- df.state$GBMstate


  plt.state.umap.list <- list()
  suppressWarnings({
    suppressMessages({

  for (i in 1:(ncol(module.scores))){

    df.ms.cur <- df.ms[ ,c("x", "y", colnames(module.scores)[i])]
    colnames(df.ms.cur) <- c("x", "y", "z")

    plt.state.umap.list[[colnames(module.scores)[i]]] <- df.ms.cur %>%
      ggplot(aes(x = x, y = y, color = z)) +
      # theme_miko(center.title = T) +
      geom_point(size = autoPointSize(n.points = nrow(df.ms.cur))) +
      theme_miko(center.title = T) +
      xlab("UMAP 1") + ylab("UMAP 2") +
      # scale_color_gradient(low = "grey95", high = "tomato") +
      viridis::scale_color_viridis() +
      scale_color_gradient2(high = scales::muted("red"), low = scales::muted("blue")) +
      labs(title = colnames(module.scores)[i])
  }
    })
  })


  plt.umap.gbm.state <- df.ms %>%
    ggplot(aes(x = x, y = y, color = state)) +
    geom_point(size = 1) +
    theme_miko(legend = T, center.title = T, color.palette = "ptol") +
    xlab("UMAP 1") + ylab("UMAP 2") +
    labs(title = "GBM State")

  miko_message("Complete!", verbose = verbose)

  return(list(
    plt.metascores = plt.metascores,
    plt.umap.scores = plt.state.umap.list,
    plt.umap.state = plt.umap.gbm.state,
    df.scores = df.ms,
    df.states = df.state
  ))

}


#' Check pubmed citations for genes
#'
#' For given topic (search query), number of publications that mention specified genes (gene.query) are tallied and returned. Based on RISmed package by Stephanie Kovalchik.
#'
#' @param gene.query Genes to query. Character vector.
#' @param search.query Search query. Boolean expression are supported. Character vector. If search.query is not specified, the number of citations for each gene will be returned without constraining the tally to a specified topic.
#' @param delay Delay between queries. If database is queried too fast, error is thrown. Default delay is 0.5.
#' @param mindate Minimum of date range for search results (examples: 2002; 2002/01/01); must be supplied with maxdate.
#' @param maxdate Maximum of date range for search results; must be supplied with mindate.
#' @param verbose Print progress (T or F). Default is T.
#' @param ... additional arguments passed to \code{\link{EUtilsSummary}}
#' @name citationCheck
#' @seealso \code{\link{EUtilsSummary}}
#' @author Nicholas Mikolajewicz
#' @examples
#' m8.citations <- citationCheck(gene.query = gene.list$m8, search.query = "hematopoietic|hematogenic|HSC", delay = 0.5,
#' mindate=2000, maxdate=2021, verbose = T)
citationCheck <- function(gene.query, search.query = NULL, delay = 0.5, mindate=2000, maxdate=2021, verbose = T,  ...){

  require(RISmed, quietly = T)

  if (!is.null(search.query)){
    res <- EUtilsSummary(query = search.query, type="esearch", db="pubmed", datetype='pdat',
                         mindate=mindate, maxdate=maxdate, ...)

    miko_message(paste0("Retrieving ",  res@count, " citations that match '", search.query, "' query..."), verbose = verbose)
    res <- EUtilsSummary(query = search.query, type="esearch", db="pubmed", datetype='pdat',  retmax=res@count,
                         mindate=mindate, maxdate=maxdate, ...) #, retmax=500
  } else {
    res <- NULL
  }


  df.pubmed <- NULL
  miko_message(paste0("Retrieving citations for ", length(gene.query) , " gene queries..."),verbose =  verbose)
  for (i in 1:length(gene.query)){

    try({

    res.gene <- EUtilsSummary(query = gene.query[i], type="esearch", db="pubmed", datetype='pdat',
                              mindate=mindate, maxdate=maxdate)

    miko_message(paste0("Retrieving ",  res.gene@count, " citations for'", gene.query[i], "' (", i, "/", length(gene.query), ")"), verbose = verbose)

    if (res.gene@count > 1000){
      res.gene <- EUtilsSummary(query = gene.query[i], type="esearch", db="pubmed", datetype='pdat',
                                mindate=mindate, maxdate=maxdate, retmax=res.gene@count) #, retmax=500
    }


    if (!is.null(search.query)){
    n.intersect <- intersect(unique(res@PMID), unique(res.gene@PMID))
    } else {
      n.intersect <- unique(res.gene@PMID)
    }


    df.pubmed <- bind_rows(df.pubmed, data.frame(
      gene = gene.query[i],
      n.publications = length(n.intersect)
    ))

    Sys.sleep(delay)


    })
  }

  plt.pubmed <- df.pubmed %>%
    ggplot(aes(x = (n.publications), y = reorder(gene, n.publications)))  +
    geom_bar(stat = "identity") +
    theme_miko() +
    scale_x_log10() +
    xlab("N Publications") +
    ylab("Gene")

  miko_message("Complete!",verbose =  verbose)
  return(
    list(
      # gene.search = res,
      string.search = res,
      df.pubmed = df.pubmed,
      plt.pubmed = plt.pubmed
    )
  )


}


#' Identify artifact genes
#'
#' For each gene in each dataset, compute how many cells in that dataset have more UMI than a specified threshold (umi.count.threshold) of that gene. A plot comparing the largest number across data sets with the third largest number is generated. For the majority of genes, these values are expected to be similar, and therefore lie on a diagonal. Genes that exhibit differences (difference.threshold) between the largest and third-largest number of cells are flagged and returned. This approach is adopted from work by Lause, Berens, Kobak (2021) BioRxiv (See Figure S6)
#'
#' @param object seurat object or named list of seurat objects.
#' @param assay assay containing count matrix. If unspecified, "RNA" assay is used. If "RNA" assay is missing, the default assay is used.
#' @param features genes used for analysis. If unspecified, variable genes are used.
#' @param meta.feature feature in meta data that provides dataset grouping information. Default is "Barcode".
#' @param umi.count.threshold cells that exceed this UMI count threshold are counted.
#' @param difference.threshold fold difference between largest and 3rd largest number of cells that is used to flag artifact genes for omission.
#' @param group.specific.is.artefact include genes that are only expressed in one dataset as artifacts. Default is FALSE.
#' @param verbose Print progress (T or F). Default is F.
#' @name findArtifactGenes
#' @author Nicholas Mikolajewicz
#' @examples
#' ag.res <-  findArtifactGenes(object = so, assay = NULL, features = NULL, meta.feature = "Barcode", umi.count.threshold = 5, difference.threshold = 100, verbose = T)
findArtifactGenes <- function(object, assay = NULL, features = NULL, meta.feature = "Barcode", umi.count.threshold = 5, nth.max = 2,
                              difference.threshold = 30, group.specific.is.artefact = F, verbose = F){


  if (class(object) == "list"){
    which.method <- 1
    object.list <- object
    object <- object[[1]]
    if (is.null(features)) stop("If object is list, features must be specified")
  } else {
    which.method <- 2
  }
  stopifnot(class(object) == "Seurat")

  if (is.null(assay)) {
    all.assays <- names(object@assays)
    if ("RNA" %in% all.assays) {
      assay <- "RNA"
    } else {
      assay = DefaultAssay(object)
    }
  }
  if (is.null(features)) {
    features = VariableFeatures(object)
    if (length(features) == 0){
      stop("Variable features not found Please provide features for evaluation.")
    }
  }



  if (which.method == 1){
    u.bc <- names(object.list)
    if (length(u.bc) < nth.max) stop("Insufficient number of datasets to perform artifact detection")
    cumi.mat <- matrix(nrow = length(features), ncol = length(u.bc))
    for (i in 1:length(u.bc)){
      object <- object.list[[u.bc[i]]]
      e.mat <- object@assays[[assay]]@counts
      if (is.null(dim(e.mat))) stop("Count matrix not found. Cannot perform artefact detection.")
      e.mat <- e.mat[rownames(e.mat) %in% features, ]
      cumi.mat[,i] <-  apply(e.mat, 1, function(x) sum(x>umi.count.threshold))
    }

    rownames(cumi.mat) <- rownames(e.mat);
    colnames(cumi.mat) <- u.bc


  } else if (which.method == 2){
    miko_message("Preparing expression matrix...", verbose = verbose)
    # get grouping variable
    u.bc <- unique(object@meta.data[ ,meta.feature])
    if (length(u.bc) < nth.max) stop("Insufficient number of datasets to perform artifact detection")
    all.bc <- as.character(object@meta.data[ ,meta.feature])

    # get expression matrix
    e.mat <- object@assays[[assay]]@counts
    if (is.null(dim(e.mat))) stop("Count matrix not found. Cannot perform artefact detection.")
    e.mat <- e.mat[rownames(e.mat) %in% features, ]

    miko_message(paste0("Identifying how many cells have over ", umi.count.threshold, " UMI for each gene...."), verbose = verbose)
    df.cell.umi <- NULL
    cumi.mat <- matrix(nrow = nrow(e.mat), ncol = length(u.bc))

    # flag cells with expression above threshold
    for (i in 1:length(u.bc)){
      e.subset <- e.mat[ ,all.bc %in% u.bc[i]]
      cumi.mat[,i] <-  apply(e.subset, 1, function(x) sum(x>umi.count.threshold))
    }
    rownames(cumi.mat) <- rownames(e.mat);
    colnames(cumi.mat) <- u.bc

  }

  # rank datasets
  cumi.rank <- t(apply(cumi.mat, 1, rank, ties.method = "random"))

  # get largest and 3rd largest number of cells
  df.cumi.ind <- data.frame(
    x = apply(cumi.rank, 1, function(x) which(x == length(u.bc))),
    y = apply(cumi.rank, 1, function(x) which(x == length(u.bc)-(nth.max-1)))
  )

  miko_message(paste0("Consolidating results..."), verbose = verbose)
  # condsolidate data
  df.cumi.n <- NULL
  for (i in 1:nrow(df.cumi.ind)){
    df.cumi.n <- bind_rows(
      df.cumi.n,
      data.frame(
        x =  cumi.mat[i ,df.cumi.ind$x[i]],
        y =  cumi.mat[i ,df.cumi.ind$y[i]] )

    )
  }

  df.cumi.n$gene <- rownames(df.cumi.ind)

  if (nth.max == 2){
    nth.label <- "2nd"
  } else if (nth.max == 3){
    nth.label <- "3rd"
  } else if (nth.max == 4){
    nth.label <- "4th"
  } else if (nth.max == 5){
    nth.label <- "5th"
  } else if (nth.max == 6){
    nth.label <- "6th"
  } else {
    nth.label <- paste0(nth.max, "th")
  }

  miko_message(paste0("Generating plot..."), verbose = verbose)
  # generate plot
  df.cumi.n$ratio <- (df.cumi.n$x + 1)/(df.cumi.n$y + 1)
  df.cumi.n$is.artifact <- df.cumi.n$ratio >= difference.threshold

  if (group.specific.is.artefact){
    df.cumi.n$is.artifact2 <- F
    df.cumi.n$is.artifact2[df.cumi.n$y == 0]  <- T
  } else {
    df.cumi.n$is.artifact2 <- F
  }

  p3 <- df.cumi.n %>%
    ggplot(aes(x = x + 1, y = y + 1, color = is.artifact | is.artifact2))+
    geom_point() +
    # geom_path(data = df.threshold, aes(x = x, y = y)) +
    ggrepel::geom_text_repel(data = df.cumi.n %>% dplyr::filter(is.artifact), aes(x = x+ 1, y = y + 1, label = gene), max.overlaps = Inf, size = 2) +
    scale_x_log10() +
    scale_y_log10() +
    labs(x = paste0("1 + max(#cells with > ", umi.count.threshold, " counts)"),
         color = "is.artifact",
         caption = paste0(sum(df.cumi.n$is.artifact | df.cumi.n$is.artifact2), " artifact genes detected") ,
         y = paste0("1 + ", nth.label, "-max(#cells with > ", umi.count.threshold, " counts)"),
         title = "Artifact gene detection",
         subtitle = paste0("Difference Threshold: ", difference.threshold)) +
    theme_miko(color.palette = "ptol", legend = T)



  artifact.gene <- df.cumi.n$gene[(df.cumi.n$is.artifact | df.cumi.n$is.artifact2)]

  miko_message(paste0("Complete!"), verbose = verbose)

  return(list(
    artifact.gene = artifact.gene,
    data = df.cumi.n,
    plot = p3
  ))

}


#' Identify optimal bin size for AddModuleScore() function
#'
#' Identify optimal bin size for AddModuleScore() function
#'
#' @param object Seurat object
#' @param pool List of features to check expression levels against. Defaults to rownames(x = object).
#' @param nbinNumber of bins of aggregate expression levels for all analyzed features. Initial value to begin with looking for optimal bin.
#' @param seed See a random seed. If NULL, seed is not set.
#' @name optimalBinSize
#' @author Nicholas Mikolajewicz
#' @examples
#' optimal.nbin <-  optimalBinSize(object = so.query)
optimalBinSize <- function (object, pool = NULL, nbin = 24, seed= 1023, verbose = T){


  if (!is.null(x = seed)) {
    set.seed(seed = seed)
  }

  assay.data <- GetAssayData(object = object)

  pool <- pool %||% rownames(x = object)
  data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
  data.avg <- data.avg[order(data.avg)]


  StartNBin = nbin
  data.cut <- NA
  while(class(data.cut)=="try-error" || is.na(data.cut)[1]){

    data.cut <- try({cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30,
                                n = StartNBin, labels = FALSE, right = FALSE)}, silent = T)

    if (class(data.cut)=="try-error" || is.na(data.cut)[1]){
      StartNBin = round(StartNBin-1)
    }



  };

  miko_message(paste0("Optimal bin size: ",StartNBin), verbose = verbose)

  return(StartNBin)

}



#' Convert long data frame to named list
#'
#' Convert long data frame to named list. Input long data frame consists of two columns, the first corresponding to the names within the list (group_by), and the second to the corresponding list entries (values).
#'
#' @param df.long long data frame
#' @param group_by name of data frame column to group values by
#' @param values name of data frame column containing values
#' @name longDF2namedList
#' @return named list
#' @author Nicholas Mikolajewicz
#' @examples
#'
longDF2namedList <- function(df.long, group_by, values){

  if (!("data.frame" %in% class(df.long))) {
    try({df.long <- as.data.frame(df.long)}, silent = T)
    if (!("data.frame" %in% class(df.long))) {
      stop("Input must be a data frame")
    }
  }


  n.list <- list()
  ugroup <- as.character(unique(unlist(df.long[ ,group_by])))
  ugroup <- ugroup[order(ugroup)]
  for (i in 1:length(ugroup)){

    # col.name <- colnames(df.wide)[i]
    entries <- unique(unlist(df.long[ which(as.vector(unlist(df.long[ ,group_by])) %in% ugroup[i]) , values] ))
    entries <- entries[!is.na(entries)]
    entries <- entries[entries != ""]

    if (is.factor(entries)) entries <- as.character(entries)
    n.list[[ugroup[i]]] <- entries

  }

  return(n.list)
}

