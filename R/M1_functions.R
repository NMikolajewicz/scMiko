


#' Load CellRanger preprocessed data
#'
#' Input data is output from cell ranger (10x datasets)
#'
#' @param import_set Character specifying folder containing matrix.mtx, genes.tsv (or features.tsv) and barcodes.tsv files provided by 10x.
#' @param input_organisms Character specifying species to include. This is necessary to parse ensemble ids into mouse and human. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item c("Hs", "Mm") - both species included
#' }
#' @param dir Character. folder containing import_set folder.
#' @name loadCellRanger
#' @return list containing Seurat Object and named gene vector.
#'
loadCellRanger <- function(import_set, input_organisms, dir = "") {

  import_set_path <- paste(dir, import_set, sep ="")

  # import cell ranger-processed data
  expression_matrix <- Seurat::Read10X(data.dir = import_set_path[1])

  # assign to secondary matrix
  expression_matrix2 <- expression_matrix

  # import gene and ensembl names
  load.success <- F
  try({
    feature.names = read.delim(paste(import_set_path[1], "/genes.tsv", sep = ""),
                               header = FALSE,
                               stringsAsFactors = FALSE)
    load.success <- T
  }, silent = T)
  if (!load.success){
    feature.names = read.delim(paste(import_set_path[1], "/features.tsv", sep = ""),
                               header = FALSE,
                               stringsAsFactors = FALSE)
  }


  # remove tags
  feature.names[,1] <- gsub("hg19_", "", feature.names[,1])
  feature.names[,2] <- gsub("hg19_", "", feature.names[,2])
  feature.names[,1] <- gsub("hg38_", "", feature.names[,1])
  feature.names[,2] <- gsub("hg38_", "", feature.names[,2])


  # create gene list
  gNames <- as.character( feature.names$V2 );
  names(gNames) <- as.character( feature.names$V1 );
  names(gNames) <-gsub("\\..*","",as.vector( names(gNames)))

  # convert gene symbols to ensemble (in expression matrix)
  rownames(expression_matrix2) <- names(gNames)

  # assign species
  orgs <- detectSpecies(expression_matrix2)
  # orgs <- scMiko::m1.inferSpecies(expression_matrix2, input_organisms)

  # create seurat object
  # so = CreateSeuratObject(counts = expression_matrix2,project="cell_ranger",min.cells=3,min.features=200)
  so = CreateSeuratObject(counts = expression_matrix2,project="cell_ranger",min.cells=0,min.features=0)

  # add gene symbols as meta data in seurat object
  mat_ens <- rownames(so@assays[["RNA"]])
  match.id <- match(mat_ens, names(gNames))

  # Add gene symbols as meta data that we can use later
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=names(gNames[match.id]),col.name="ENSEMBL");
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=as.vector(gNames[match.id]),col.name="SYMBOL");

  # Add in inferred organism
  so$Species <- orgs;

  # Specify barcodes
  so$Barcode <- "unspecified";

  output <- list(so, gNames)

  return(output)
}


#' Load preprocessed data from Moffat lab sciRNA-seq3 pipeline
#'
#' Load preprocessed data from Moffat lab sciRNA-seq3 pipeline. RT barcode and plate summary is stored in misc slot of resulting seurat object.
#'
#' @param import_set Character vector specifying expression matrix (import_set[1]), PCR barcodes (import_set[2]) and RC barcodes (import_set[3]). Expression matrix will be imported successfully if barcodes are omitted.
#' @param subsample_factor Numeric [0,1]. Subsampling factor
#' @param input_organisms All species included in input files. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item c("Hs", "Mm") - both species included
#' }
#' @param organism_include Species to include in downstream analysis. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item c("Hs", "Mm") - both species included
#' }
#' @param dir Character. folder containing import_set files
#' @name loadMoffat
#' @return list containing Seurat Object and named gene vector.
#'
loadMoffat <-function(import_set, subsample_factor, input_organisms, organism_include, dir) {

  # load gene count matrix
  import_set_path <- paste(dir, import_set, sep ="")
  load(import_set_path[1])

  # assign row (genes) and column (cells) names
  gene_count2 <- gene_count

  # filter out incorrect species genes immediately.
  all.genes <- rownames(gene_count2)
  include.which.genes <- rep(FALSE, length(all.genes))
  if ((length(input_organisms) > length(organism_include)) | (length(organism_include) == 1)){
    for (i in 1:length(organism_include)){
      if (organism_include[i] == "Mm"){
        include.which.genes[grepl("MUSG", all.genes)] <- T
      } else if (organism_include[i] == "Hs"){
        include.which.genes[!(grepl("MUSG", all.genes))] <- T
      }
    }

    df_gene <- df_gene[include.which.genes, ]
    gene_count <- gene_count[include.which.genes, ]
    gene_count2 <- gene_count2[include.which.genes, ]
  }

  # Infer organism
  # orgs <- scMiko::inferSpecies(gene_count2, input_organisms)
  orgs <-  detectSpecies(gene_count2)
  df_cell$orgs <- orgs

  # only keep protein coding genes

  # assign gene and cell names
  rownames(gene_count2) <- df_gene$gene_id
  colnames(gene_count2) <- df_cell$sample


  # subsample data
  if (subsample_factor < 1){
    col_ind <- c(1:dim(gene_count2)[2])
    subsample_ind <- sample(col_ind, round(subsample_factor * dim(gene_count2)[2]), replace = FALSE)
    gene_count2 <- gene_count2[,subsample_ind]
  } else {
    subsample_ind <- c(1:dim(gene_count2)[2])
  }

  # Need names vectors to add additional metadata to Seurat object
  gNames <- as.character( df_gene$gene_name );
  names(gNames) <- as.character( df_gene$gene_id );
  names(gNames) <-gsub("\\\\..*","",as.vector( names(gNames)))

  # Add PCR barcodes
  pcr.barcode.flag <- F
  if (dir != import_set_path[2]){
    try({
      grps <- read.csv(import_set_path[2], header=T,sep="\\t",stringsAsFactors=F);
      wells <- gsub("Han_[0-9]+_([A-Z0-9]{2,3}).[ACGT]+$","\\\\1",colnames(gene_count2));
      myGrps <- grps$ConditionGroup[ match(wells,grps$sampleWell) ];
      names(myGrps) <- colnames(gene_count2)
      myGrps <- myGrps[subsample_ind]
      pcr.barcode.flag <- T
    }, silent = T)
  }


  # Add RT barcode
  rtBC <- tryCatch({
    rtBC <- read.csv(import_set_path[3],header=T,sep="\\t",stringsAsFactors=F);
  }, error = function(e){
    rtBC <- read.csv2(import_set_path[3],header=T,sep="\t",stringsAsFactors=F)
    return(rtBC)
  })
  rtBC <- rtBC[subsample_ind, ]
  if (is.null(rtBC$Plate)) rtBC$Plate <- ""
  rtBC$plate.well.sample <- paste0("P", rtBC$Plate, ".", rtBC$Well, ".", rtBC$Sample.type)


  barcodes <- gsub(".+([ACGT]{10}$)","\\\\1",colnames(gene_count2));
  if (unique(barcodes) == "\\1")  barcodes <- gsub(".+([ACGT]{10}$)","\\1",colnames(gene_count2));
  idx <- match( barcodes,rtBC$rt.bc );
  rtGrp <- rtBC$Sample.type[idx];
  names(rtGrp) <- colnames(gene_count2);
  rtPWS <- data.frame(
    plate = rtBC$Plate[idx],
    well = rtBC$Well[idx],
    sample = rtBC$Sample.type[idx])
  rtPWS.summary <- rtPWS %>% dplyr::group_by(plate, well, sample) %>% dplyr::tally()

  # remove .* suffix in ensembl
  rownames(gene_count2)<-gsub("\\\\..*","",as.vector( rownames(gene_count2)))

  # Create seurat object with count matrix data
  so <- CreateSeuratObject(counts=gene_count2,project="Hong-sciSeq3",min.cells=0,min.features=0,names.field=2,names.delim="." )

  # Add gene symbols as meta data that we can use later
  mat_ens <- rownames(so@assays[["RNA"]])
  match.id <- match(mat_ens,  names(gNames))

  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=names(gNames)[match.id],col.name="ENSEMBL");
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=as.vector((gNames)[match.id]),col.name="SYMBOL");

  # assign meta data
  so.meta <- so@meta.data
  so.meta$ind <- seq(1, nrow(so.meta))
  so.meta$sample <- rownames(so.meta)
  so.merge <- merge(so.meta, df_cell, by = "sample", all.x = T)
  so.merge <- so.merge %>% dplyr::arrange(ind)
  if ("unmatched_rate" %in% colnames(df_cell))  so[["unmatched.rate"]] <- so.merge$unmatched_rate
  so[["Species"]] <- so.merge$orgs

  if (length(import_set) > 1){ # same as above, if barecode and well info. provided, assign to so structure
    # Add group #s for mapping cells to PCR condition groups
    if (pcr.barcode.flag) {
      so[["group"]] <- myGrps;
    } else {
      so[["group"]] <- "unspecified";
    }

    # Add in cell line group info
    so$Barcode <- rtGrp;

  }

  try({so@misc[["gene.info"]] <- df_gene}, silent = T)
  try({so@misc[["cell.info"]] <- df_cell}, silent = T)
  try({so@misc[["plate.summary"]] <- rtPWS.summary}, silent = T)

  # return output
  output <- list(so, gNames)

  return(output)
}


#' Load gene x cell count matrix into seurat object
#'
#' Load gene x cell count matrix into seurat object. For intended performance, first column of spreadsheet must contain gene names, and first row of spreadsheet must contain cell names.
#'
#' @param file Matrix file name. (.csv, .tsv supported)
#' @param dir Character. folder containing import_set file
#' @name loadMat
#' @return list containing Seurat Object and named gene vector.
#'
loadMat <- function(file, dir) {

  import_set_path <- paste(dir, file, sep ="")

  if (grepl(".rds", import_set_path[1])){
    expression_matrix <- readRDS(import_set_path[1])
    expression_matrix[is.na(expression_matrix)] <- 0
    feature.names <- rownames(expression_matrix)
    expression_matrix$GENE <- feature.names
  } else {
    expression_matrix <- data.table::fread(import_set_path[1], header = T, stringsAsFactors = F, showProgress = T)
    first_col <- colnames(expression_matrix)[1]
    expression_matrix <- col2rowname(expression_matrix, first_col)
    expression_matrix$GENE <- rownames(expression_matrix)
    expression_matrix[is.na(expression_matrix)] <- 0
    feature.names <- as.vector(expression_matrix$GENE)
  }


  # average duplicate rows (i.e., duplicate genes)
  if (length(unique(feature.names)) < length(feature.names)){

    success.attempt <- F
    try({

      expression_matrix.col <- expression_matrix
      expression_matrix.col <- dplyr::select(expression_matrix.col, -c("GENE"))
      expression_matrix.col$rowID <- seq(1, nrow(expression_matrix.col))
      expression_matrix.col.noDup <-  WGCNA::collapseRows(expression_matrix.col, rowGroup = expression_matrix$GENE, rowID = expression_matrix.col$rowID, method = "Average")
      new.mat <- as.data.frame(expression_matrix.col.noDup[["datETcollapsed"]])
      new.mat <- dplyr::select(new.mat, -c("rowID"))
      new.mat$GENE <- rownames(new.mat)
      expression_matrix2 <- new.mat
      rm("new.mat")
      success.attempt <- T

    }, silent = T)

    if (!success.attempt){
      which.duplicate <- duplicated(feature.names)
      expression_matrix2 <- expression_matrix[!(which.duplicate), ]
    }
  } else {
    # assign to secondary matrix
    expression_matrix2 <- expression_matrix
  }


  feature.names <- as.vector(expression_matrix2$GENE)
  expression_matrix2 <- dplyr::select(expression_matrix2, -c("GENE"))
  rownames(expression_matrix2) <- feature.names

  # if (length(input_organisms) != 1){
  ens.sum <-  sum(grepl("ENS", feature.names))
  ens.mus.sum <-  sum(grepl("ENSMUS", feature.names))
  hi.cap.sum <-  sum(feature.names == toupper(feature.names))
  lo.cap.sum <-  sum(feature.names == firstup(feature.names))

  df.rep <-as.data.frame(t(data.frame(
    ens.sum = ens.sum,
    ens.mus.sum = ens.mus.sum,
    hi.cap.sum = hi.cap.sum,
    lo.cap.sum = lo.cap.sum
  )))

  which.rep <- rownames(df.rep)[which.max(df.rep[, 1])]

  if ((which.rep %in% c("ens.sum", "ens.mus.sum")) & (ens.sum > ens.mus.sum)) {
    species <- "Hs"
    rep_gene <- "ENS"
  } else if ((which.rep %in% c("ens.sum", "ens.mus.sum")) & (ens.sum <= ens.mus.sum)) {
    species <- "Mm"
    rep_gene <- "ENS"
  } else if (which.rep == "hi.cap.sum") {
    species <- "Hs"
    rep_gene <- "SYMBOL"
  } else if (which.rep == "lo.cap.sum") {
    species <- "Mm"
    rep_gene <- "SYMBOL"
  }

  input_organisms <- species
  # }

  if (rep_gene == "SYMBOL"){
    g2eNames <- sym2ens(my.symbols =  feature.names, my.species = input_organisms)
    gNames <- g2eNames$SYMBOL
    names(gNames) <- g2eNames$ENSEMBL
  } else if (rep_gene == "ENS"){
    e2gNames <- ensembl2sym(my.ensembl =  feature.names, my.species = input_organisms)
    gNames <- e2gNames$SYMBOL
    names(gNames) <- e2gNames$ENSEMBL
  }
  gNames <- gNames[!is.na(gNames)]
  names(gNames) <-gsub("\\..*","",as.vector( names(gNames)))


  # assign organism
  orgs <- rep(input_organisms,ncol(expression_matrix2));
  names(orgs) <- colnames(expression_matrix2);


  # create seurat object
  so = CreateSeuratObject(counts = expression_matrix2, min.cells=0, min.features=0)

  # add gene symbols as meta data in seurat object
  mat_ens <- rownames(so@assays[["RNA"]])
  match.id <- match(mat_ens, gNames)
  gNames_filtered <- names(gNames)[match.id]

  # Add gene symbols as meta data that we can use later
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=gNames_filtered,col.name="ENSEMBL");
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=as.vector((gNames)[match.id]),col.name="SYMBOL");

  # Add in inferred organism
  so$Species <- orgs;

  # Specify barcodes
  so$Barcode <- "unspecified";

  return(list(object = so,
              genes = gNames))
}



#' Subset and assign labels to seurat
#'
#' Subset and assign barcode labels (in Seurat object) to "subset_group" metadata field, as specified by which.strata parameter.
#'
#' @param so Seurat Object
#' @param which.strata Barcode labeling parameter. If NA, "subset_group" metadata field is set to "pooled".
#' @name barcodeLabels
#' @return Seurat Object (with updated metadata)
#'
barcodeLabels <- function(so, which.strata) {
  # set Seurat subset labels for cells of interest
  if (!is.na(which.strata)){
    pattern <- paste(which.strata, collapse="|")
    # so <- Seurat::SubsetData(object = so, cells = (grepl(pattern, so@meta.data[["Barcode"]])))
    so <- so[ , grepl(pattern, so@meta.data[["Barcode"]])]

    for (i in 1:length(which.strata)){
      so@meta.data[["subset_group"]][grepl(which.strata[i], so@meta.data[["Barcode"]])] <- which.strata[i]
    }
  } else {
    so@meta.data[["subset_group"]] <- "pooled"
  }

  return(so)

}


#' Calculate mitochondrial content
#'
#' Calculate cell-level mitochondrial content based on mapped transcripts and store in Seurat "percent.mt' metadata field.
#'
#' @param so Seurat Object
#' @param gNames Named gene list; entries are Symbols, names are Ensemble.
#' @param assay assay to use.
#' @param omit.na Logical specifying whether to omit NA entries (present when unfiltered 10x dataset is used). Default is True.
#' @name getMitoContent
#' @return Seurat Object (with updated metadata)
#'
getMitoContent <- function(so, gNames = NULL, assay = NULL, omit.na = T) {

  if (!("Seurat" %in% class(so))) stop("'object' does not belong to Seurat class")

  assay <- DefaultAssay(so)
  if (!is.null(gNames)){
    f <- intersect(names(gNames)[ grep("^(MT|mt)-",gNames) ],rownames(so[[assay]]))
  } else {
    f <- intersect(rownames(so)[ grep("^(MT|mt)-",rownames(so)) ],rownames(so[[assay]]))
  }

  if (length(f) == 0){
    so[["percent.mt"]] <- 0
  } else {
    so[["percent.mt"]] <- PercentageFeatureSet(so, features=f)
  }


  if (omit.na) {
    # omit nan entries - they occur when cell ranger output data is unfiltered.
    so <-so[ , !(is.nan(so@meta.data[["percent.mt"]]))]
  }

  return(so)

}


#' Apply QC filters to Seurat Object
#'
#' Filter scRNAseq data by gene/cell, unmatch rate and mitochondrial content.
#'
#' @param so Seurat Object
#' @param RNA.upperlimit Upper limit threshold for genes/cell
#' @param RNA.lowerlimit Lower limit threshold for genes/cell
#' @param mt.upperlimit Numeric [0, 100]. Upper limit threshold for microchondiral content.
#' @param unmatch.low Numeric [0,1]. Lower limit threshold for unmatch rate. Default is 0. Ignored if data is not from sciRNA-seq3 pipeline.
#' @param unmatch.high Numeric [0,1]. Upper limit threshold for unmatch rate. Default is 1. Ignored if data is not from sciRNA-seq3 pipeline.
#' @param set_names Character specifying dataset name. Optional.
#' @name filterSeurat
#' @author Nicholas Mikolajewicz
#' @return List containing seurat object, filter summary statistics, and filter breakdown list
#'
filterSeurat <- function(so, RNA.upperlimit = 9000, RNA.lowerlimit = 200, mt.upperlimit = 60, unmatch.low = 0, unmatch.high = 1, set_names = NULL) {

  if (!("Seurat" %in% class(so))) stop("'object' does not belong to Seurat class")
  if (!("percent.mt" %in% colnames(so@meta.data))) so <- getMitoContent(so)

  # determine unfiltered UMI count
  original_count <- length(so@meta.data[["nCount_RNA"]])

  # breakdown of filter
  df.criteria <- data.frame(cells = colnames(so), mito.content = so$percent.mt, gene.count = so$nFeature_RNA)
  if ("unmatched.rate" %in% names(so@meta.data)) {
    df.criteria$u.rate <- so@meta.data[["unmatched.rate"]]
  } else {
    df.criteria$u.rate <- 0
  }

  filter.list <- list(
    mito.content = df.criteria$cells[df.criteria$mito.content > mt.upperlimit],
    gene.count = df.criteria$cells[(df.criteria$gene.count > RNA.upperlimit) | (df.criteria$gene.count < RNA.lowerlimit)],
    unmatched.rate = df.criteria$cells[(df.criteria$u.rate > unmatch.high) | (df.criteria$u.rate < unmatch.low)]
  )


  # filter dataset
  if ((unmatch.low == 0) & (unmatch.high == 1)){
    so <- subset(so, subset = ((nFeature_RNA <= RNA.upperlimit) &
                                 (nFeature_RNA >= RNA.lowerlimit) &
                                 (percent.mt <= mt.upperlimit)))
  } else {
    so <- subset(so,
                 subset = ((nFeature_RNA < RNA.upperlimit) &
                             (nFeature_RNA >= RNA.lowerlimit) &
                             (percent.mt <= mt.upperlimit) &
                             (unmatched.rate >= unmatch.low) &
                             (unmatched.rate <= unmatch.high)))
  }

  # determine filtered UMI count
  filtered_count = length(so@meta.data[["nCount_RNA"]])
  percent_remaining <- 100*filtered_count/original_count

  filter_summary <- data.frame(pre.filtering = original_count, post.filtering = filtered_count, per_remaining = percent_remaining)

  if (is.null(set_names)) set_names <- "scRNAseq"

  # generate pre- post- filter statistics plot
  filter_summary_list <- filter_summary
  filter_summary_list$set <- replicate(dim(filter_summary_list)[1], (set_names))
  filter_summary <- filter_summary_list[, c(4, 1, 2, 3)]
  filter_summary_gathered <- gather(filter_summary, key = "filter", value = "count", c("pre.filtering", "post.filtering"))
  filter_summary_gathered$filter <- factor(filter_summary_gathered$filter, levels = c("pre.filtering", "post.filtering"))
  per.omitted <- signif((100 - as.numeric(unique(filter_summary_gathered$per_remaining))), 3)
  plt.filter_pre_post <- ggplot(filter_summary_gathered, aes(x = set, y = count, fill = filter )) +
    geom_bar(stat="identity", position=position_dodge()) +
    ylab("Cell Count") +
    ggtitle(paste(per.omitted, "% cells omitted", sep = ""))

  output <- list(seurat = so,
                 plot = plt.filter_pre_post,
                 filter.breakdown = filter.list)
  return(output)

}


#' Normalize and Scale scRNAseq Data
#'
#' Applies one of two normalization/scaling approaches supported by Seurat.
#'
#' @param object Seurat Object
#' @param method Character specifying data normalization and scaling method. One of:
#' \itemize{
#' \item "NFS" - Seurat's NormalizeData, FindVariableFeatures, ScaleData workflow. Parameters are set to use LogNormalization method with a scale.factor of 1000. Variable features are selceted using 'mvp' method, and var2regress is regressed out during data scaling.
#' \item "SCT" - Default. SCTransform workflow; Uses regularized negative binomial regression to normalize UMI count data. All genes are returned (not only variable), and residual variace cutoff is set to 1.3; var2regress is included in model.
#' }
#' @param vars2regress Character vector specifying which variables to regress out during data scaling.
#' @param enable.parallelization Logical specifying whether to enable parallelization. Default is T.
#' @param n.works Number of works to used during parallel processing. Default is 1.
#' @param max.memory Max memory to use during parallel processing. Default is 20480 * 1024^2
#' @param variable.features.n If method = SCT, integer that specifies number of variable features to retrieve. If specified, overides variable feature threshold specified by variable.features.rv.th.
#' @param variable.features.rv.th If method = SCT, numerical that specifies residual variance theshold for variable features. Default is 1.3.
#' @param return.only.var.genes If method = SCT, logical that specifies whether only variable genes are retrieved.
#' @param mean.cutoff If method = NFS, a two-length numeric vector with low- and high-cutoffs for feature means.
#' @param dispersion.cutoff If method = NFS, a two-length numeric vector with low- and high-cutoffs for feature dispersions.
#' @param conserve.memory If set to TRUE the residual matrix for all genes is never created in full; useful for large data sets, but will take longer to run; this will also set return.only.var.genes to TRUE; default is FALSE.
#' @param assay Name of assay to pull the count data from; default is 'RNA'
#' @param verbose print progress report. Default is FALSE.
#' @param ... additional parameters passed to SCTransform.
#' @name scNormScale
#' @return Seurat Object
#'
scNormScale <- function(object,  method = "SCT", vars2regress = NULL, enable.parallelization = T, n.workers = 1, max.memory = (20480 * 1024^2), variable.features.n = NULL,
                           variable.features.rv.th = 1.3, return.only.var.genes = F, mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf), conserve.memory = FALSE, assay = "RNA", verbose = F, ...){


  stopifnot("'object' is not Seurat object"= ("Seurat" %in% class(object)))

  # enable parallelization
  if (enable.parallelization){
    future::plan(strategy = "multisession", workers = n.workers)
    options(future.globals.maxSize = max.memory)
  }


  # Normalize and scale data
  if (method == "NFS"){

    stopifnot(length(mean.cutoff) == 2)
    stopifnot(length(dispersion.cutoff) == 2)

    # Normalize data
    object <- NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, assay = assay)
    # Find variable features
    object <- FindVariableFeatures(object = object, selection.method = 'mvp', mean.cutoff = mean.cutoff, dispersion.cutoff = dispersion.cutoff, assay = assay)
    # Scale data
    if (is.null(vars2regress)){
      object <- ScaleData(object, features = rownames(object), assay = assay)
    } else {
      object <- ScaleData(object, features = rownames(object), vars.to.regress = vars2regress, assay = assay)
    }


  } else if (method == "SCT"){

    # apply sctransform (regularized negative binomial regression)
    # also removes confounding source of variation (i.e., mitochonrdial mapping percentage)
    if (is.null(vars2regress)){
      object <- SCTransform(object = object,
                        verbose = verbose,
                        return.only.var.genes = return.only.var.genes,
                        variable.features.n = variable.features.n,
                        variable.features.rv.th = variable.features.rv.th,
                        assay = assay,
                        method = "glmGamPoi",
                        conserve.memory = conserve.memory,
                        ...)
    } else {
      object <- SCTransform(object = object,
                        vars.to.regress = vars2regress,
                        verbose = verbose,
                        return.only.var.genes = return.only.var.genes,
                        variable.features.n = variable.features.n,
                        variable.features.rv.th = variable.features.rv.th,
                        assay = assay,
                        method = "glmGamPoi",
                        conserve.memory = conserve.memory,
                        ...)
    }


  }

  return(object)
}



#' Infer species from list of Ensemble ids
#'
#' For given gene expression matrix, assign species to each cell (column) based on Ensemble IDs.
#'
#' @param exp.mat Expression matrix (genes x cell); row names are ensemble IDs, col names are cell IDs.
#' @param expected.species Character vector specifying which species are expected in expression matrix.
#' \itemize{
#' \item "Mm" - Mouse
#' \item "Hs" - Human
#' \item c("Mm", "Hs") - Mouse and Human.
#' }
#' @param rep.ens.method Method used to sample representative Ensemble ID for each cell (i.e., column).
#' \itemize{
#' \item "orig" - Original method; computation slower than alternative.
#' \item "alt" - Default. Alternative method;
#' }
#' @name inferSpecies
#' @return Character vector specifying species of each cell in expression  matrix.
#'
inferSpecies <- function(exp.mat, expected.species, rep.ens.method = "orig"){

  # Infer organism
  if (length(expected.species) > 1) {

    # get representative Ensemble for each cell
    if (rep.ens.method == "orig"){
      orgIDs <- rownames(exp.mat)[ apply( exp.mat,2,which.max )]
    } else if (rep.ens.method == "alt"){

      # SOME ISSUES WITH ALT METHOD. EITHER REMOVE OR TROUBLESHOOT (180720). ORIG WORKS.
      n <- 1000
      nc <- ncol(exp.mat)
      z <- split(colnames(exp.mat), rep(1:ceiling(nc/n), each=n, length.out=nc))
      orgIDs <- character()
      for (i in names(z)) {
        chunk <- exp.mat[,z[[i]]]
        b <- rownames(chunk)[apply(chunk,2,which.max)]
        orgIDs <- c(orgIDs,b)
        rm(list=c("chunk","b"))
      }
      rm(list=c("n","nc","z"))
    }

    orgs <- rep("Hs",ncol(exp.mat));
    orgs[grep("^ENSMUS",orgIDs) ] <- "Mm";
    names(orgs) <- colnames(exp.mat);
  } else {
    orgs <- rep(expected.species,ncol(exp.mat));
    names(orgs) <- colnames(exp.mat);
  }

  return(orgs)

  }

