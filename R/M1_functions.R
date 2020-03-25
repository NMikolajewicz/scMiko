#' Load parameter specifications for Module 1
#'
#' Load parameter specification list for Module 1 (preprocessing and QC). If not specified, default settings are applied.
#'
#' @param which.data Character. Specify input data (input file specification are stored in separate data.frame. See Module 1 for details)
#' @param which.strata Character. Specify how input data are subset. Leave unspecified to keep all data.
#' @param organism.filter.flag Logical specifying whether inputs should be filtered by species. Recommended if multiple species are present in input.
#' @param organism.include Character. Which species to include. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item "Hs, Mm" - both species included
#' }
#' @param save.flag Logical specifying whether preprocessed data are saved as RData file.
#' @param save.filename Character. Output filename (must include .Rdata suffix)
#' @param save.directory Character. output directory (Default is "Preprocessed Datasets/")
#' @param set_names Optional character. Label for dataset.
#' @param data.imputed.flag Logical specifying wehther data were imputed. Default is F.
#' @param vars2regress Character vector specifying which parameters to regress during data scaling. Default is "percent.mt"
#' @param subsample_factor Numeric [0,1]. Factor used to subsample data. Recommended during initial data exploration, especially if dataset is large.
#' @param plt.log.flag Logical specifying whether QC violin plots are on log scale.
#' @param RNA.upperlimit Numeric. Upper limit for number of genes per cell. Default is 200.
#' @param RNA.lowerlimit Numeric. Lower limit for number of genes per cell. Default is 9000.
#' @param mt.upperlimit Numeric [0,100]. Upper limit for mitochondrial percentage. Default is 60.
#' @param cluster.resolution Resolution used for cell clustering.
#' @param print.inline Logical specifying wheter data are printed in R notebook. Default is F - recommended if generated flexdashboard output.
#' @name m1.analysisParameters
#' @return list
#'
m1.analysisParameters <- function (which.data, which.strata = NULL, organism.filter.flag, organism.include, save.flag = T, save.filename, save.directory = "Preprocessed Datasets/",
                                   set_names = NULL, data.imputed.flag = F, vars2regress = "percent.mt", subsample_factor = 1, plt.log.flag = TRUE,
                                   RNA.upperlimit = 9000, RNA.lowerlimit = 200, mt.upperlimit = 60, cluster.resolution = 0.4, print.inline = FALSE){

  # Assertions
  if (is.null(which.data)) stop("Input not specified")
  if (is.null(which.strata)) which.strata <- NA
  if (is.null(set_names)) set_names <- which.data
  if (!is.logical(data.imputed.flag)) stop("data.imputed.flag must be logical")
  stopifnot(organism.include %in% c("Mm", "Hs"))
  if (!is.logical(print.inline)) stop("print.inline must be logical")
  stopifnot(is.logical(plt.log.flag) & length(plt.log.flag) == 1)
  if (RNA.upperlimit < 0 ) stop("RNA.upperlimit must be positive value")
  if (RNA.lowerlimit < 0 ) stop("RNA.lowerlimit must be positive value")
  if (mt.upperlimit < 0 | mt.upperlimit > 100) stop("mt.upperlimit must be value between 0 and 100")
  stopifnot(is.numeric(subsample_factor))
  stopifnot(subsample_factor <= 1 | subsample_factor >= 0)

  if (data.imputed.flag) RNA.upperlimit <- Inf

  # assign parameters
  analysis.parameters <- list(
    which.data = which.data,                          # specify input data files
    which.strata = which.strata,
    organism.filter.flag = organism.filter.flag,      # REQUIRED; logical
    organism.include = organism.include,              # character; options: "Hs", "Mm"
    save.flag = save.flag,                            # OPTIONAL; logical (default = T)  # save results
    save.filename = save.filename,                    # string; e.g., filename.Rdata
    save.directory = save.directory,
    set_names = set_names,
    data.imputed.flag = data.imputed.flag,
    vars2regress = vars2regress,
    subsample_factor = subsample_factor,
    plt.log.flag = plt.log.flag,
    RNA.upperlimit = RNA.upperlimit,                  # OPTIONAL; positive numerical (default = 9000)
    RNA.lowerlimit = RNA.lowerlimit,                  # OPTIONAL; positive numerical (default = 200)
    mt.upperlimit = mt.upperlimit,                    # OPTIONAL; positive numerical (default = 60)
    cluster.resolution = cluster.resolution,          # cluster.resolution
    print.inline = print.inline                      # print figures in R
  )

  return(analysis.parameters)
}


#' Load CellRanger preprocessed data
#'
#' Input data is output from cell ranger (10x datasets)
#'
#' @param import_set Character specifying folder containing matrix.mtx, gense.tsv and barcodes.tsv files provided by 10x.
#' @param subsample_factor Numeric [0,1]. Subsampling factor
#' @param all_input_organisms Character specifying species to include. This is necessary to parse ensemble ids into mouse and human. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item c("Hs", "Mm") - both species included
#' }
#' @param dir Character. folder containing import_set folder.
#' @name m1.loadCellRanger
#' @return list containing Seurat Object and named gene vector.
#'
m1.loadCellRanger <- function(import_set, subsample_factor, all_input_organisms, dir) {

  import_set_path <- paste(dir, import_set, sep ="")

  # import cell ranger-processed data
  expression_matrix <- Read10X(data.dir = import_set_path[1])

  # assign to secondary matrix
  expression_matrix2 <- expression_matrix

  # import gene and ensembl names
  feature.names = read.delim(paste(import_set_path[1], "/genes.tsv", sep = ""),
                             header = FALSE,
                             stringsAsFactors = FALSE)

  # remove hg19 tags
  feature.names[,1] <- gsub("hg19_", "", feature.names[,1])
  feature.names[,2] <- gsub("hg19_", "", feature.names[,2])

  # create gene list
  gNames <- as.character( feature.names$V2 );
  names(gNames) <- as.character( feature.names$V1 );
  names(gNames) <-gsub("\\..*","",as.vector( names(gNames)))

  # convert gene symbols to ensemble (in expression matrix)
  rownames(expression_matrix2) <- names(gNames)

  # assign species
  orgs <- scMiko::m1.inferSpecies(expression_matrix2, all_input_organisms)
  # if (length(all_input_organisms) > 1) {
  #   orgIDs <- rownames(expression_matrix2)[ apply( expression_matrix2,2,which.max )];
  #   orgs <- rep("Hs",ncol(expression_matrix2));
  #   orgs[grep("^ENSMUS",orgIDs) ] <- "Mm";
  #   names(orgs) <- colnames(expression_matrix2);
  # } else {
  #   orgs <- rep(all_input_organisms,ncol(expression_matrix2));
  #   names(orgs) <- colnames(expression_matrix2);
  # }


  # create seurat object
  # so = CreateSeuratObject(counts = expression_matrix2,project="cell_ranger",min.cells=3,min.features=200)
  so = CreateSeuratObject(counts = expression_matrix2,project="cell_ranger",min.cells=0,min.features=0)

  # add gene symbols as meta data in seurat object
  mat_ens <- rownames(so@assays[["RNA"]])
  match.id <- match(mat_ens, names(gNames))
  # gNames_filtered <- gNames[match.id]

  # Add gene symbols as meta data that we can use later
  # so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=gNames_filtered,col.name="gene_name");
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=names(gNames[match.id]),col.name="ENSEMBL");
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=as.vector(gNames[match.id]),col.name="SYMBOL");

  # Add in inferred organism
  so$Organism <- orgs;

  # Specify barcodes
  so$Barcode <- "unspecified";

  output <- list(so, gNames)

  return(output)
}


#' Load preprocessed data from Moffat lab sciRNA-seq3 pipeline
#'
#' Load preprocessed data from Moffat lab sciRNA-seq3 pipeline
#'
#' @param import_set Character vector specifying expression matrix (import_set[1]), rt barcodes (import_set[2]) and pcr barcodes (import_set[3]). Expression matrix will be imported successfully if barcodes are omitted.
#' @param subsample_factor Numeric [0,1]. Subsampling factor
#' @param all_input_organisms All species included in input files. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item c("Hs", "Mm") - both species included
#' }
#' @param organism.include Species to include in downstream analysis. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' \item c("Hs", "Mm") - both species included
#' }
#' @param dir Character. folder containing import_set files
#' @name m1.loadMoffat
#' @return list containing Seurat Object and named gene vector.
#'
m1.loadMoffat <- function(import_set, subsample_factor, all_input_organisms, organism.include, dir) {

  # load gene count matrix
  import_set_path <- paste(dir, import_set, sep ="")
  load(import_set_path[1])

  # assign row (genes) and column (cells) names
  gene_count2 <- gene_count

  # Infer organism
  orgs <- scMiko::m1.inferSpecies(gene_count2, all_input_organisms)
  df_cell$orgs <- orgs
  # if (length(all_input_organisms) > 1) {
  #   orgIDs <- rownames(gene_count2)[ apply( gene_count2,2,which.max )];
  #   orgs <- rep("Hs",ncol(gene_count2));
  #   orgs[grep("^ENSMUS",orgIDs) ] <- "Mm";
  #   names(orgs) <- colnames(gene_count2);
  # } else {
  #   orgs <- rep(all_input_organisms,ncol(gene_count2));
  #   names(orgs) <- colnames(gene_count2);
  # }


  # filter out incorrect species genes immediately.
  all.genes <- rownames(gene_count2)
  include.which.genes <- rep(FALSE, length(all.genes))
  if ((length(all_input_organisms) > length(organism.include)) | (length(organism.include) == 1)){
    for (i in 1:length(organism.include)){
      if (organism.include[i] == "Mm"){
        include.which.genes[grepl("MUSG", all.genes)] <- T
      } else if (organism.include[i] == "Hs"){
        include.which.genes[!(grepl("MUSG", all.genes))] <- T
      }
    }

    df_gene <- df_gene[include.which.genes, ]
    gene_count <- gene_count[include.which.genes, ]
    gene_count2 <- gene_count2[include.which.genes, ]
  }

  # only keep protein coding genes

  # assign gene and cell names
  rownames(gene_count2) <- df_gene$gene_id
  colnames(gene_count2) <- df_cell$sample


  # subsample data
  if (subsample_factor < 1){
    row_ind <- c(1:dim(gene_count2)[1])
    subsample_ind <- sample(row_ind, round(subsample_factor * dim(gene_count2)[1]), replace = FALSE)
    gene_count2 <- gene_count2[subsample_ind,]
  }

  # Need names vectors to add additional metadata to Seurat object
  gNames <- as.character( df_gene$gene_name );
  names(gNames) <- as.character( df_gene$gene_id );
  names(gNames) <-gsub("\\..*","",as.vector( names(gNames)))



  # Add PCR barcodes
  pcr.barcode.flag <- F
  try({
    grps <- read.csv(import_set_path[2], header=T,sep="\t",stringsAsFactors=F);
    wells <- gsub("Han_[0-9]+_([A-Z0-9]{2,3}).[ACGT]+$","\\1",colnames(gene_count2));
    myGrps <- grps$ConditionGroup[ match(wells,grps$sampleWell) ];
    names(myGrps) <- colnames(gene_count2)
    pcr.barcode.flag <- T
  }, silent = T)

  # Add RT barcodes
  rtBC <- read.csv(import_set_path[3],header=T,sep="\t",stringsAsFactors=F);

  barcodes <- gsub(".+([ACGT]{10}$)","\\1",colnames(gene_count2));
  idx <- match( barcodes,rtBC$rt.bc );
  rtGrp <- rtBC$Sample.type[idx];
  names(rtGrp) <- colnames(gene_count2);

  # remove .* suffix in ensembl
  rownames(gene_count2)<-gsub("\\..*","",as.vector( rownames(gene_count2)))

  # Create seurat object with count matrix data
  # so <- CreateSeuratObject(counts=gene_count2,project="Hong-sciSeq3",min.cells=3,min.features=200,names.field=2,names.delim="." )
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
  so[["Organism"]] <- so.merge$orgs



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

  # return output
  output <- list(so, gNames)

  return(output)
}



#' Load preprocessed TPM data (e.g., Neftel 2019 datasets)
#'
#' Load preprocessed TPM data (e.g., Neftel 2019 datasets)
#'
#' @param import_set Character specifying TPM matrix file name.
#' @param subsample_factor Numeric [0,1]. Subsampling factor
#' @param all_input_organisms All species included in input files. If multiple species are provided, error will thrown; only one expected. One of:
#' \itemize{
#' \item "Hs" - Human
#' \item "Mm" - Mouse
#' }
#' @param dir Character. folder containing import_set file
#' @name m1.loadTPM
#' @return list containing Seurat Object and named gene vector.
#'

m1.loadTPM <- function(import_set, subsample_factor, all_input_organisms, dir) {


  import_set_path <- paste(dir, import_set, sep ="")

  # import expression matrix data
  expression_matrix<-read.table(import_set_path[1],header=TRUE)

  # if (!exists("expression_matrix")) expression_matrix<-read.table(file = import_set[1], sep = '\t')

  # get gene and ensembl names
  feature.names <- as.vector(expression_matrix$GENE)

  # average duplicate rows (i.e., duplicate genes)
  if (length(unique(feature.names)) < length(feature.names)){
    expression_matrix.col <- expression_matrix
    expression_matrix.col <- dplyr::select(expression_matrix.col, -c("GENE"))
    expression_matrix.col$rowID <- seq(1, nrow(expression_matrix.col))
    expression_matrix.col.noDup <-  WGCNA::collapseRows(expression_matrix.col, rowGroup = expression_matrix$GENE, rowID = expression_matrix.col$rowID, method = "Average")
    new.mat <- as.data.frame(expression_matrix.col.noDup[["datETcollapsed"]])
    new.mat$GENE <- rownames(new.mat)
    expression_matrix2 <- new.mat
  } else {
    # assign to secondary matrix
    expression_matrix2 <- expression_matrix
  }

  feature.names <- as.vector(expression_matrix2$GENE)

  expression_matrix2 <- dplyr::select(expression_matrix2, -c("GENE"))
  rownames(expression_matrix2) <- feature.names
  stopifnot(length(all_input_organisms) == 1)
  g2eNames <- sym2ens(my.symbols =  feature.names, my.species = all_input_organisms)
  gNames <- g2eNames$SYMBOL
  names(gNames) <- g2eNames$ENSEMBL
  gNames <- gNames[!is.na(gNames)]
  names(gNames) <-gsub("\\..*","",as.vector( names(gNames)))

  # assign organism
  orgs <- rep(all_input_organisms,ncol(expression_matrix2));
  names(orgs) <- colnames(expression_matrix2);


  # create seurat object
  so = CreateSeuratObject(counts = expression_matrix2,project=import_set[1],min.cells=0,min.features=0)
  # so = CreateSeuratObject(counts = expression_matrix2,project=import_set[1],min.cells=3,min.features=200)

  # add gene symbols as meta data in seurat object
  mat_ens <- rownames(so@assays[["RNA"]])
  match.id <- match(mat_ens, gNames)
  gNames_filtered <- names(gNames)[match.id]

  # Add gene symbols as meta data that we can use later
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=gNames_filtered,col.name="ENSEMBL");
  so[["RNA"]] <- AddMetaData( object=so[["RNA"]],metadata=as.vector((gNames)[match.id]),col.name="SYMBOL");

  # Add in inferred organism
  so$Organism <- orgs;

  # Specify barcodes
  so$Barcode <- "unspecified";

  output <- list(so, gNames)

  return(output)
}



#' Assign barcodes labels to seurat metadata
#'
#' Assign barcode labels (in Seurat object) to "subset_group" metadata field, as specified by which.strata parameter.
#'
#' @param so Seurat Object
#' @param which.strata Barcode labeling parameter. If NA, "subset_group" metadata field is set to "pooled".
#' @name m1.barcodeLabels
#' @return Seurat Object (with updated metadata)
#'
m1.barcodeLabels <- function(so, which.strata) {
  # set Seurat subset labels for cells of interest
  if (!is.na(which.strata)){
    pattern <- paste(which.strata, collapse="|")
    so <- Seurat::SubsetData(object = so, cells = (grepl(pattern, so@meta.data[["Barcode"]])))

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
#' @param omit.na Logical specifying whether to omit NA entries (present when unfiltered 10x dataset is used). Default is True.
#' @name m1.getMitoContent
#' @return Seurat Object (with updated metadata)
#'
m1.getMitoContent <- function(so, gNames, omit.na = T) {
  hs <- grep("^ENSG",rownames(so[["RNA"]]),value=T)
  mm <- grep("^ENSMUSG",rownames(so[["RNA"]]),value=T)
  pctHS <- PercentageFeatureSet(so, features=hs)
  pctMM <- PercentageFeatureSet(so, features=mm)

  # find relevant entries (genes that were detected from master gene list)
  f.mm <- intersect(names(gNames)[ grep("^mt-",gNames) ],rownames(so[["RNA"]]))
  f.hs <- intersect(names(gNames)[ grep("^MT-",gNames) ],rownames(so[["RNA"]]))
  f <- intersect(names(gNames)[ grep("^(MT|mt)-",gNames) ],rownames(so[["RNA"]]))
  so[["percent.mt"]] <- PercentageFeatureSet(so, features=f)

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
#' @name m1.filterSeurat
#' @return List containing seurat object and filter summary statistics (as data.frame)
#'
m1.filterSeurat <- function(so, RNA.upperlimit = 9000, RNA.lowerlimit = 200, mt.upperlimit = 60, unmatch.low = 0, unmatch.high = 1, set_names = NULL) {

  # determine unfiltered UMI count
  original_count = length(so@meta.data[["nCount_RNA"]])

  # filter dataset
  if ((unmatch.low == 0) & (unmatch.high == 1)){
    so <- subset(so, subset = ((nFeature_RNA < RNA.upperlimit) & (nFeature_RNA > RNA.lowerlimit) & (percent.mt < mt.upperlimit)))
  } else {
    so <- subset(so,
                 subset = ((nFeature_RNA < RNA.upperlimit) & (nFeature_RNA > RNA.lowerlimit) & (percent.mt < mt.upperlimit) & (unmatched.rate > unmatch.low) & (unmatched.rate < unmatch.high)))
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

  output <- list(so, plt.filter_pre_post)
  return(output)

}


#' Normalize and Scale Data
#'
#' Applies one of two normalization/scaling approaches supported by Seurat.
#'
#' @param so Seurat Object
#' @param gNames Named vector of genes. Names are ensemble, entries are Symbols.
#' @param method Character specifying data normalization and scaling method. One of:
#' \itemize{
#' \item "NFS" - Seurat's NormalizeData, FindVariableFeatures, ScaleData workflow. Parameters are set to use LogNormalization method with a scale.factor of 1000. Variable features are selceted using 'mvp' method, and var2regress is regressed out during data scaling.
#' \item "SCT" - Default. SCTransform workflow; Uses regularized negative binomial regression to normalize UMI count data. All genes are returned (not only variable), and residual variace cutoff is set to 1.3; var2regress is included in model.
#' }
#' @param var2regress Character vector specifying which variables to regress out during data scaling.
#' @param enable.parallelization Logical specifying whether to enable parallelization. Default is T.
#' @param n.works Number of works to used during parallel processing. Default is 3.
#' @param max.memory Max memory to use during parallel processing. Default is 20480 * 1024^2
#' @name m1.scNormScale
#' @return Seurat Object
#'
m1.scNormScale <- function(so, gNames, method = "SCT", var2regress = NULL, enable.parallelization = T, n.workers = 3, max.memory = (20480 * 1024^2)){


  # enable parallelization
  # plan(strategy = "multisession", workers = parallel::detectCores())

  if (enable.parallelization){
    plan(strategy = "multisession", workers = n.workers)
    options(future.globals.maxSize = max.memory)
  }


  # Normalize and scale data
  if (method == "NFS"){

    # Normalize data
    so <- NormalizeData(so, normalization.method = "LogNormalize", scale.factor = 10000)
    # Find variable features
    so <- FindVariableFeatures(object = so, selection.method = 'mvp', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
    # Scale data
    so <- ScaleData(so, features = rownames(so), vars.to.regress = vars2regress)

  } else if (method == "SCT"){

    # apply sctransform (regularized negative binomial regression)
    # also removes confounding source of variation (i.e., mitochonrdial mapping percentage)
    so <- SCTransform(so,
                      vars.to.regress = vars2regress,
                      verbose = FALSE,
                      return.only.var.genes = FALSE,
                      variable.features.n = NULL,
                      variable.features.rv.th = 1.3)
  }

  return(so)
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
#' @name m1.inferSpecies
#' @return Character vector specifying species of each cell in expression  matrix.
#'
m1.inferSpecies <- function(exp.mat, expected.species, rep.ens.method = "alt"){

  # Infer organism
  if (length(expected.species) > 1) {

    # get representative Ensemble for each cell
    if (rep.ens.method == "orig"){
      orgIDs <- rownames(exp.mat)[ apply( exp.mat,2,which.max )]
    } else if (rep.ens.method == "alt"){
      n <- 1000
      nc <- ncol(expression_matrix2)
      z <- split(colnames(expression_matrix2), rep(1:ceiling(nc/n), each=n, length.out=nc))
      orgIDs <- character()
      for (i in names(z)) {
        chunk <- expression_matrix2[,z[[i]]]
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

