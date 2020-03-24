#' M1 Load Module 1-Specific Parameter Specifications
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
  if (length(all_input_organisms) > 1) {
    orgIDs <- rownames(expression_matrix2)[ apply( expression_matrix2,2,which.max )];
    orgs <- rep("Hs",ncol(expression_matrix2));
    orgs[grep("^ENSMUS",orgIDs) ] <- "Mm";
    names(orgs) <- colnames(expression_matrix2);
  } else {
    orgs <- rep(all_input_organisms,ncol(expression_matrix2));
    names(orgs) <- colnames(expression_matrix2);
  }

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
  if (length(all_input_organisms) > 1) {
    orgIDs <- rownames(gene_count2)[ apply( gene_count2,2,which.max )];
    orgs <- rep("Hs",ncol(gene_count2));
    orgs[grep("^ENSMUS",orgIDs) ] <- "Mm";
    names(orgs) <- colnames(gene_count2);
  } else {
    orgs <- rep(all_input_organisms,ncol(gene_count2));
    names(orgs) <- colnames(gene_count2);
  }
  df_cell$orgs <- orgs

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

  # assign to secondary matrix
  expression_matrix2 <- expression_matrix

  # get gene and ensembl names
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





