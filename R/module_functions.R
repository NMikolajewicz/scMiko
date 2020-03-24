#' M9. Import Gene List
#'
#' Import gene lists for subsequent analysis in Module 9.
#'
#' @param query.format Numerical code specifying gene list imput format. One of:
#' \itemize{
#' \item 1 - Character vector of genes. E.g., tnf.signalling = c("Tnfrsf1a", "Tnfrsf1b", "Tnfaip3")
#' \item 2- csv file input.
#' }
#' @param input.file Character specifying module name. Required if query.format = 2.
#' @param markers.of.interest Character vector specifying list of gene(s). Required if query.format = 1.
#' @param which.col Character specifying which columsn to take gene lists from. Required if query.format = 2.
#' \itemize{
#' \item "first" - Include only first csv column.
#' \item "all" - Include all csv columns
#' }
#' @param which.species Character specifying species.
#' @name m9.importMarkers
#' @return List
#'
m9.importMarkers <- function(query.format, input.file = NULL, markers.of.interest = NULL, which.col = "first", which.species) {

  stopifnot(which.col %in% c("first", "all"))

  if (query.format == 1){

    if (which.species == "Hs"){
      markers.of.interest <-  toupper(markers.of.interest)
    } else {
      markers.of.interest <-  firstup(markers.of.interest)
    }
    marker_set_name <- ""

  } else if (query.format == 2){
    spreadsheet_format <- unlist(strsplit(input.file, '[.]'))[2]

    if (spreadsheet_format == "csv") {
      markers.of.interest <- read.csv(input.file, header = TRUE)
    } else if (spreadsheet_format == "xlsx") {
      markers.of.interest <- readxl::read.xlsx(input.file, header = TRUE)
    }
    marker_set_name <- colnames(markers.of.interest)
    if (which.col == "first"){
      spread_sheet_col <- 1
      marker_set_name <- marker_set_name[spread_sheet_col]
      markers.of.interest <- as.data.frame(as.vector(markers.of.interest[ , spread_sheet_col]))
      colnames(markers_of_interest) <- marker_set_name
    } else if (which.col == "all"){
      marker_set_name <- ""
    }

    for (i in 1:dim(markers.of.interest)[2]){
      if (which.species == "Hs"){
        markers.of.interest[ ,i] <-  toupper(as.vector(markers.of.interest[ ,i]))
      } else {
        markers.of.interest[ ,i] <-  firstup(as.vector(markers.of.interest[ ,i]))
      }

    }

  }

  output <- list(markers.of.interest, marker_set_name)
  names(output) <- c("markers.of.interest", "marker_set_name")
  return(output)
}



#' M23. Stratify gene expression into high and low expression groups
#'
#' Creates new metadata field stratifying gene expression (for queried gene) into high and low expresion groups. If return.plt is true, histogram illustrating threshold is generated; Otherwise seurat object with new metadata entry is returned.
#'
#' @param so Seurat Object
#' @param query.gene Character specifying gene of interest.
#' @param quartile.threshold Quantile used to threshold expression data. E.g., 0.25 means that <25th percentile and >75th percentile are flagged as low and high expression groups, respecitvely.
#' @param return.plt Logical indicating whether histogram with thresholds is returned. Otherwise Seurat object is returned.
#' @name m23.binarizeExpression
#' @return Seurat object (unless return.plt is T, in which case ggplot handle)
#'
m23.binarizeExpression <- function(so, query.gene, quantile.threshold = 0.25, return.plt = F){


  which.match <- which(rownames(so@assays[[DefaultAssay(so)]]) %in% query.gene)

  df.exp <- data.frame(data = so@assays[[DefaultAssay(so)]]@data[which.match, ],
                       scale = so@assays[[DefaultAssay(so)]]@scale.data[which.match, ])

  df.exp$cellID <- seq(1, ncol(so@assays[[DefaultAssay(so)]]@data))

  df.exp.filtered <- df.exp[df.exp$data > median(df.exp$data), ]

  df.exp.filtered$rank <- percent_rank(df.exp.filtered$data)

  df.exp.filtered$expression.grp <- "other"
  df.exp.filtered$expression.grp[df.exp.filtered$rank > (1-quantile.threshold)] <- "high"
  df.exp.filtered$expression.grp[df.exp.filtered$rank < (quantile.threshold)] <- "low"

  expression.group <- paste(query.gene, ".Expression", sep  = "")
  so@meta.data[[expression.group]] <- "other"
  so@meta.data[[expression.group]][df.exp.filtered$cellID[df.exp.filtered$expression.grp == "high"]] <- "high"
  so@meta.data[[expression.group]][df.exp.filtered$cellID[df.exp.filtered$expression.grp == "low"]] <- "low"


  if (return.plt){
    plt.hist.filtered <-   df.exp.filtered %>%
      ggplot(aes(x = data)) +
      geom_histogram(aes(fill =  df.exp.filtered$expression.grp)) +
      xlab("Expression") + ggtitle(paste(query.gene, "Expression")) +
      scale_fill_manual(values = c( "red", 'blue', 'grey') )  +
      theme_classic()  + labs(fill = "Expression")

    return(plt.hist.filtered)
  } else {
    return(so)
  }


}


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
#' \item c("Hs", "Mm") - both species included
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




