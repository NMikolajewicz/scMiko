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





