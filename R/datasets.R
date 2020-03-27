#' List of scMiko gene sets
#'
#' Reference gene sets consolidated from numerous sources to faciliate cell annotation and module scoring.
#'
#' @docType data
#' @name geneSets
#'
#' @usage data(geneSets)
#'
#' @format A list of data.frames - each data.frame is an independent dataset.
#'
#' @keywords genesets, cell annotation, module scoring
#'
#' @references Moore et al. (2013) Genetics 195:1077-1086
#' (\href{https://www.ncbi.nlm.nih.gov/pubmed/23979570}{PubMed})
#'
#' @source \href{https://phenome.jax.org/projects/Moore1b}{QTL Archive}
#'
#' @examples
#' data(geneSets)
#' CancerSEA_Hs <- geneSets[["CancerSEA_Hs"]]
#'
"geneSets"
