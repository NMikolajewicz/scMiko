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
