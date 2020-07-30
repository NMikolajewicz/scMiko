#' Load Module-Specific Packages
#'
#' Load packages required for each analysis module
#'
#' @param module.number Numeric.
#' @param load.default Logical flag specifying wther to load default packages. If module.number is specified, load.default is ignored.
#' @name modulePackages
#' @seealso \code{\link{library}}
#' @return
#'
modulePackages <- function (module.number = NULL, load.default = T){

  if (!is.null(module.number)){


    if (module.number == 1){

      load.default <- F

      # List of packages to load
      packages2load <- c("Seurat", "sctransform",
                         "plyr", "dplyr", "tidyr", "reshape2", "Matrix", "RColorBrewer", "ggplot2", "gridExtra",
                         "DT", "flexdashboard", "future")

      # load packages
      lapply(packages2load, library, character.only = TRUE)

    } else if (module.number == 2){
    } else if (module.number == 3){
    } else if (module.number == 4){
    } else if (module.number == 5){
    } else if (module.number == 6){
    } else if (module.number == 7){
    } else if (module.number == 8){
    } else if (module.number == 9){

      load.default <- F

      # List of packages to load
      # packrat::disable(project = getwd(), restart = TRUE)
      packages2load <- c("Seurat", "plyr", "dplyr", "tidyr", "reshape2", "RColorBrewer", "gridExtra",
                         "DT", "flexdashboard", "ggpmisc", "ggExtra", "grid", "ggrepel", "ddpcr",
                         "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db", "fgsea", "plotly", "ggplot2", "reactome.db")

      # load packages
      lapply(packages2load, library, character.only = TRUE)

      # Note: ensure that Plotly version 4.8.0 is used only
      #   require(devtools)
      #   install_version("plotly", version = "4.8.0", repos = "http://cran.us.r-project.org")
      session.info <- sessionInfo()
      stopifnot(session.info[["otherPkgs"]][["plotly"]][["Version"]] == "4.8.0")

    } else if (module.number == 10){

      load.default <- F

      # List of packages to load
      packages2load <- c("Seurat", "sctransform",
                         "plyr", "dplyr", "tidyr", "reshape2", "RColorBrewer", "ggplot2", "gridExtra",
                         "DT", "flexdashboard", "ggpmisc", "ggpmisc", "ggExtra", "grid", "ddpcr",  "future",
                         "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db", "scales")


      library("viridis")

      # load packages
      lapply(packages2load, library, character.only = TRUE)


    } else if (module.number == 11){
    } else if (module.number == 12){
    } else if (module.number == 13){
    } else if (module.number == 14){
    } else if (module.number == 15){
    } else if (module.number == 16){
    } else if (module.number == 17){
    } else if (module.number == 18){
    } else if (module.number == 19){
    } else if (module.number == 20){
    } else if (module.number == 21){
    }

  }

  if (load.default){
    # List of packages to load
    # packrat::disable(project = getwd(), restart = TRUE)
    packages2load <- c("Seurat", "plyr", "dplyr", "tidyr", "reshape2", "RColorBrewer", "gridExtra",
                       "DT", "flexdashboard", "ggpmisc", "ggExtra", "grid", "ggrepel", "ddpcr",
                       "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db", "fgsea", "plotly", "ggplot2", "reactome.db")

    # load packages
    lapply(packages2load, library, character.only = TRUE)

    # Note: ensure that Plotly version 4.8.0 is used only
    #   require(devtools)
    #   install_version("plotly", version = "4.8.0", repos = "http://cran.us.r-project.org")
    session.info <- sessionInfo()
    stopifnot(session.info[["otherPkgs"]][["plotly"]][["Version"]] == "4.8.0")
  }

}

#' Initiate analysis log
#'
#' Create analysis log data.frame template
#'
#' @param module.name Character specifying module name
#' @name initiateLog
#' @seealso \code{\link{addLogEntry}}
#' @return data.frame
#' @examples
#'
#' # intiate data log
#' df.log <- initiateLog("log name")
#'
#' # add new entry to log
#' df.log <- addLogEntry("Query File (.Rdata)", (input.file), df.log, "input.file")
#'
initiateLog <- function (module.name = ""){

  # Module
  df.log <- data.frame()
  df.log[nrow(df.log)+1, 1] <- as.character("Module")
  df.log[nrow(df.log), 2] <- as.character("")
  df.log[nrow(df.log), 3] <- as.character(module.name)
  colnames(df.log) <- c("Description", "Variable Name", "Value")

  # Date
  df.log[nrow(df.log)+1, 1] <- as.character("Date")
  df.log[nrow(df.log), 2] <- as.character("Sys.time()")
  df.log[nrow(df.log), 3] <- as.character(Sys.time())

  return(df.log)
}


#' Add entry to analysis log
#'
#' Add entry to analysis log
#'
#' @param entry.name name of entry. A character.
#' @param entry entry.
#' @param df.log existing analysis log to add entry to. A data.frame.
#' @param entry.variable.name name of variable storing entry (optional). A character.
#' @name addLogEntry
#' @seealso \code{\link{initiateLog}}
#' @author Nicholas Mikolajewicz
#' @return datalog data.frame
#' @examples
#'
#' # intiate data log
#' df.log <- initiateLog("log name")
#'
#' # add new entry to log
#' df.log <- addLogEntry("Query File (.Rdata)", (input.file), df.log, "input.file")
#'
addLogEntry <- function(entry.name, entry, df.log, entry.variable.name = ""){

  df.log[nrow(df.log)+1, 1] <- as.character(entry.name)
  df.log[nrow(df.log), 2] <- as.character(entry.variable.name)
  df.log[nrow(df.log), 3] <- paste(entry, collapse=", ")

  return(df.log)
}



