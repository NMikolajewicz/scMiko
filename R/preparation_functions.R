
#' Fix barcode labels
#'
#' Rename metadata entry 'CellTypes' to 'Barcode'. This is a fix implemented to correct an error from an earlier analysis pipeline. Analyses post January 2020 do not require this fix.
#'
#'
#' @param so Seurat Object
#' @name fixBarcodeLabel
#' @return Seurat object
#'
fixBarcodeLabel <- function (so){
  # merge CellType and Barcode, if necessary
  meta.data.names <- names(so@meta.data)

  if (("CellType" %in% meta.data.names) & ("Barcode" %in% meta.data.names)){
    if (DefaultAssay(so) == "integrated"){
      barcode <- so@meta.data[["Barcode"]]
      celltype <- so@meta.data[["CellType"]]
      barcode[is.na(barcode)] <- celltype[is.na(barcode)]
    } else {
      barcode <- so@meta.data[["CellType"]]
    }
  } else if (!("CellType" %in% meta.data.names) & ("Barcode" %in% meta.data.names)) {
    barcode <- so@meta.data[["Barcode"]]
  } else if (("CellType" %in% meta.data.names) & !("Barcode" %in% meta.data.names)) {
    barcode <- so@meta.data[["CellType"]]

  } else {stop("Problem with CellType/Barcode metadata detected. Troubleshooting required")}

  so@meta.data[["Barcode"]] <- barcode

  return(so)
}

#' Set cluster resolution
#'
#' Set 'Seurat_Clusters' metadata entry to quieried cluster resolution
#'
#'
#' @param so Seurat Object
#' @param cluster.resolution Seurat Object
#' @name setResolution
#' @return Seurat object
#'
setResolution <- function (so, cluster.resolution){

  so <- FindClusters(object = so, resolution = cluster.resolution, verbose = 0, algorithm = 1, modularity.fxn = 1)

  return(so)
}


#' prep Gene List
#'
#' Ensure is available and represented correctly.
#'
#' @param so Seurat Object
#' @param global.enviroment global.enviroment <- object()
#' @name prepGeneList
#' @return Named vector of available genes
#'
#'
prepGeneList <- function (so, global.enviroment){

  # global.enviroment <- object()


  if (("gNames.list_master" %in% global.enviroment)){
    gNames.list <- NULL
    for (i in 1:length(gNames.list_master)){
      gNames.list <- c(gNames.list, gNames.list_master[[i]] )
    }

    gNames.df <-  data.frame(n = gsub("\\..*","",as.vector(names(gNames.list))), g = as.vector(gNames.list))
    gNames.df <- unique(gNames.df)
    gNames.list <- as.vector(gNames.df$g)
    names(gNames.list) <- as.vector(gNames.df$n)
  }  else {

    # check if gene-ensemble pair are present in meta-data
    av.meta <- so@assays[["RNA"]]@meta.features

    if (all(c("SYMBOL", "ENSEMBL") %in% colnames(av.meta))){
      gNames.list <- as.vector(av.meta$SYMBOL)
      names(gNames.list) <- as.vector(av.meta$ENSEMBL)
    }

  }

  # ensure gene list is available
  stopifnot(exists("gNames.list"))

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

#' prep Seurat
#'
#' Preprocess using fixBarcodeLabel() and UpdateSeuratObject() Functions.
#'
#' @param so Seurat objects
#' @name prepSeurat
#' @return Seurat object
#'
prepSeurat <- function (so){

  so <- fixBarcodeLabel(so)
  so <- UpdateSeuratObject(so) # required after Seurat 3.1.2 update

  return(so)
}


#' Subset Seurat Object
#'
#' Subset Seurat object according to specific metadata field. Only specified metadata entries are retained, while remaining of data is omitted.
#'
#'
#' @param so Seurat object
#' @param subset.df Data.frame specifying which field to subset on, and which field entires to retain.
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
subsetSeurat <- function (so, subset.df){

  # check if subset input is validd
  if (is.na(subset.df$field)){
    subset.flag <- FALSE
  } else if (subset.df$field %in% names(so.query@meta.data)) {
    subset.flag <- TRUE
  } else {
    subset.flag <- FALSE
  }

  # subset data
  if (subset.flag){
    pattern <- paste( "^", as.vector(subset.df$subgroups), "$", collapse="|")
    pattern <- gsub(" ", "", pattern)
    cur.field <- as.vector(unique(subset.df$field))
    so.query <- SubsetData(object = so.query, cells = (grepl(pattern, so.query@meta.data[[cur.field]])))
  }

  return(so)
}


#' Initiate analysis log
#'
#' Create analysis log data.frame template
#'
#' @param module.name Character specifying module name
#' @name inititateLog
#' @return data.frame
#'
inititateLog <- function (module.name = ""){

  # Module
  df.log <- data.frame()
  df.log[nrow(df.log)+1, 1] <- as.character("Module")
  df.log[nrow(df.log), 2] <- as.character("")
  df.log[nrow(df.log), 3] <- as.character(module.name)
  colnames(df.log) <- c("Description", "Variable Name", "Value")

  # User
  df.log[nrow(df.log)+1, 1] <- as.character("User")
  df.log[nrow(df.log), 2] <- as.character("Sys.getenv('USERDOMAIN')")
  df.log[nrow(df.log), 3] <- as.character(Sys.getenv("USERDOMAIN"))

  # Date
  df.log[nrow(df.log)+1, 1] <- as.character("Date")
  df.log[nrow(df.log), 2] <- as.character("Sys.time()")
  df.log[nrow(df.log), 3] <- as.character(Sys.time())

  return(df.log)
}


#' Load Module-Specific Packages
#'
#' Load packages required for each analysis module
#'
#' @param module.number Numeric.
#' @param load.default Logical flag specifying wther to load default packages. If module.number is specified, load.default is ignored.
#' @name modulePackages
#' @return
#'
modulePackages <- function (module.number = NULL, load.default = T){

  if (!is.null(module.number)){


    if (module.number == 1){
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

