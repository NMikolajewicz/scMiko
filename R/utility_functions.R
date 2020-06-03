#' Uppercase first letter and lowercase rest
#'
#' Function used to naively convert gene list to mouse representation; i.e., Convert Hs to Mm gene symbols.
#'
#'
#' @param x Character vector
#' @name firstup
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
#' @return Character specifying ensembl or symbol
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
#' @param reference.genes Named vector of genes; names are ENSEMBL, entries are SYMBOL.
#' @param query.genes vector of query genes to check representation
#' @name ens2sym.so
#' @return Seurat object with specified gene representation.
#'
ens2sym.so <- function(so, gNames.list, convert.RNA = TRUE){

    # var features
    so_temp <- so@assays[["SCT"]]@var.features
    so@assays[["SCT"]]@var.features <- as.vector((gNames.list[so_temp]))

    # scale data
    so_temp <- so@assays[["SCT"]]@scale.data
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@assays[["SCT"]]@scale.data <- so_temp

    # data
    so_temp <- so@assays[["SCT"]]@data
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@assays[["SCT"]]@data <- so_temp

    # metadata
    so_m <- so@assays[["SCT"]]@meta.features
    rownames(so_m) <-  make.names(as.vector((gNames.list[rownames(so_m)])), unique = T)
    so@assays[["SCT"]]@meta.features <- so_m

    # dimnames
    so_temp <- so@assays[["SCT"]]@counts@Dimnames[[1]]
    so@assays[["SCT"]]@counts@Dimnames[[1]] <- as.vector((gNames.list[so_temp]))

    # pca feature loading
    so_temp <-  so@reductions[["pca"]]@feature.loadings
    row.names(so_temp) <-  as.vector((gNames.list[row.names(so_temp)]))
    so@reductions[["pca"]]@feature.loadings <- so_temp

    # RNA ASSAY

    if (convert.RNA == TRUE){
      # var features
      so_temp <- so@assays[["RNA"]]@var.features
      so@assays[["RNA"]]@var.features <- as.vector((gNames.list[so_temp]))

      # data
      so_temp <- so@assays[["RNA"]]@data
      row.names(so_temp) <-   so@assays[["RNA"]]@meta.features[["gene_name"]]
      so@assays[["RNA"]]@data <- so_temp

      # counts
      so_temp <- so@assays[["RNA"]]@counts
      row.names(so_temp) <-   so@assays[["RNA"]]@meta.features[["gene_name"]]
      so@assays[["RNA"]]@counts <- so_temp

      # Dimnames
      so@assays[["RNA"]]@counts@Dimnames[[1]] <- so@assays[["RNA"]]@meta.features[["gene_name"]]
    }

    if (DefaultAssay(so) == "integrated"){

      # var features
      so_ens <- so@assays[["integrated"]]@var.features
      so@assays[["integrated"]]@var.features <-as.vector((gNames.list[so_ens]))

      # scale data
      so_sd <- so@assays[["integrated"]]@scale.data
      row.names(so_sd) <-  as.vector((gNames.list[row.names(so_sd)]))
      so@assays[["integrated"]]@scale.data <- so_sd

      # data
      so_d <- so@assays[["integrated"]]@data
      row.names(so_d) <-  as.vector((gNames.list[row.names(so_d)]))
      so@assays[["integrated"]]@data <- so_d
    }


    return(so)
}


#' Convert gene symbol representation to Hs or Mm
#'
#' Convert gene symbol to approriate species (as specified), while ensuring that gene is present in reference gene list.
#'
#' @param query.gene Character. Query gene symbol
#' @param gene.list Character(s). Vector of all genes under consideration
#' @param which.species Character. Species
#' @name speciesConvert
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


#' Create heatmap object
#'
#' Uses gplots::heatmap.2 function to generate heatmap object. Run in chunk specifying include = F to suppress plot and only output object.
#'
#' @param mat Matrix. Input matrix for heatmap.
#' @param hmcol Heatmap colors
#' @param scale.limit Numeric. Color scale limit.
#' @param main Character. Name of plot.
#' @param xlab Character. X axis label.
#' @param ... additional parameters passed to gplots::heatmap.2(...)
#' @name getHeat
#' @return Heatmap object
#'
getHeat <- function(mat, hmcol = NULL, scale.limit = NULL, main = NULL, xlab = NULL, ...){

  # get color palette
  if (is.null(hmcol)) hmcol <- colorRampPalette(RColorBrewer::brewer.pal(9, "RdBu"))(100)

  # get scale limits
  if (is.null(scale.limit)) scale.limit <- max(abs(mat))

  # get plot title
  if (is.null(main)) main <- "Gene Exp Matrix"

  # get x label
  if (is.null(xlab)) xlab <- "Cluster ID"

  #threw this into a function to suppress plot that is automatically generated
  heat.object <- NULL
  try({
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      breaks =seq(-scale.limit, scale.limit, by = ((2*scale.limit/100))),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = main,
                                      xlab = xlab,
                                      ...)

  }, silent = T)

  if (is.null(heat.object)){
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = "Gene Exp Matrix",
                                      xlab = "Cluster ID",
                                      ...)
  }

  return(heat.object)
}


#' Get list of available files
#'
#' Get list of files in specified directory. Default directory is "Preprocessed Datasets/".
#'
#' @param directory Character. Directory name.
#' @name getAvailableFiles
#' @return list of files
#'
getAvailableFiles <- function(directory = "Preprocessed Datasets/"){
return(list.files(directory))
}


#' Normalize and scale data within appropriate assay
#'
#' Ensures data are properly normalized and scaled. If intergrated dataset provided, Seurat's NormalizeData, FindVariableFeatures, ScaleData workflow is applied and default assay is set to 'RNA'. Otherwise, if  non-integrated dataset is provided, SCT transform is assumed and default assay is set to 'SCT'. If non-intergrated data have not been process with SCT workflow, see m1.scNormScale() function.
#'
#'
#' @param so Seurat Object
#' @name prepExpression
#' @return Seurat Object
#'
prepExpression <- function(so){

  if (DefaultAssay(so) == "integrated"){
    DefaultAssay(so) <- "RNA"
    so <-NormalizeData(so, verbose = FALSE)
    so <- ScaleData(so, verbose = FALSE)
    so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = 10000)
  } else {
    if ("SCT" %in% names(so@assays)) DefaultAssay(so) <- "SCT"
  }

  return(so)

}




#' Check if gene is avaialble in Seurat Object
#'
#' Return logical flag indicated whether query gene is present in seurat expression matrix.
#'
#' @param so Seurat Object
#' @param query.gene Character. Gene of interest.
#' @param reference.genes Named vector of genes; names are ENSEMBL, entries are SYMBOL.
#' @name isGeneAvailable
#' @return Logical
#'
isGeneAvailable <- function(so, query.gene, reference.genes){
  all.genes <- rownames(so@assays[[DefaultAssay(so)]])
  gene.rep <-  checkGeneRep(reference.genes, all.genes)
  stopifnot(gene.rep == "symbol")
  geneAvailable <- query.gene %in% all.genes

  return(geneAvailable)
}



#' reload scMiko
#'
#' Function that detachs and attaches scMiko package.
#'
#' @name scMikoReload
#' @examples
#'
#' # reload scMiko package
#' scMikoReload()
#'
scMikoReload <- function(){

  detach("package:scMiko", unload=TRUE)
  library(scMiko)

}



#' Merge list of seurat objects
#'
#' Merges list of seurat objects without any normalization of batch correction
#'
#' @param so.list List of seurat objects
#' @name mergeSeuratList
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
#' Set 'Seurat_Clusters' metadata entry to specified cluster resolution [0, inf]. See Seurat::FindClusters() for details.
#'
#'
#' @param so Seurat Object
#' @param cluster.resolution Numeric [0, inf] specifying cluster resolution. Values [0.1,1] typically perform well.
#' @name setResolution
#' @return Seurat object
#'
setResolution <- function (so, cluster.resolution){

  so <- FindClusters(object = so, resolution = cluster.resolution, verbose = 0, algorithm = 1)

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
subsetSeurat <- function (so, subset.df){

  # check if subset input is validd
  if (is.na(unique(subset.df$field))){
    subset.flag <- FALSE
  } else if ( unique(subset.df$field) %in% names(so@meta.data)) {
    subset.flag <- TRUE
  } else {
    subset.flag <- FALSE
  }


  # subset data
  if (subset.flag){
    pattern <- paste( "^", as.vector(subset.df$subgroups), "$", collapse="|")
    pattern <- gsub(" ", "", pattern)
    cur.field <- as.vector(unique(subset.df$field))
    so <- subset(x = so, cells = which(grepl(pattern, as.character(so@meta.data[[cur.field]]))))
  }

  return(so)
}


#' Returns list of annotations for given Entrez gene IDs
#'
#' Returns list of Reactome or GO annotations for given Entrez gene IDs
#'
#' @param query.genes Entrez IDs of query genes
#' @param db Database to retrieve annotations from. One of:
#' \itemize{
#' \item "Bader" - Default. List of pathway annotations curated by Bader lab (http://baderlab.org/GeneSets)
#' \item "Reactome"
#' \item "GO"
#' }
#' @param ontology GO ontologies to retrieve if GO db is selected. One of:
#' \itemize{
#' \item "BP" - Default. Biological processes
#' \item "MF" - Molecular functions
#' \item "CC" - Cellular components
#' }
#' @param species "Mm" or "Hs". Default is "Hs".
#' @name getAnnotationPathways
#' @return Named list of vectors with gene sets (Entrez format).
#'
getAnnotationPathways <- function(query.genes, db = c("Bader"), ontology = c("BP"), species = c("Hs")){

  if (db == "GO"){

    which.ontology <- ontology

    if (species == "Hs"){
      go.e2g <- org.Hs.egGO
      go.g2e <- as.list(org.Hs.egGO2EG)
    } else if (species == "Mm"){
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
#' @param which.assay Seurat assay to get data frome. Default is DefaultAssay(so).
#' @param which.data Specify which data to use (refers to slots in Seurat object assay). One of:
#' \itemize{
#' \item "scale" - Default
#' \item "data"
#' }
#' @param use.additional.genes Character vector of additional genes to include (in addition to varibale, if variable flag is specificed). Default is NA.
#' @name getExpressionMatrix
#' @return gene x cell expression matrix
#'
getExpressionMatrix <- function(so, only.variable = F, which.assay = NULL, which.data = "scale", use.additional.genes = NA){

  # specify assay
  if (is.null(which.assay)) which.assay <- DefaultAssay(so)

  # get complete matrix
  if (which.data == "scale"){
    exp.mat.complete <- so@assays[[which.assay]]@scale.data
  } else if (which.data == "data"){
    exp.mat.complete <- as.matrix(so@assays[[which.assay]]@data)
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
#' @return Character vector of cell barcodes/ids
#'
getExpressingCells <- function(so, query, expression.threshold = 0, which.data = "data"){
  exp.data.all<- getExpressionMatrix(so, which.data = which.data)
  exp.data.query <- exp.data.all[rownames(exp.data.all) %in% query, ]
  expressing.cells <- names(exp.data.query)[exp.data.query > expression.threshold]

  return(expressing.cells)
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
#' @return data.frame
#'
addLogEntry <- function(entry.name, entry, df.log, entry.variable.name = ""){

  df.log[nrow(df.log)+1, 1] <- as.character(entry.name)
  df.log[nrow(df.log), 2] <- as.character(entry.variable.name)
  df.log[nrow(df.log), 3] <- paste(entry, collapse=", ")

  return(df.log)
}


#' Create/update list of gene sets for scMiko package.
#'
#' Takes genessets stored in Excel sheets, and converts them to dataframes stored in lists.
#'
#' @param input.file Excel file (input). Must have ".xlsx" suffix. A character.
#' @param output.file Rdata file (output). Must have ".rda" suffix. A character.
#' @param dir Directory of input and output file (same folder). A character.
#' @param dev.directory.flag Logical indicating whether to use developer specific director. Default is False. If true dir is ignored.
#' @name updateGeneSets
#' @return List of data.frames saved as Rdata file.
#' @examples
#' \dontrun{
#' updateGeneSets("geneSets_MASTER_260420update.xlsx", "geneSets.rda", dev.directory.flag = T)
#' }
#'
updateGeneSets <- function(input.file, output.file, dir = "", dev.directory.flag = F){

  geneSets <- list()

  if (dev.directory.flag){
    input.dir <- "D:/Users/Nick/Dropbox/PDF Projects - JM/Data/scRNA-seq/01_sci-RNA-seq3_Hong_Kevin_Jason/Nick/scRNAseq/Reference Datasets/"
    output.dir <- "D:/Users/Nick/Dropbox/PDF Projects - JM/R Packages/scMiko/data/"
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
#' @param which.species Species
#' @name cleanFilterGenes
#' @return Character vector of genes
#'
cleanFilterGenes <- function(genes, so, which.species){

  # clean dataset and include only those available in seurat object
  cur.features <- genes
  cur.features <- cur.features[!is.na(cur.features)]
  # n.input <- length(cur.features)
  # cur.features <- lapply(cur.features,  gene.species.filter, rownames(so@assays[[DefaultAssay(so)]]@data), which.species)
  cur.features <- lapply(cur.features,  speciesConvert, rownames(so@assays[[DefaultAssay(so)]]@data), which.species)
  cur.features <-as.vector(unlist(cur.features))
  cur.features <- cur.features[!is.na(cur.features)]
  # n.output <- length(cur.features)

  return(cur.features)

}


#' Downsample single cell data
#'
#' Downsample data in Seurat object by specified factor
#'
#' @param so Seurat Object
#' @param subsample.factor Numeric [0,1]. Factor to downsample data by.
#' @name downsampleSeurat
#' @return Seurat Object
#'
downsampleSeurat <- function(so, subsample.factor){

  if (subsample.factor < 1){
    if (subsample.factor<0) stop("subsample.factor must be numeric between 0 and 1")
    n.subset <- round(subsample.factor *ncol(so))
    cell.ind <- sample(x = seq(1, ncol(so)), size = n.subset, replace = FALSE, prob = NULL)
    so <- SubsetData(so , cells = cell.ind)

  }
  return(so)

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
#' Get summary group expression in Seurat object. Can include mean, median, fraction (of expressing cells), sd, or cv.
#'
#' @param so Seurat Object
#' @param which.data Character specfying which data slot. Default is "data".
#' @param which.center Character indicating which summary measure to use. Must be one of "mean", "median", "fraction", "sd", or "cv". If unspecified, default is "mean".
#' @param which.group Character specfying group field in Seurat metadata. Default is "seurat_clusters".
#' @param do.parallel Logical specifying whether to perform computations in parallel. Default is F. Uses future.apply package.
#' @name avgGroupExpression
#' @return data.frame (gene rows, group columns)
#'
avgGroupExpression <- function(so, which.data = "data", which.center = "mean", which.group = "seurat_clusters", do.parallel = F){
  # which.center options: "mean", "fraction", "median", "sd", "cv"

  # inititate parallel processes
  if (do.parallel){
    library(future.apply)
    plan(multisession) ## Run in parallel on local computer
  }

  # entire matrix
  exp.mat.complete <- getExpressionMatrix(so, which.data = which.data)

  # group ID vector
  cluster.membership <- so@meta.data[[which.group]]

  # gene list
  gene.list <- rownames(exp.mat.complete)

  # ordered vector of unique groups
  u.clusters <- getOrderedGroups(so, which.group, is.number = F)

  if ((which.center == "fraction") & (which.data != "data")){
    warning("Data from 'data' slot used to compute expressing fraction")
    which.data <- "data"
    exp.mat.complete <- getExpressionMatrix(so, which.data = which.data)
    gene.list <- rownames(exp.mat.complete)
  }

  # clear some memory
  rm(so)
  gc()

  # compute measure of centrality
  avg.mat <- matrix(nrow = length(gene.list), ncol = length(u.clusters))
  for (i in 1:length(u.clusters)){
    if (which.center == "mean"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x))
      } else {
        avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x))
      }
    } else if (which.center == "median"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) median(x))
      } else {
        avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) median(x))
      }
    } else if (which.center == "sd"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x))
      } else {
        avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x))
      }
    } else if (which.center == "cv"){
      if (do.parallel){
        sd.cur <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x))
        av.cur <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x))
      } else {
        sd.cur <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sd(x))
        av.cur <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) mean(x))
      }
      avg.mat[,i] <- sd.cur / abs(av.cur)
    } else if (which.center == "fraction"){
      if (do.parallel){
        avg.mat[,i] <- future_apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sum(x>0)/length(x))
      } else {
        avg.mat[,i] <- apply(exp.mat.complete[ ,cluster.membership %in% u.clusters[i]], 1, function(x) sum(x>0)/length(x))
      }
    } else {
      stop("which.center must be specified as 'mean', 'median', 'fraction', 'sd', or 'cv'")
    }
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




#' Over-representation enrichment of GO terms using weighted fisher method
#'
#' Over-representation enrichment of GO terms using weighted fisher method
#'
#' @param gene.candidates Vector of candidate genes.
#' @param gene.universe Vector of all genes.
#' @param which.species Character specifying species. "Mm" or "Hs".
#' @param which.ontology Character specifying ontology. "BP", "MF", or "CC".
#' @param p.threshold p value threshold. Default is 0.001.
#' @param padj.threshold adjusted p value threhsold (BH). Default is 0.05.
#' @param topGO.object topGO object (Optional). If unspecified, new topGO object is new created. If specified, gene list is used to updated existing object.
#' @name enrichGO.fisher
#' @import topGo
#' @return list of 2 data.frames (unadjusted.results (p<0.001) and adjusted.results (padj < 0.05)) and a topGO object.
#'
enrichGO.fisher <- function(gene.candidates, gene.universe, which.species , which.ontology = "BP", p.threshold = 0.001,
                            padj.threshold = 0.05, topGO.object = NULL){

  if (!(which.species %in% c("Hs", "Mm"))) stop("Species incorrectly specified. Must be either Hs or Mm")

  if (which.species == "Hs"){
    library(org.Hs.eg.db)
    db <- "org.Hs.eg.db"
  } else if (which.species == "Mm"){
    library(org.Mm.eg.db)
    db <- "org.Mm.eg.db"
  }


  allGO2genes <- topGO::annFUN.org(whichOnto=which.ontology, feasibleGenes=NULL, mapping=db, ID="symbol")
  # all.genes <- colnames(datExpr)

  # make named factor showing which genes are of interest
  geneList=factor(as.integer(gene.universe %in% gene.candidates))
  names(geneList)= gene.universe

  if (is.null(topGO.object)){
    # if no GO object provided, create new
    GOdata <- new("topGOdata",
                  ontology=which.ontology,
                  allGenes=geneList,
                  annot=topGO::annFUN.GO2genes,
                  GO2genes=allGO2genes,
                  nodeSize=10)
  } else {
    # if GO object exists, update gene list
    GOdata <- topGO::updateGenes(topGO.object, geneList)
  }

  # define test using the weight01 algorithm (default) with fisher
  res.wfisher <- runTest(GOdata, algorithm='weight01', statistic='fisher')

  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO <- usedGO(GOdata)
  res.table <- GenTable(GOdata, weightFisher=res.wfisher, orderBy='weightFisher', topNodes=length(allGO))

  # ensure p values are numeric
  res.table$weightFisher <- as.numeric(res.table$weightFisher)

  #performing BH correction on our p values
  p.adj <- round(p.adjust(res.table$weightFisher,method="BH"),digits = 4)

  # create the file with all the statistics from GO analysis
  res.table.adj <- cbind(res.table,p.adj)
  res.table.adj=res.table.adj[order(res.table.adj$p.adj),]

  #get list of significant GO before multiple testing correction
  results.table.p <-  res.table.adj[which(res.table.adj$weightFisher<=p.threshold),]

  #get list of significant GO after multiple testing correction
  results.table.bh <- res.table.adj[which(res.table.adj$p.adj<=padj.threshold),]


  return(list(unadjusted.results = results.table.p, adjusted.results = results.table.bh, topGo.object = GOdata))
}




#' Filter seurat object according to specified meta data entries
#'
#' Using a mapping.list, entries from an existing meta data field are relabeled and mapped to new meta data field, and cells which are not specified in this mapping are omitted from the seurat object.
#'
#' @param so Seurat object.
#' @param old.field Existing seurat meta field that will be relabelled.
#' @param new.field New seurat meta field that will be created with relablled entries from old.field.
#' @param mapping.list Named list specifying how to relabel entries in old.field (mapping.list entries) to new.field (mapping.list names) in seurat meta data. Note that entries in mapping list do not have to match old.field entires exactly; entries are used as pattern arguemtn to grepl() function.
#' @name mapSubsetSeurat
#' @return Seurat object in which old field was mapped to new field, according to mapping specified in mapping.list. Any cells that are not mapped are omitted from seurat object.
#' @examples
#'
#' # define mapping lists
#' mapping.list.1 <- list(NSG = "NSG", BALBc = "BALB")
#' mapping.list.2 <- list(Early = "Early", Mid = "Mid", Late = "Late")
#'
#' # Create new field called "Mouse" using relabled entries from "Condition" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query, old.field = "Condition", new.field = "Mouse", mapping.list = mapping.list.1)
#'
#' # Create new field called "Time" using relabled entries from "Group" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query.2, old.field = "Group", new.field = "Time", mapping.list = mapping.list.2)
#'
mapSubsetSeurat <- function(so, old.field, new.field, mapping.list){

  keep.this <- NULL
  for (i in 1:length(mapping.list)){

    cur.group <- names(mapping.list)[order(names(mapping.list))][i]
    cur.pattern <- (mapping.list)[order(names(mapping.list))][[i]]

    if (is.null(keep.this)){
      keep.this <- grepl(cur.pattern, so@meta.data[[old.field]])
    } else {
      keep.this <- (keep.this | grepl(cur.pattern, so@meta.data[[old.field]]))
    }

    so@meta.data[[new.field]][grepl(cur.pattern, so@meta.data[[old.field]])] <- cur.group

  }


  # ensure groups are ordered and that order is maintained throughout analysis.
  u.groups <- as.character(unique(so@meta.data[[new.field]]))
  u.groups <- u.groups[order(u.groups)]
  so@meta.data[[new.field]] <- factor(so@meta.data[[new.field]], levels = u.groups)
  so <-so[, keep.this]

  return(so)
}


#' Perform differential expression analysis
#'
#' Give two existing fields within a Seurat Object, data are stratified by the first, and pairwise comparisons between groups in the second are performed to identify differential gene expression. Uses Seurat's FindMarkers() function.
#'
#' @param so Seurat object.
#' @param ordered.levels Vector of one or two seurat meta fields. If one provided, the second is set to "seurat_clusters".
#' @param which.assay Seurat assay to use. If unspecified, set to DefaultAssay(so).
#' @param ... Additional arguments passed to Seurat's FindMarkers() function.
#' @name multiLevel.FindMarkers
#' @return list of data.frames where each data.fram contains results from an individual pairwise comparison.
#' @examples
#'
#' # define mapping lists
#' mapping.list.1 <- list(NSG = "NSG", BALBc = "BALB")
#' mapping.list.2 <- list(Early = "Early", Mid = "Mid", Late = "Late")
#'
#' # Create new field called "Mouse" using relabled entries from "Condition" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query, old.field = "Condition", new.field = "Mouse", mapping.list = mapping.list.1)
#'
#' # Create new field called "Time" using relabled entries from "Group" - relabelling based on mapping list
#' so.query.2 <- mapSubsetSeurat(so.query.2, old.field = "Group", new.field = "Time", mapping.list = mapping.list.2)
#'
#' # Perform multilevel differential gene expression.
#' # Data are stratified by time, and then pairwise comparison between "Mouse" groups are performed.
#' # in this example, 3 total pairwise comparisons were made: BALBc vs NSG (Early), BALBc vs NSG (Mid) and BALBc vs NSG (Late)
#' deg.list <- multiLevel.FindMarkers(so.query.2, ordered.levels = c("Time", "Mouse"),
#'                                    logfc.threshold = 0,
#'                                    min.pct = 0,
#'                                    test.use = "wilcox")
#'
multiLevel.FindMarkers <- function(so, ordered.levels, which.assay = NULL, ...){

  if (length(ordered.levels) == 1) ordered.levels <- c(ordered.levels, "seurat_clusters")
  if (is.null(which.assay)) which.assay <- DefaultAssay(so)

  # get unique level factors
  level.1 <- unique(as.character(so@meta.data[[ordered.levels[1]]]))
  level.2 <-  unique(as.character(so@meta.data[[ordered.levels[2]]]))

  # set idents to level.2
  Idents(so) <- ordered.levels[2]

  # get all level 2 combinations
  level.2.combinations <- t(combinat::combn(level.2, 2))

  #N comparisons
  n.comparisons <- ncol(level.2.combinations) * length(level.1)

  #initiate results list
  deg.list <- list()
  for (i in 1:length(level.1)){

    # get current level 1 strata
    level.1.cur <- level.1[i]

    # subset seurat according to level 1
    keep.this <- as.character(so@meta.data[[ordered.levels[1]]]) %in% level.1.cur
    so.cur <-so[, keep.this]

    for (j in 1:nrow(level.2.combinations)){

      # define comparison level 2 comparison group
      group1 <- level.2.combinations[j, 1]
      group2 <- level.2.combinations[j, 2]

      comparison.label <- paste(level.1.cur, "|" , group1, ".vs.", group2, sep = "")

      suppressWarnings({
        deg.list[[comparison.label]] <- FindMarkers(so.cur,
                                                    assay = which.assay,
                                                    ident.1 = group1,
                                                    ident.2 = group2,
                                                    verbose = F,
                                                    ...)

      })

    }

  }

  return(deg.list)

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

  # compute connectivity
  # wi <- apply(w.mat, 1, function(x) sum(x))

  # store in data.frame
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
runWGCNA <- function(e.mat, s.mat = NULL, cor.metric = "rho_p", soft.power = 2, use.TOM = T, network.type = "signed", TOM.type = "unsigned", rescale.adjacency = F, ...){

  # similarity matrix - using proportionality metric for scRNAseq data.
  if (is.null(s.mat)){
    s.mat <-  dismay::dismay(e.mat, metric = cor.metric)
  }

  # adjacency matrix
  a.mat <-  sim2adj(s.mat, soft.power, network.type)

  # rescale value if needed
  if (rescale.adjacency) a.mat <- recaleValues(a.mat, new.min = 0, new.max = 1)

  # compute topological overlap matix (TOM)
  if (use.TOM){

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

    print2hide <-  capture.output(w.mat <- TOMsimilarity(a.mat.tom, TOMType = TOM.type, ...))
  } else {
    w.mat <- a.mat
  }

  # assign row and col names
  rownames(w.mat) <- rownames(a.mat)
  colnames(w.mat) <- colnames(a.mat)

  # dissimilarity measure
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
#' @param ... Additional arguments passessed to flashClust {flashClust package}
#' @name dist2hclust
#' @return hclust object
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

  stopifnot(class(d.mat) == "matrix")

  tree <- flashClust(as.dist(d.mat), method = method, ...)

  return(tree)
}


#' Cut clustered tree at varying heights and overlay with dendrogram to find optimal parameter set.
#'
#' WGCNA::cutreeDynamic is run for varying levels of deepSplit parameter (0:4), and cluster membership is assigned for each parameter set.
#'
#' @param tree h.clust object generated by dist2hclust.
#' @param d.mat distance matrix used to generate tree.
#' @param genes vector of gene names corresponding rows/col of distance matrix (d.mat). If specified, additional "genes" column is provided in output.
#' @param ... additional arguments passed to WGCNA::cutreeDynamic()
#' @name optimalDS
#' @return mColorh, matrix specfying module membership at varying deep split parameter specifications (0:4)
#' @examples
#'
#' geneTree <- dist2hclust(d.mat)
#' # determine number of modules based on refrence dataset
#' print2hide <- capture.output(mColorh <- optimalDS(tree = geneTree, d.mat = d.mat, genes  = rownames(a.mat)))
#'
optimalDS <- function(tree, d.mat, genes = NULL, ...){

  mColorh = NULL
  for (ds in 0:4){
    cut.tree <- cutreeDynamic(dendro = tree,distM= d.mat, cutHeight = 0.998, deepSplit=ds, pamRespectsDendro = FALSE, ...)
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



#' Wrapper to run  GO enrichment, using fisher enrichment method
#'
#' Wrapper to run  GO enrichment, using fisher enrichment method. Builds on TopGo functionality.
#'
#' @param gene.list named list of gene sets to query, where name specify name of gene set, and entry is vector of gene symbols.
#' @param gene.universe background genes
#' @param species Species: Hs or Mm
#' @name runEnrichment
#' @return list of data.frames. results.table.p contains unadjusted results, results.table.bh contains BH-adjusted results.
#' @import topGO
#' @examples
#'
#' enrichment.list <- runEnrichment(module.list, gene.universe = all.genes, species = "Hs")
#' results.p <- enrichment.list$results.table.p
#' results.bh <- enrichment.list$results.table.bh
#'
runEnrichment <- function(gene.list, gene.universe, species){


  #get list of significant GO before multiple testing correction
  results.table.p <- NULL

  #get list of significant GO after multiple testing correction
  results.table.bh <- NULL

  # enrich each module
  topGO.data <- NULL
  for (i in 1:length(gene.list)){

    if (species == "Hs"){

      enrich.list <- enrichGO.fisher(gene.candidates = toupper(gene.list[[i]]),
                                     gene.universe = toupper(gene.universe),
                                     which.species = species,
                                     p.threshold = 0.01,
                                     padj.threshold = 0.05,
                                     topGO.object = topGO.data)

    } else if (species == "Mm"){

      enrich.list <- enrichGO.fisher(gene.candidates = gene.list[[i]],
                                     gene.universe = gene.universe,
                                     which.species = species,
                                     p.threshold = 0.01,
                                     padj.threshold = 0.05,
                                     topGO.object = topGO.data)

    }

    # get results
    results.table.p.cur <- enrich.list$unadjusted.results
    results.table.bh.cur <- enrich.list$adjusted.results

    # get and reuse topGo object by updating gene list.
    topGO.data <- enrich.list$topGo.object


    if (nrow(results.table.p.cur) > 0) {
      results.table.p.cur$module <- names(gene.list)[i]
      results.table.p <- bind_rows(results.table.p, results.table.p.cur)
    }

    if (nrow(results.table.bh.cur) > 0) {
      results.table.bh.cur$module <-  names(gene.list)[i]
      results.table.bh <- bind_rows(results.table.bh, results.table.bh.cur)
    }
  }

  output <- list(
    results.table.p = results.table.p,
    results.table.bh = results.table.bh
  )

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
  # w.graph.top <- igraph::graph_from_data_frame(df.data, directed = F, vertices = NULL)
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
#' @examples
#'
#' # specify named list
#' my.list <- list(group.1 = c("a", "b", "c"), group.2 = c("d", "e", "f"))
#'
#' # convert named list to long data.frame
#' my.df <- namedList2longDF(my.list)
#'
namedList2longDF <- function(my.list, name.header = NULL, value.header = NULL){

  # unlist and assign values to data.frame
  my.vector <- unlist(my.list)
  my.df <- data.frame(as.vector(my.vector))

  # get data.frame column names
  if (is.null(value.header)) value.header <- "value"
  if (is.null(name.header)) name.header <- "name"
  colnames(my.df) <- value.header

  # assign names to corresponding values
  my.df[ ,name.header] <- NA
  for (i in 1:length(my.list)){
    my.df[  my.df[ ,value.header] %in% my.list[[i]] ,name.header] <- names(my.list)[i]
  }

  # rearrange column order
  my.df <- my.df[ , c(name.header, value.header)]

  # return long data.frame
  return(my.df)
}


#' Convert named list to wide data.frame
#'
#' Convert named list to wide data.frame. Resulting dataframe will have the same number of columns as there are names in the list.
#'
#' @param my.list named list
#' @name namedList2wideDF
#' @return wide data.frame
#' @examples
#'
#' # specify named list
#' my.list <- list(group.1 = c("a", "b", "c"), group.2 = c("d", "e", "f"))
#'
#' # convert named list to wide data.frame
#' my.df <- namedList2wideDF(my.list)
#'
namedList2wideDF <- function(my.list){

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
#' @return Numeric vector of rescaled values
#' @examples
#'
#' # rescale values
#' values <- rescaleValues(values)
#'
rescaleValues <- function(values, new.min = 0, new.max = 1){

  # set lower bound to zero
  old.min <- min(values)
  if (old.min < 0) {
    values <- values + abs(old.min)
  } else if (old.min > 0) {
    values <- values - abs(old.min)
  }
  stopifnot( min(values) == 0)

  # set upper bound to one
  old.max <- (max(values))
  values <- values/old.max
  stopifnot( max(values) == 1)

  new.range <- new.max - new.min
  values <- values * new.range
  values <- values + new.min

  stopifnot(min(values) == new.min)
  stopifnot(max(values) == new.max)

  return(values)

}



#' Bayesian Correlation algorithm
#'
#' Bayesian correlation scheme that assigns low similarity to genes that have low confidence expression estimates. Shown to be more reproducible tahn Pearson correlations. Source: https://www.biorxiv.org/content/10.1101/714824v1
#'
#' @param X Expression matrix
#' @name recaleValues
#' @return Correlation matrix
#' @examples
#'
#' #create a matrix (or load your own)
#' X <- matrix(1:1000, ncol=20)
#'
#' #compute the Bayesian correlation matrix
#' B <- BaCo(X)
#'
BaCo <- function(X){

  alpha0 <- rep(1/nrow(X),ncol(X))
  beta0=1-alpha0
  nrowsX <- nrow(X)
  k <- ncol(X)
  cs <- colSums(X)
  alphas <- alpha0 + X
  betas  <- matrix(rep(beta0,nrowsX), nrow=nrowsX, byrow=TRUE) + matrix(rep(cs,nrowsX), nrow=nrowsX, byrow=TRUE) - X
  alphasPLUSbetas <- alphas + betas
  Psi <- alphas/alphasPLUSbetas - matrix(rep(rowSums(alphas/alphasPLUSbetas)/k, k), ncol=k, byrow=FALSE)
  var_vec <- as.matrix( ( rowSums( (alphas*betas)/( (alphasPLUSbetas^2)*(alphasPLUSbetas+1) ) ) + rowSums(Psi^2) )/k )
  cov_mtrx <- (Psi %*% t(Psi))/k
  Bcorrvals <- cov_mtrx / sqrt( var_vec %*% t(var_vec) )
  diag(Bcorrvals) <- 1
  Bcorrvals
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
#' @name getSoftThreshold.v2
#' @return named list containing power estimates, r2 estimates, distribution plots, optimiation plot and results data.frame.
#' @examples
#'
#' # determine optimal soft threshold
#' sft <- getSoftThreshold.v2(s.mat, power =c(seq(0.5,5, by = 0.5), seq(6,10)),
#' network.type = "signed", rescale.adjacency = F)
#'
#' # visualize optimization plot
#' print(sft$optimization.plot)
#'
#' # visualize node-linkage density plots
#' cowplot::plot_grid(plotlist = sft$distribution.plot, ncol = 5)

getSoftThreshold.v2 <- function(s.mat, power =c(seq(0.5,5, by = 0.5), seq(6,10)), network.type = "signed", nBreaks = 20, removeFirst = T, rescale.adjacency = F){


  powers <- c(seq(0.5,5, by = 0.5), seq(6,10))
  plt.sf.list <- list()

  r2.sf <- c()

  for (i in 1:length(powers)){

    power.cur <- powers[i]
    a.cur <-  sim2adj(s.mat, soft.power = power.cur, network.type = network.type)

    if (rescale.adjacency)  a.cur <- recaleValues(a.cur, new.min = 0, new.max = 1)


    # get connectivity for specified power
    net.connectivity.df <- getConnectivity(a.cur, rownames(a.cur), flag.top.n = 20)

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
    # lm.summary <- summary(lm1)

    r2.sf[i] <- summary(lm1)[["r.squared"]]

    df.sf <- data.frame(x = log.dk, y = log.p.dk)

    # store node linkage distribution plot
    plt.sf.list[[as.character(power.cur)]] <- df.sf %>%
      ggplot(aes(x=x, y=y)) +
      geom_smooth(method = "lm", color = "tomato", fill = "tomato") +
      geom_point(size = 3) +
      xlab("N Links (Log)") +
      ylab("N Nodes (Log)") +
      ggtitle(paste0("Node Linkages\nSoft Power = ", power.cur)) +
      theme_classic() +
      stat_fit_glance(method = "lm",
                      label.y = "bottom",
                      method.args = list(formula = y ~ x),
                      mapping = aes(label = sprintf('r^2~"="~%.3f~~italic(P)~"="~%.2g',
                                                    stat(r.squared), stat(p.value))),
                      parse = TRUE)

  }

  # store powers and r2
  df.r2.sf <- data.frame(sf = powers, r2 = r2.sf)

  # optimizatio plot
  plt.opt.sf <- df.r2.sf %>%
    ggplot(aes(x = sf, y = r2)) +
    geom_smooth(color = "tomato", fill = "tomato") +
    geom_point(size = 3) +
    geom_hline(yintercept = df.r2.sf$r2[which.max(df.r2.sf$r2)]) +
    geom_vline(xintercept = df.r2.sf$sf[which.max(df.r2.sf$r2)]) +
    xlab("Soft Power") +
    ylab("R2 (Scale Free Topology)") +
    theme_classic() +
    ggtitle(paste0("Soft Power Optimization\nbest power = ", df.r2.sf$sf[which.max(df.r2.sf$r2)], ", r2 = ", signif(df.r2.sf$r2[which.max(df.r2.sf$r2)], 3)))

  # store results
  output <- list(
    powerEstimate = df.r2.sf$sf[which.max(df.r2.sf$r2)],
    r2Estimate = df.r2.sf$r2[which.max(df.r2.sf$r2)],
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
#' Remove duplicate genes from Seurat Object
#'
#' @param so Seurat Object
#' @name rmDuplicateGenes
#' @return seurat object
#' @examples
#'
#' so.query <- rmDuplicateGenes(so.query)
#'
rmDuplicateGenes <- function(so){

  # get gene names
  all.genes <- rownames(so)

  # if duplicates exist, subset seurat object
  if (sum(duplicated(all.genes)) > 0){
    # which.dup <- all.genes[duplicated(all.genes)]
    which.unique <- all.genes[!duplicated(all.genes)]
    so <- subset(so, features = which.unique)
  }

  return(so)

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
#' @name getJaccard
#' @return Jaccard similarity score (numeric)
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
getJaccard <- function(x1, x2){
  I <- length(intersect(x1, x2))
  S <- I/(length(x1) + length(x2) - I)
  return(S)
}

#' Jaccard Similarity Matrix
#'
#' Computes Jaccard similarity mtrix for list of genesets
#'
#' @param gene.sets named list of genesets, where names specify name of gene set, and entries are character vectors specifying genes belongs to the respective set.
#' @name jaccardSimilarityMatrix
#' @return Jaccard similarity matrix
#' @examples
#'
#' # compute jaccard similarity matrix for (named) list of genesets.
#' j.mat <- jaccardSimilarityMatrix(gene.sets)
#'
#' # generate heatmap
#' pheatmap::pheatmap(j.mat, show_colnames = F, main = "Jaccard Similarity")
#'
jaccardSimilarityMatrix <- function(gene.sets){
  n.sets <- length(gene.sets)
  j.mat <- matrix(nrow = n.sets, ncol = n.sets)
  for (i in 1:n.sets){
    for (j in 1:n.sets){
      i.name <- names(gene.sets)[i]
      j.name <- names(gene.sets)[j]
      j.mat[i, j] <- scMiko::getJaccard(gene.sets[[i.name]], gene.sets[[j.name]])
    }
  }

  rownames(j.mat) <- names(gene.sets)
  colnames(j.mat) <- names(gene.sets)

  return(j.mat)
}



