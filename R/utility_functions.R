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
#' @name getHeat
#' @return Heatmap object
#'
getHeat <- function(mat, hmcol, scale.limit){

  #threw this into a function to suppress plot that is automatically generated
  heat.object <- NULL
  try({
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      breaks =seq(-scale.limit, scale.limit, by = ((2*scale.limit/100))),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = "Gene Exp Matrix",
                                      xlab = "Cluster ID")

  }, silent = T)

  if (is.null(heat.object)){
    heat.object <- gplots::heatmap.2( mat, labCol= colnames(mat),
                                      trace="none",
                                      col= rev(hmcol),
                                      distfun = function(x) as.dist(1-cor(t(x))),
                                      hclustfun = function(x) hclust(x, method="average"),
                                      main = "Gene Exp Matrix",
                                      xlab = "Cluster ID")
  }

  return(heat.object)
}


#' Get list of available files
#'
#' Get list of files in specified directory. Default directory is "Preprocessed Datasets/".
#'
#' @param directory Character. Directory name.
#' @param scale.limit Numeric. Color scale limit.
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
#' Function that detachs and attached scMiko package.
#'
#' @name scMikoReload
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
#' \item "Reactome" - Default
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
#' @return Named list of vectors with gene sets.
#'
getAnnotationPathways <- function(query.genes, db = c("Reactome"), ontology = c("BP"), species = c("Hs")){

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
    which.match <- names(entrez2go.list) %in% all.genes.entrez

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
    pathways <- reactomePathways(all.genes.entrez)

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


