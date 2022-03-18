#' List of scMiko gene sets
#'
#' Reference gene sets consolidated from numerous sources to faciliate cell annotation and module scoring.
#'
#' @docType data
#' @name geneSets
#'
#' @usage data(geneSets)
#'
#' @format A list of data.frames containing the following genesets:
#' \describe{
#' \item{CancerSEA_Hs (\href{http://biocc.hrbmu.edu.cn/CancerSEA/home.jsp}{link})}{Cancer single-cell state atlas}
#' \item{CellMarker_Hs_Zhang2019 (\href{http://biocc.hrbmu.edu.cn/CellMarker/}{link})}{Cell marker atlas: Manually-curated celltype-specific gene sets from over 100,000 published papers. (Zhang et al. (2019) Nucleic Acids Res.)}
#' \item{coreCTL}{Core cytotoxic T lymphocyte kill genes. (\href{https://www.nature.com/articles/s41586-020-2746-2}{Lawson et al. (2020)}}
#' \item{cyclingGenes_Hs_Seurat (\href{https://science.sciencemag.org/content/352/6282/189}{link})}{Cell-cycle markers (Tirosh et al. (2016) Science)}
#' \item{GBM_Hs_Neftel2019 (\href{https://www.ncbi.nlm.nih.gov/pubmed/31327527}{link})}{Glioblastoma-subtype specific markers (Neftel et al. (2019) Cell)}
#' \item{Immune_Hs_Nirmal2018 (\href{https://www.ncbi.nlm.nih.gov/pubmed/30266715}{link})}{Immune cell signatures for 7 immune cell types identified in tumors using scRNAseq (Nirmal et el. (2018) Cancer Immunol Res.)}
#' \item{MCA_Mm_top50 (\href{https://www.ncbi.nlm.nih.gov/pubmed/30787437/}{link})}{Top 50 differentially-markers for each cell in Mouse organogenesis cell atlas characterized using scRNAseq. (Cao et al. (2019) Nature)}
#' \item{MCA_Mm_top100 (\href{https://www.ncbi.nlm.nih.gov/pubmed/30787437/}{link})}{Top 100 differentially-markers for each cell in Mouse organogenesis cell atlas characterized using scRNAseq. (Cao et al. (2019) Nature)}
#' \item{neuralDif_Mm_HH }{Neural differential markers provided by Hong Han.}
#' \item{Panglao_Mm (\href{https://panglaodb.se/}{link})}{Murine cell-type specific markers consolidated from multiple studies through unified framework. (Franzen et al. (2019) Database (Oxford))}
#' \item{Panglao_Hs (\href{https://panglaodb.se/}{link})}{Human cell-type specific markers consolidated from multiple studies through unified framework. (Franzen et al. (2019) Database (Oxford))}
#' \item{Renal_Mm_PMID29622724 (\href{https://www.ncbi.nlm.nih.gov/pubmed/29622724}{link})}{Murine kidney cell-specific markers (Park et al. (2018) Science)}
#' \item{TRRUSTv2_Mm (\href{https://www.grnpedia.org/trrust/}{link})}{Murine-specific trancription factor genesets; TRRUSTv2 database. (Han et al. (2018) Nucleic Acid Res.)}
#' \item{TRRUSTv2_Hs (\href{https://www.grnpedia.org/trrust/}{link})}{Human-specific trancription factor genesets; TRRUSTv2 database. (Han et al. (2018) Nucleic Acid Res.)}
#' \item{TRRUSTv2_regulation_Mm (\href{https://www.grnpedia.org/trrust/}{link})}{Murine-specific trancription factor genesets, with regulatory annotation (i.e., activating vs repressing); TRRUSTv2 database. (Han et al. (2018) Nucleic Acid Res.)}
#' \item{TRRUSTv2_regulation_Hs (\href{https://www.grnpedia.org/trrust/}{link})}{Human-specific trancription factor genesets, with regulatory annotation (i.e., activating vs repressing); TRRUSTv2 database. (Han et al. (2018) Nucleic Acid Res.)}
#' \item{VastDB_Mm (\href{http://vastdb.crg.eu/wiki/Main_Page}{link})}{Bulk RNAseq-derived cell-specific markers (Tapial et al. (2017) Genome Res.)}
#' \item{Thorsson_immuneSubtypes (\href{https://gdc.cancer.gov/about-data/publications/panimmune}{link})}{Tumor immunesubtype markers (Thorsson et al. (2018) Immunity)}
#' \item{immune_literaturePanel}{Immune marker panels obtained from literature review}
#' \item{universal_literaturePanel}{universal marker panels obtained from literature review}
#' \item{cancer_literaturePanel}{cancer marker panels obtained from literature review}
#' \item{AnimalTFDB (\href{http://bioinfo.life.hust.edu.cn/AnimalTFDB/#!/}{link})}{Animal transcription factor database. List of human and murine transcription factors}
#' \item{Richards_NatureCancer_2021_sc (\href{https://www.nature.com/articles/s43018-020-00154-9#Sec47/}{link})}{Glioblastoma stem cell injury response and developmental genesets. Single cell derived. }
#' \item{Richards_NatureCancer_2021_bulk (\href{https://www.nature.com/articles/s43018-020-00154-9#Sec47/}{link})}{Glioblastoma stem cell injury response and developmental genesets. Bulk RNA derived. }
#' \item{Cahoy_JNeurosci}{Neural cell gene panels}
#' \item{Liddelow_Nature_2017}{Astrocyte gene panels, specifically A1 and A2 reactive astrocyte markers}
#' \item{Zhong_Nature_2018}{Neural cell gene panels}
#' \item{Verhaak_CancerCell_2010}{Verhaak GBM-subtype genesets}
#' \item{Combes2022_Immune_TableS4 (\href{https://pubmed.ncbi.nlm.nih.gov/34963056/}{link})}{Immune-related gene markers reported in Table S4 of Combes 2022}
#' \item{Combes2022_Immune_TableS5 (\href{https://pubmed.ncbi.nlm.nih.gov/34963056/}{link})}{Immune-related gene markers reported in Table S5 of Combes 2022}
#' \item{Cell_Catalog (\href{https://nmikolajewicz.github.io/scMiko/articles/Marker_Catalog.html}{link})}{Cell-type markers derived from public scRNAseq atlases.}
#' \item{CancerSCEM (\href{https://ngdc.cncb.ac.cn/cancerscem/documents/#cellTypeDisL}{link})}{Cell-type specific marker genes used by CancerSCEM.}
#' \item{Stem_Catalog}{Stemness markers, literature curated. Sources include (\href{http://193.136.227.155/stemchecker/stemsethome.jsf}{StemChecker}, \href{https://www.sciencedirect.com/science/article/pii/S0092867418303581#mmc1}{Malta}, \href{https://bmccancer.biomedcentral.com/articles/10.1186/s12885-021-08351-0}{Hong}, \href{https://www.nature.com/articles/nature07056}{Mikkelsen}, \href{https://link.springer.com/article/10.1186/gb-2012-13-8-r71#ref-CR32}{Palmer}, \href{https://ashpublications.org/blood/article/103/8/2956/18022/Gene-expression-in-human-embryonic-stem-cell-lines}{Bhattacharya}, \href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6313194/}{Dong}, \href{https://www.thelancet.com/journals/ebiom/article/PIIS2352-3964(19)30114-8/fulltext#supplementaryMaterial}{Pece}, \href{https://translational-medicine.biomedcentral.com/articles/10.1186/s12967-020-02527-1}{Feng}, \href{https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215452#pone.0215452.s002}{Schulten})}
#' \item{HALLMARK (\href{https://www.cell.com/fulltext/S2405-4712(15)00218-5}{link})}{50 “hallmark” gene sets from the Molecular Signature Database (MSigDB)}
#' }
#'
#' @keywords genesets, cell annotation, module scoring
#'
#' @examples
#' data(geneSets)
#'
#' # retrieve geneset
#' cell_catalog <- geneSets[["Cell_Catalog"]]
#'
#' # convert wide dataframe format to named list format
#' call_catalog.list <- wideDF2namedList(cell_catalog)
#'
"geneSets"



#' Pathway annotations from Bader Lab
#'
#' Updated list of gene sets used for enrichment analysis. Consolidated by Bader Lab, includes KEGG, Msigdb, NCI, IOB, NetPath, HumanCyc, Reactome, GO, Panther. http://baderlab.org/GeneSets
#'
#' @docType data
#' @name baderPathways
#'
#' @usage {
#' data(baderPathways)
#'
#'
#' # To update list:
#'
#' # specify location of input data
#' # directory <- "D:/Users/Nick/Dropbox/PDF Projects - JM/Data/CRISPR screen/GIN Data/"
#' # file.pathways <- list(
#' #   Hs.entrez = "Human_GOBP_AllPathways_no_GO_iea_April_01_2020_entrezgene.gmt",
#' #   Hs.symbol = "Human_GOBP_AllPathways_no_GO_iea_April_01_2020_symbol.gmt",
#' #   Mm.entrez = "Mouse_GOBP_AllPathways_no_GO_iea_April_01_2020_entrezgene.gmt",
#' #   Mm.symbol = "Mouse_GOBP_AllPathways_no_GO_iea_April_01_2020_symbol.gmt"
#' # )
#' #
#' # pathway.list <- list()
#' # for (i in 1:length(file.pathways)){
#' #
#' #   df.pathway.2 <- qusage::read.gmt(paste0(directory, file.pathways[[names(file.pathways)[i]]]))
#' #   pathways <- df.pathway.2
#' #   names(pathways) <- gsub("%.*", "", names(pathways))
#' #
#' #   pathway.list[[names(file.pathways)[i]]] <- pathways
#' #
#' # }
#' #
#' # baderPathways <- pathway.list
#' #
#' # output.dir <- "D:/Users/Nick/Dropbox/PDF Projects - JM/R Packages/scMiko/data/"
#' #
#' # save(baderPathways, file=getLoadPath("baderPathways.rda", output.dir))
#' }
#'
#' @format A list of pathways
#'
#' @keywords genesets, pathway annotation
#'
#' @examples
#' data(baderPathways)
#'
"baderPathways"


#' Ligand-Receptor Database
#'
#' List of ligand-receptor pairs consolidated from iTalk, SCSR, Phantom5, NATMI and Moffat Lab (in house).
#'
#' @docType data
#' @name LR.db
#' @author Magali Aguilera Uribe
#' @usage {
#' data(LR.db)
#'
#'
#' # To update list:
#'
#' # import data
#' load("C:/Users/Owner/Dropbox/PDF Projects - JM/Data/scRNA-seq/01_sci-RNA-seq3_Hong_Kevin_Jason/NM_HH/Data/Gene_Sets/LR.db_Magali_301120.RData")
#'
#' # parse databases
#' list.source <- strsplit(LR.db$Which_db, " ") # LR.db$Which_db
#' names(list.source) <- LR.db$Pair
#'
#' all.db <- unique(unlist(list.source))
#' db.list <- list()
#' for (i in 1:length(all.db)){
#'   db.list[[all.db[i]]] <- names(list.source)[unlist(lapply(list.source, function(x) all.db[i] %in% x))]
#' }
#'
#' # upset plot
#' plt.upset <- upset.Plot(gene.sets = db.list, row.title = "", column.title = "")
#'
#' LR.df <- LR.db
#' LR.list <- db.list
#' LR.upset.plot <- plt.upset
#' LR.db <- list(
#'   LR.df = LR.df,
#'   LR.list = LR.list,
#'   LR.upset.plot = LR.upset.plot
#' )
#'
#' output.dir <- "C:/Users/Owner/Dropbox/PDF Projects - JM/R Packages/scMiko/data/"
#' save(LR.db, file=getLoadPath("LR.db.rda", output.dir))
#' }
#'
#' @format A list containing LR database
#'
#' @keywords ligand-receptor pairs, communication network
#'
#' @examples
#' data(LR.db)
#'
"LR.db"
