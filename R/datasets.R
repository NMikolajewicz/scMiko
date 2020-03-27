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
#' \item{coreCTL}{Core cytotoxic T lymphocyte kill genes. (Lawson et al. (2020) in revision.)}
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
#' }
#'
#' @keywords genesets, cell annotation, module scoring
#'
#' @examples
#' data(geneSets)
#' CancerSEA_Hs <- geneSets[["CancerSEA_Hs"]]
#'
"geneSets"
