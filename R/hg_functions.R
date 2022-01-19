
#' Summarize hypergeometric enrichment results
#'
#' Summarize hypergeometric enrichment results. Collapses HG result list into dataframe, and generates plots to visualize enrichments for each geneset.
#'
#' @param hg.res RunHG output list
#' @param fdr.filter FDR filter threshold. Default 1 (i.e., include all)
#' @param do.plot Logical to generate enrichment plots for each geneset.
#' @param show.n Number of top enrichments to visualize in plots. Ignored if do.plot == F.
#' @param genesets Geneset list used for HG enrichment (i.e., RunHG input). Optional.
#' @param str.wrap.width positive integer giving target line width in characters. A width less than or equal to 1 will put each word on its own line. Ignored if do.plot == F.
#' @param plt.title Character specifying plot title. Ignored if do.plot == F.
#' @param pathway.name.size Font size for pathway names in plots. Ignored if do.plot == F.
#' @name summarizeHG
#' @concept enrichment
#' @author Nicholas Mikolajewicz
#' @return list of summarized results
summarizeHG <- function(hg.res, fdr.filter = 1, do.plot = T, show.n = 5, genesets = NULL,
                        str.wrap.width = 25, col.pal = NULL, plt.title = "", pathway.name.size = 7){

  for (i in 1:length(hg.res)){
    hg.res[[i]]$geneset <- names(hg.res)[i]
  }
  res.df <- bind_rows(hg.res)

  if (fdr.filter < 1){
    res.df <- res.df[res.df$padj< fdr.filter,]
  }

  if (is.null(col.pal)){
    col.pal <- categoricalColPal(labels = names(hg.res))
  }

  plt.enrich.list <- list()
  if (do.plot){

    for (i in 1:length(hg.res)){

      current.mod <- names(hg.res)[i]
      res.df.cur <- (res.df %>% dplyr::filter(geneset %in% current.mod) %>% dplyr::arrange(padj))[1:show.n, ]

      if (!is.null(genesets)){
        current.gene.set <- genesets[[current.mod]]
      }


      res.df.cur2 <- bind_rows(
        data.frame(
          xx = res.df.cur$overlap / res.df.cur$size,
          y = gsub("_", " ", res.df.cur$pathway),
          z = res.df.cur$padj,
          ov = res.df.cur$overlap / res.df.cur$size,
          set = "Overlap",
          threshold = 0
        ),
        data.frame(
          xx = log10(res.df.cur$padj),
          y = gsub("_", " ", res.df.cur$pathway),
          z = res.df.cur$padj,
          set = "log10(FDR)",
          ov = res.df.cur$overlap / res.df.cur$size,
          threshold = log10(0.05)
        )
      )

      res.df.cur2$set <- factor(res.df.cur2$set, levels = c("log10(FDR)","Overlap"))

      res.df.cur2 <- res.df.cur2 %>% dplyr::arrange(z)
      res.df.cur2$y <- factor(res.df.cur2$y, levels = unique(res.df.cur2$y))

      plt.enrich <- res.df.cur2 %>%
        ggplot(aes(y = reorder(stringr::str_wrap(y,str.wrap.width),-z), x = xx, fill = -log10(z))) +
        geom_bar(stat = "identity", alpha = 0.8, color = "black") +
        theme_miko(legend = T, center.title = T) +
        scale_fill_gradient(high = col.pal[current.mod], low = "white") +
        ylab("") +
        geom_vline(aes(xintercept = threshold), linetype = "dashed") +
        xlab("Enrichment") +
        labs(title = plt.title, subtitle = current.mod, fill = "-log(FDR)") +
        facet_wrap(~ set, scales = "free_x") +
        scale_x_continuous(expand = c(0, 0),  labels = function(x) signif(abs(x), 3)) +
        theme(panel.spacing.x = unit(0, "mm"),
              axis.text.y = element_text(size = pathway.name.size ))

      plt.enrich.list[[current.mod]] <- plt.enrich

    }

    return(list(
      results = res.df,
      plots = plt.enrich.list
    ))

  }



}
