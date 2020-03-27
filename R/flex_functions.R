#' Generate multi-tab datatable list for flexdashboard
#'
#' Prepares list of data.frames for multi-tab table presentation in flexdashboards.
#'
#' @param df.list Named list of data.frames.
#' @param df.list.name Character specifying name of df.list.
#' @name flex.multiTabTables
#' @return flexdashboard compatable list of tables.
#' @examples
#'
#' Tab Header 1
#' =====================================
#'
#' Row {.tabset}
#' -------------------------------------
#'
#' '''{r multi tab table}
#' out <- flex.multiTabTables(df.list, "df.list")
#' '''
#'
#' `r paste(knitr::knit(text = paste(out, collapse = '\n')))`
#'
flex.multiTabTables <- function(df.list, df.list.name) {

  out <- lapply(seq_along(df.list), function(i) {

    # get get name of list
    # list.name <- deparse(substitute(df.list))

    s1 <- paste(df.list.name, "[[", i, "]]", sep = "")

    s4 <- paste("datatable(", s1, ",
    filter = 'top',
    extensions = 'Buttons',
    options = list(pageLength = 50,
    dom = 'Bfrtip',
    buttons = c('copy', 'csv', 'pdf')))", sep = "")

    # get current table name
    table.name <- paste(names(df.list)[i], "")

    a1 <- knitr::knit_expand(text = sprintf("### %s\n", table.name))  # tab header
    a2 <- knitr::knit_expand(text = sprintf("\n```{r %s}", paste(df.list.name,".",table.name,".", i, sep = "")))  # start r chunk
    a3 <- knitr::knit_expand(text = sprintf("\n %s",s4))
    a4 <- knitr::knit_expand(text = "\n```\n") # end r chunk

    paste(a1, a2, a3, a4, collapse = '\n') # collapse together all lines with newline separator

  })

  return(out)

}




#' Generate multi-tab analysis log list for flexdashboard
#'
#' Prepares list of prior analysis logs for multi-tab table presentation in flexdashboards.
#'
#' @param module.logs Character vector specifying name of prior analysis logs available in global enviroment. Must have "df.log_Module_" stem.
#' @name flex.multiTabLogs
#' @return flexdashboard compatable list of analysis logs
#' @examples
#'
#' '''{r multi tab log}
#' out <- flex.multiTabLogs(module.logs)
#' '''
#'
#' `r paste(knitr::knit(text = paste(out, collapse = '\n')))`
#'
flex.multiTabLogs <- function(module.logs) {

  out <- lapply(seq_along(module.logs), function(i) {

    module.n <- as.numeric(gsub("[^\\d]+", "", module.logs[i], perl=TRUE))

    a1 <- knitr::knit_expand(text = sprintf("\nLog (Module %s)", paste(module.n)))
    a2 <- knitr::knit_expand(text = "\n=====================================")
    a3 <- knitr::knit_expand(text = sprintf("\n```{r %s}", paste("mod_", i, sep = ""))) # start r chunk
    a4<- knitr::knit_expand(text = sprintf("\nknitr::kable(%s)",module.logs[i]))
    a5 <- knitr::knit_expand(text = "\n```\n") # end r chunk

    paste(a1, a2, a3, a4, a5, collapse = '\n') # collapse together all lines with newline separator
  })

  return(out)

}


#' Outputs datatable with print buttom options
#'
#' Outputs datatable with print buttom options
#'
#' @param df data.frame
#' @name flex.asDT
#' @return data.table
#'
flex.asDT <- function(df) {

  if (class(df) == "data.frame"){
    dt <-  datatable(df,
                     filter = 'top',
                     extensions = 'Buttons',
                     options = list(pageLength = 50,
                                    dom = 'Bfrtip',
                                    buttons = c('copy', 'csv', 'pdf')))
  } else {
    dt <- NULL
  }


  return(dt)

}


#' Generate multi-tab ggplot handle list for flexdashboard
#'
#' Prepares list of ggplot handles for multi-tab plot presentation in flexdashboards.
#'
#' @param plt.list list of ggplot handles
#' @param plt.list.name list name
#' @param fig.width Numeric. Figure width. Default is 5.
#' @param fig.height Numeric. Figure width. Default is 5.
#' @name flex.multiTabPlot
#' @return flexdashboard compatable list of plots
#'
flex.multiTabPlot <- function(plt.list, plt.list.name, fig.width = 5, fig.height = 5) {

  out <- lapply(seq_along(plt.list), function(i) {
    a1 <- knitr::knit_expand(text = sprintf("### %s\n", names(plt.list)[i])) # tab header
    a2 <- knitr::knit_expand(text = sprintf("\n```{r %s, fig.width=%d, fig.height=%d}", paste(i, plt.list.name), fig.width, fig.height)) # start r chunk
    a3 <- knitr::knit_expand(text = sprintf("\nprint(plt.list[[%d]])", i))
    a4 <- knitr::knit_expand(text = "\n```\n") # end r chunk
    paste(a1, a2, a3, a4, collapse = '\n') # collapse together all lines with newline separator

  })
  return(out)
}



#' Generate multi-tab list of plotly figures for flexdashboard
#'
#' Prepares list of ggplot handles as plotly figures for multi-tab presentation in flexdashboards.
#'
#' @param plt.list list of ggplot handles
#' @param plt.list.name list name
#' @param fig.width Numeric. Figure width. Default is 5.
#' @param fig.height Numeric. Figure width. Default is 5.
#' @name flex.multiTabPlotly
#' @return data.table
#'
flex.multiTabPlotly <- function(plt.list, plt.list.name, fig.width = 5, fig.height = 5) {
  out <- lapply(seq_along(plt.list), function(i) {

    s1 <- paste("plotly::config(ggplotly(plt.list[[", i, "]]), toImageButtonOptions = list(format = 'svg', filename = 'myplot', width = 600, height = 700))", sep = "")

    a1 <- knitr::knit_expand(text = sprintf("### %s\n", names(plt.list)[i])) # tab header
    a2 <- knitr::knit_expand(text = sprintf("\n```{r %s,  fig.width=%d, fig.height=%d}", paste(i, plt.list.name), fig.width, fig.height)) # start r chunk
    a3 <- knitr::knit_expand(text = sprintf("\n %s", s1))
    a4 <- knitr::knit_expand(text = "\n```\n") # end r chunk
    paste(a1, a2, a3, a4, collapse = '\n') # collapse together all lines with newline separator

  })
  return(out)
}
