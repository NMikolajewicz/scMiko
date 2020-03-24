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
