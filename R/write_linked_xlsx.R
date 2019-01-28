#' Write Excel XLSX file with links to CSVs
#' 
#' Create directory and write CSV files with stats from \code{stats.tab} per pathway, and write an Excel XLSX file 
#' using package \code{writexl} that links to the CSVs.
#' 
#' @param name A name for the folder and Excel file. If \code{NA}, a directory is not created and no file is written.
#' @param res A table of pathway results.
#' @param index Named list of Indices of pathways in G to write.
#' @param stats.tab A table of feature (e.g. gene) statistics that the Excel file can link to.
#' @param n.toptabs Number of pathway toptables to write to CSV and link to from Excel.
#' @return Invisibly, the data frame that's written to Excel.

write_linked_xlsx <- function(name, res, index, stats.tab, n.toptabs){
  stopifnot(!is.null(name), nrow(res) > 0, !is.null(names(index)), nrow(stats.tab) > 0, is.numeric(n.toptabs),
            n.toptabs > 0)
  if (!requireNamespace("writexl", quietly = TRUE)){
    stop("Install 'writexl' package.", call. = FALSE)
  }
  
  if (n.toptabs > nrow(res)) n.toptabs <- nrow(res)
  pwys <- rownames(res)[1:n.toptabs]
  #don't allow invalid names in pwys, which are written as filenames
  pwys <- clean_filenames(pwys)
  pwy.nms <- substr(pwys, 1, 150)
  
  urls <- paste0('pathways/', pwy.nms, '.csv')
  xl_links <- writexl::xl_hyperlink(url=urls, name = )
  xl <- data.frame(xl_links, res)
  #1st col is rownames
  colnames(xl)[1] <- ""
  
  if (!is.na(name)){
    dir.create(name)
    dir.create(paste0(name, '/pathways'))
    names(index) <- clean_filenames(names(index))
    for(pwy in pwy.nms){
      stat <- stats.tab[index[[pwy]], ]
      stat <- stat[order(combine_pvalues(stat)), ]
      utils::write.csv(stat, paste0(name, '/pathways/', pwy, '.csv'))
    }
    writexl::write_xlsx(x=xl, path = paste0(name, "/", name, ".xlsx"))
  }#end if
  return(invisible(xl))
}