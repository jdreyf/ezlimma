#' Write Excel XLSX file with links to CSVs
#' 
#' Create directory and write CSV files with stats from \code{feat.tab} per pathway, and write an Excel XLSX file 
#' using package \code{writexl} that links to the CSVs.
#' 
#' @param pwy.tab A table of pathway results.
#' @inheritParams prop_changed
#' @inheritParams roast_contrasts
#' @return Invisibly, the data frame that's written to Excel.

write_linked_xlsx <- function(pwy.tab, feat.lst, feat.tab, name=NA, n.toptabs=Inf){
  stopifnot(!is.null(name), nrow(pwy.tab) > 0, !is.null(names(feat.lst)), nrow(feat.tab) > 0, is.numeric(n.toptabs),
            n.toptabs > 0)
  if (!requireNamespace("writexl", quietly = TRUE)){
    stop("Install 'writexl' package.", call. = FALSE)
  }
  
  if (n.toptabs > nrow(pwy.tab)) n.toptabs <- nrow(pwy.tab)
  pwys <- rownames(pwy.tab)[1:n.toptabs]
  #don't allow invalid names in pwys, which are written as filenames
  pwys <- clean_filenames(pwys)
  pwy.nms <- substr(pwys, 1, 150)
  
  urls <- paste0('pathways/', pwy.nms, '.csv')
  xl_links <- writexl::xl_hyperlink(url=urls, name = )
  xl <- data.frame(xl_links, pwy.tab)
  #1st col is rownames
  colnames(xl)[1] <- ""
  
  if (!is.na(name)){
    dir.create(name)
    dir.create(paste0(name, '/pathways'))
    names(feat.lst) <- clean_filenames(names(feat.lst))
    for(pwy in pwy.nms){
      stat <- feat.tab[feat.lst[[pwy]], ]
      stat <- stat[order(combine_pvalues(stat)), ]
      utils::write.csv(stat, paste0(name, '/pathways/', pwy, '.csv'))
    }
    writexl::write_xlsx(x=xl, path = paste0(name, "/", name, ".xlsx"))
  }#end if
  return(invisible(xl))
}