#' Write Excel XLSX file with links to CSVs
#' 
#' Create directory and write CSV files with stats from \code{feat.tab} per pathway, and write an Excel XLSX file 
#' using package \code{writexl} that links to the CSVs.
#' 
#' @param pwy.tab A table of pathway results.
#' @inheritParams prop_changed
#' @inheritParams roast_contrasts
#' @return Invisibly, the data frame that's written to Excel.

write_top_xl <- function(pwy.tab, feat.lst, feat.tab, name=NA, n.toptabs=Inf){
  stopifnot(!is.null(name), nrow(pwy.tab) > 0, !is.null(feat.lst), !is.null(names(feat.lst)), nrow(feat.tab) > 0, 
            is.numeric(n.toptabs), n.toptabs > 0)
  if (!requireNamespace("writexl", quietly = TRUE)){
    stop("Install 'writexl' package.", call. = FALSE)
  }
  
  tx <- top_xl(pwy.tab=pwy.tab, n.toptabs=n.toptabs)
  
  if (!is.na(name)){
    dir.create(name)
    dir.create(paste0(name, '/pathways'))
    names(feat.lst) <- clean_filenames(names(feat.lst))
    for(pwy in tx$pwy.csv.nms){
      stat <- feat.tab[feat.lst[[pwy]], ]
      stat <- stat[order(combine_pvalues(stat)), ]
      utils::write.csv(stat, paste0(name, '/pathways/', pwy, '.csv'))
    }
    writexl::write_xlsx(x=tx$xl, path = paste0(name, "/", name, ".xlsx"))
  }#end if
  return(invisible(tx$xl))
}