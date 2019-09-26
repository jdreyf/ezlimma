#' Write Excel XLSX file with links to CSVs
#' 
#' Create directory and write CSV files with stats from \code{feat.tab} per pathway, and write an Excel XLSX file 
#' using package \code{writexl} that links to the CSVs.
#' 
#' @param pwy.tab A table of pathway results.
#' @inheritParams prop_changed
#' @inheritParams roast_contrasts
#' @return Invisibly, the data frame that's written to Excel.

write_linked_xl <- function(pwy.tab, feat.lst, feat.tab, name){
  stopifnot(!is.null(name), nrow(pwy.tab) > 0, !is.null(feat.lst), !is.null(names(feat.lst)), nrow(feat.tab) > 0,
            rownames(pwy.tab) %in% names(feat.lst))
  if (!requireNamespace("writexl", quietly = TRUE)){
    stop("Install 'writexl' package.", call. = FALSE)
  }
  
  feat.lst <- feat.lst[rownames(pwy.tab)]
  tx <- xl_pwys(pwy.tab=pwy.tab)
  names(feat.lst) <- rownames(tx)
  
  if (file.exists(name)) unlink(name, recursive = TRUE)
  
  dir.create(name)
  dir.create(paste0(name, '/pathways'))
  
  for(pwy in rownames(tx)){
    stat <- feat.tab[feat.lst[[pwy]],, drop=FALSE]
    if (any(grepl(pattern="p|PValue", x=colnames(stat)))) stat <- stat[order(combine_pvalues(stat)),, drop=FALSE]
    utils::write.csv(stat, paste0(name, '/pathways/', pwy, '.csv'))
  }
  writexl::write_xlsx(x=tx, path = paste0(name, "/", name, ".xlsx"))
  
  return(invisible(tx))
}