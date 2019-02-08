#' Make XL data frame for \code{writexl}
#' 
#' Make XL data frame for \code{writexl} using top pathways table.
#' 
#' @inheritParams roast_contrasts
#' @inheritParams write_top_xl
#' @return A data frame with xl links
#' \describe{
#'   \item{\code{xl}}{}
#'   \item{\code{pwy.csv.nms}}{vector of pathway names to write to csv files}
#' }

top_xl <- function(pwy.tab){
  stopifnot(nrow(pwy.tab) > 0, ncol(pwy.tab) > 0)

  pwys <- rownames(pwy.tab)
  #don't allow invalid names in pwys, which are written as filenames
  pwys <- clean_filenames(pwys)
  pwys <- substr(pwys, 1, 150)
  
  urls <- paste0('pathways/', pwys, '.csv')
  xl_links <- writexl::xl_hyperlink(url=urls, name = pwys)
  xl <- data.frame(xl_links, pwy.tab, stringsAsFactors = FALSE)
  rownames(xl) <- pwys
  #1st col is rownames
  colnames(xl)[1] <- ""
  
  return(xl)
}