#' Make XL data frame for \code{writexl}
#' 
#' Make XL data frame for \code{writexl} using top pathways table.
#' 
#' @inheritParams roast_contrasts
#' @inheritParams write_top_xl
#' @return A list with 2 elements
#' \describe{
#'   \item{\code{xl}}{data frame with xl links}
#'   \item{\code{pwy.csv.nms}}{vector of pathway names to write to csv files}
#' }

top_xl <- function(pwy.tab, n.toptabs=Inf){
  stopifnot(nrow(pwy.tab) > 0, ncol(pwy.tab) > 0, is.numeric(n.toptabs), n.toptabs > 0)
  if (n.toptabs > nrow(pwy.tab)) n.toptabs <- nrow(pwy.tab)
  pwys.csv <- rownames(pwy.tab)[1:n.toptabs]
  #don't allow invalid names in pwys, which are written as filenames
  pwys.csv <- clean_filenames(pwys.csv)
  pwy.csv.nms <- substr(pwys.csv, 1, 150)
  
  urls <- paste0('pathways/', pwy.csv.nms, '.csv')
  xl_links <- writexl::xl_hyperlink(url=urls, name = pwy.csv.nms)
  if (n.toptabs < nrow(pwy.tab)){
    pwy.noncsv.nms <- rownames(pwy.tab)[-(1:n.toptabs)]
    xl_links <- c(xl_links, pwy.noncsv.nms)
  }
  xl <- data.frame(xl_links, pwy.tab)
  #1st col is rownames
  colnames(xl)[1] <- ""
  
  return(list(xl=xl, pwy.csv.nms=pwy.csv.nms))
}