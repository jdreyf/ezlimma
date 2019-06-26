#' Make XL data frame for \code{writexl}
#' 
#' Make pathways table into XL data frame for \code{writexl} by adding hyperlink column with cleaned filenames for URLs.
#' 
#' @inheritParams roast_contrasts
#' @inheritParams write_linked_xl
#' @return A data frame with xl links
#' \describe{
#'   \item{\code{xl}}{}
#'   \item{\code{pwy.csv.nms}}{vector of pathway names to write to csv files}
#' }

xl_pwys <- function(pwy.tab){
  stopifnot(nrow(pwy.tab) > 0, ncol(pwy.tab) > 0)

  pwys <- rownames(pwy.tab)
  # don't allow invalid names in pwys, which are written as filenames
  pwys.clean <- clean_filenames(pwys)
  
  urls <- paste0("pathways/", pwys.clean, ".csv")
  xl_links <- writexl::xl_hyperlink(url=urls, name = pwys)
  xl <- data.frame(xl_links, pwy.tab, stringsAsFactors = FALSE)
  rownames(xl) <- pwys.clean
  # 1st col is rownames
  colnames(xl)[1] <- ""
  return(xl)
}