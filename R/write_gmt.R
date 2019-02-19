#' Write gene set file with GMT extension from a list
#' 
#' Write Gene Matrix Transposed (gmt) file from a list. The first column of
#' the GMT are gene set names, the second are descriptions, then are genes.
#' 
#' @param gmt.lst List of gene sets. Each element has a gene set \code{name}, \code{description}, and \code{genes}.
#' @param file.gmt Name of the GMT file to write to. \code{.gmt} is appended.
#' @return Invisibly, \code{TRUE}.
#' @details This function was adapted from \code{gmtlist2file} in the \pkg{cogena} package.

write_gmt <- function(gmt.lst, file.gmt){
  stopifnot(length(gmt.lst) > 0, is.list(gmt.lst))
  
  for (pwy.ind in seq_along(gmt.lst)){
    cat(gmt.lst[[pwy.ind]]$name, gmt.lst[[pwy.ind]]$description, gmt.lst[[pwy.ind]]$genes, file=file.gmt, append=TRUE, 
        sep = "\t")
    cat("\n", append=TRUE, file=file.gmt)
  }
  return(invisible(TRUE))
}
