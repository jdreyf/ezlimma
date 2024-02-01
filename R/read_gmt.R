#' Read gene set file with GMT extension to a list
#' 
#' Read Gene Matrix Transposed (gmt) file and return a list. The first column of
#' the GMT are gene set names, the second are descriptions, then are genes.
#' 
#' @param file.path Path to the GMT file.
#' @param strip.nm Regular expression which is removed from gene set names. Case is ignored.
#' @details There are many functions in R for reading GMT files; this one was partially adapted from \code{gmt2list} 
#' in the \pkg{cogena} package.
#' @return List of gene sets. Each element has a gene set \code{name}, \code{description}, and \code{genes}.
#' @export

read_gmt <- function(file.path, strip.nm="REACTOME_"){
  x00 <- scan(file.path, what="", sep="\n", quiet=TRUE)
  # sometimes have literal \" inside character string name
  x0 <- gsub(pattern='\"', replacement = '', x=x00, fixed = TRUE)
  x1 <- strsplit(x=x0, split="\t")
  x2 <- lapply(x1, FUN=function(v){
    list(name=gsub(pattern = strip.nm, replacement = "", x = v[1], ignore.case=TRUE), description=v[2], genes=v[c(-1,-2)])
  })
  names(x2) <- lapply(x2, FUN=function(vv) vv$name)
  return(x2)
}
