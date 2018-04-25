#' Create object with indices of genes per gene set
#'
#' Create object with indices of genes per gene set.
#' 
#' @param G a gene set list as returned from \code{\link{read_gmt}}.
#' @param object A matrix-like data object containing log-ratios or
#'  log-expression values for a series of arrays, with rows corresponding to
#'  genes and columns to samples.
#' @param min.ngenes minimum number of genes needed in a gene set for testing.
#' @param max.ngenes maximum number of genes needed in a gene set for testing.
#' @return Index of acceptable gene sets in \code{G}.
#' @details This function is not meant to be called directly by the user.

g_index <- function(G, object, min.ngenes, max.ngenes){
  idx <- lapply(G, function(g) rownames(object)[rownames(object) %in% g$genes])
  names(idx) <- sapply(G, function(g) g$name)
  #remove gene sets of the wrong size
  idx <- idx[sapply(idx, function(x) length(x) >= min.ngenes & length(x) <= max.ngenes)]
  if (length(idx)==0) stop("No gene sets are of the right size.")
  return(idx)
}