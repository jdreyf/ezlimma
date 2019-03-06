#' Create object with indices of genes per gene set
#'
#' Create object with indices of genes per gene set.
#' 
#' @inheritParams limma_contrasts
#' @inheritParams roast_contrasts
#' @return Named list of indices of acceptable gene sets in \code{G}.

g_index <- function(G, object, min.nfeats, max.nfeats){
  idx <- lapply(G, function(g) rownames(object)[rownames(object) %in% g$genes])
  names(idx) <- sapply(G, function(g) g$name)
  # remove gene sets of the wrong size
  idx <- idx[sapply(idx, function(x) length(x) >= min.nfeats & length(x) <= max.nfeats)]
  if (length(idx)==0) stop("No gene sets are of the right size.")
  return(idx)
}