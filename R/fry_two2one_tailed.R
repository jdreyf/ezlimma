#' Convert 2-tailed into 1-tailed p-values from fry
#' 
#' Convert 2-tailed into 1-tailed p-values from \code{limma} function \code{fry}.
#' 
#' @param tab Table with statistics from \code{fry}.
#' @param p.col Name or index of p-value column.
#' @param dir.col Name or index of column giving direction gene set has changed.
#' @param direction Direction of gene set change. Can be \code{"Up"} or \code{"Down"}.
#' @return p-values.

fry_two2one_tailed <- function(tab, p.col="PValue", dir.col="Direction", direction="Up"){
  stopifnot(c(tab[,dir.col], direction) %in% c("Up", "Down"), tab[,p.col]<=1, tab[,p.col]>=0)
  
  new_pv <- stats::setNames(tab[, p.col]/2, nm=rownames(tab))
  if (any(tab[,dir.col]!=direction)){
    opp.ind <- which(tab[,dir.col]!=direction)
    new_pv[opp.ind] <- 1 - new_pv[opp.ind]
  }
  return(new_pv)
}