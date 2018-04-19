#'Convert 2-tailed into 1-tailed p-values from fry
#'
#'Convert 2-tailed into 1-tailed p-values from \code{limma} function \code{fry}.
#'
#'@param tab table with statistics from \code{fry}.
#'@param pv.col name or index of p-value column.
#'@param dir.col name or index of column giving direction gene set has changed.
#'@param direction direction of gene set change. Can be \code{"Up"} or \code{"Down"}.
#'@return p-values.

fry_two2one_tailed <- function(tab, pv.col='PValue', dir.col='Direction', direction='Up'){
  stopifnot(c(tab[,dir.col], direction) %in% c("Up", "Down"), tab[,pv.col]<=1, tab[,pv.col]>=0)
  
  new_pv <- stats::setNames(tab[, pv.col]/2, nm=rownames(tab))
  if (any(tab[,dir.col]!=direction)){
    opp.ind <- which(tab[,dir.col]!=direction)
    new_pv[opp.ind] <- 1 - new_pv[opp.ind]
  }
  return(new_pv)
}