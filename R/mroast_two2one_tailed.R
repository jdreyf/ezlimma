#'Convert 2-tailed into 1-tailed p-values from mroast
#'
#'Convert 2-tailed into 1-tailed p-values from \code{limma} function \code{mroast}.
#'
#'@param tab table with statistics from \code{mroast}.
#'@param pv.col name or index of p-value column.
#'@param dir.col name or index of column giving direction gene set has changed.
#'@param direction direction of gene set change. Can be \code{"Up"} or \code{"Down"}.
#'@return Vector of p-values.

mroast_two2one_tailed <- function(tab, pv.col='PValue', dir.col='Direction', direction='Up',
                                  nrot = 9999){
  stopifnot(c(tab[,dir.col], direction) %in% c("Up", "Down"), tab[,pv.col]<=1, 
            tab[,pv.col]>=0, min(tab[,pv.col]) >= 1/(nrot+1) )
  
  #pv=(b+1)/(nrot+1) where b = number of rotations giving a more extreme stat
  #if dir.col==direction, then by symmetry half of the previously extreme rotations 
  #are still extreme, so pv'=(b/2 + 1)/(nrot+1)
  #if dir.col!=direction, then nrot-b/2 of the previously extreme rotations are extreme,
  #which is the equivalent of 1-p/2, so pv' = (nrot - b/2 + 1) / (nrot+1)
  #this is also consistent with 0.5<=pv'<=1
  
  b_vec <- tab[,pv.col]*(nrot+1)-1
  #initialize as if dir.col==direction
  new_pv <- setNames((b_vec/2+1) / (nrot+1), nm=rownames(tab))
  if (any(tab[,dir.col]!=direction)){
    opp.ind <- which(tab[,dir.col]!=direction)
    new_pv[opp.ind] <- (nrot - b_vec[opp.ind]/2 + 1)/(nrot+1)
  }
  return(new_pv)
}