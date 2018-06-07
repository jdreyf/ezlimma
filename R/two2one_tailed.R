#' Convert 2-tailed into 1-tailed p-values
#' 
#' Convert 2-tailed into 1-tailed p-values from matrix. Alternative "greater" corresponds to positive association, 
#' "less" to negative association.
#' 
#' @param tab table with statistics from \code{fry}.
#' @param stat.col name or index of column with signed statistics.
#' @param p.col name or index of p-value column.
#' @param alternative direction of change. Can be a vector. \code{"greater"} or \code{"less"}.
#' @return p-values.
#' @details This function is not meant to be called directly by the user.

two2one_tailed <- function(tab, stat.col=1, p.col=2, alternative=c("greater", "less")){
  alternative <- match.arg(alternative)
  stopifnot(alternative %in% c("less", "greater"), tab[,p.col]>=0, tab[,p.col]<=1)
  if (alternative=='less'){ tab[,stat.col] <- -1*tab[,stat.col] }
  new_pv <- stats::setNames(tab[,p.col]/2, nm=rownames(tab))
  if (any(tab[,stat.col] < 0)){
    new_pv[tab[,stat.col] < 0] <- 1 - new_pv[tab[,stat.col] < 0]
  }
  return(new_pv)
}