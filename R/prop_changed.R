#' Calculate proportion of changed features per pathway
#' 
#' Calculate proportion of features up & down at p=0.05 per pathway.
#' 
#' @param feat.tab Table of feature (e.g. gene) statistics.
#' @param feat.lst Named list, with each element containing names or indices of rows of \code{feat.tab}.
#' @param stat.col Name or index of column with signed statistics.
#' @param p.col Name or index of p-value column.

prop_changed <- function(feat.tab, feat.lst, stat.col=1, p.col=2){
  stopifnot(feat.tab[,p.col]>=0, feat.tab[,p.col]<=1, ncol(feat.tab) > 1, nrow(feat.tab) > 1, 
            any(feat.tab[,stat.col] < 0))
  
  prop.tab <- t(sapply(feat.lst, FUN=function(x){
    sig <- feat.tab[x, p.col] <= 0.05
    prop.up <- mean(sig & feat.tab[x, stat.col] > 0)
    prop.down <- mean(sig & feat.tab[x, stat.col] < 0)
    return(c(PropUpP05=prop.up, PropDownP05=prop.down))
  }))
  return(prop.tab)
}