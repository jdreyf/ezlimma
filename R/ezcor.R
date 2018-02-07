#'Test correlation of each row of an object to a phenotype vector
#'
#'Test correlation of each row of an object to a phenotype vector using one of 
#'several correlation methods.
#'
#'@param object A matrix-like data object containing log-ratios or 
#'  log-expression values, with rows corresponding to 
#'  genes and columns to samples.
#'@param phenotype Vector of phenotypes of the samples. Should be same length as
#'  \code{ncol(object)}. If the vector is named, names should match 
#'  \code{colnames(object)}.
#'@param method a character string indicating which correlation coefficient is
#'  to be used for the test. One of \code{"pearson"}, \code{"kendall"}, or
#'  \code{"spearman"}, can be abbreviated, see \code{\link[stats]{cor.test}}.
#'@param reorder.rows logical, should rows be reordered by p-value or be left in 
#'  the same order as \code{object}?
#'@param prefix character string to add to beginning of column names.
#'@param adjust.method method used to adjust the p-values for multiple testing.
#'@param alternative indicates the alternative hypothesis and must be one of 
#'  \code{"two.sided"}, \code{"greater"}, or \code{"less"}. You can specify just
#'  the initial letter. \code{"greater"} corresponds to positive association, 
#'  \code{"less"} to negative association. See \code{\link[stats]{cor.test}}.
#'@return Dataframe.
#'@export

ezcor <- function(object, phenotype, method="pearson", reorder.rows=TRUE, 
                  prefix=NULL, adjust.method='BH', alternative='two.sided'){
  stopifnot(length(phenotype)==ncol(object), names(phenotype)==colnames(object))
  
  #res.tmp$estimate has cor, res.tmp$statistic has t-stat
  res.mat <- t(apply(X=object, MARGIN=1, FUN=function(v){ 
    res.tmp <- cor.test(v, phenotype, method=method, alternative=alternative)
    return(c(res.tmp$estimate, p=res.tmp$p.value))
  }))
  #adjust p-vals and name column by adjustment
  q <- p.adjust(p=res.mat[,'p'], method=adjust.method)
  res.mat <- cbind(res.mat, adj.P.Val=q)
  if (adjust.method %in% c('BH', 'fdr')){
    colnames(res.mat) <- sub('adj.P.Val', 'FDR', colnames(res.mat))
  } else {
    colnames(res.mat) <- sub('adj.P.Val', adjust.method, colnames(res.mat))
  }
  if (reorder.rows){ res.mat <- res.mat[order(res.mat[,'p']),] }
  if (!is.null(prefix)){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep='.') }
  return(res.mat)
}
