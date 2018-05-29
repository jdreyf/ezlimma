#' Test partial regression of each row of an object to a phenotype vector
#' 
#' Test partial regression of each row of an object to a phenotype vector given covariates in a linear 
#' regression of \code{phenotype ~ 1 + x + covar} for \code{x} being each row of \code{object}.
#' 
#' @param object A matrix-like data object containing log-ratios or 
#'   log-expression values, with rows corresponding to 
#'   genes and columns to samples.
#' @param phenotype Vector of phenotypes of the samples. Should be same length as
#'   \code{ncol(object)}. If the vector is named, names should match 
#'   \code{colnames(object)}.
#' @param covar A matrix or dataframe with one or more covariates as columns.
#' @param reorder.rows logical, should rows be reordered by p-value or be left in 
#'   the same order as \code{object}?
#' @param prefix character string to add to beginning of column names.
#' @param adjust.method method used to adjust the p-values for multiple testing.
#' @return Dataframe.
#' @import stats

ezpreg <- function(object, phenotype, covar, reorder.rows=TRUE, prefix=NULL, adjust.method='BH'){
  stopifnot(length(phenotype)==ncol(object), names(phenotype)==colnames(object), nrow(covar)==ncol(object))
  
  #res.tmp$estimate has cor, res.tmp$statistic has t-stat
  res.mat <- t(apply(X=object, MARGIN=1, FUN=function(v){ 
    form.mat <- data.frame(phenotype, v=v, covar)
    fm <- stats::lm(phenotype ~ ., data=form.mat)
    sfm <- summary(fm)$coefficients
    return(c(slope=sfm["v", "Estimate"], p=sfm["v", "Pr(>|t|)"]))
  }))
  #adjust p-vals and name column by adjustment
  q <- stats::p.adjust(p=res.mat[,'p'], method=adjust.method)
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
