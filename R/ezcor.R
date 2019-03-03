#' Test correlation of each row of an object to a phenotype vector
#' 
#' Test correlation of each row of an object to a phenotype vector using one of 
#' several correlation methods.
#' 
#' @param method Character string indicating which correlation coefficient to be used for the test. One of 
#' \code{"pearson"}, \code{"kendall"}, or \code{"spearman"}, can be abbreviated, see \code{\link[stats]{cor.test}}.
#' @param alternative Alternative hypothesis, and must be one of\code{"two.sided"}, \code{"greater"}, or \code{"less"}. 
#' You can specify just the initial letter. \code{"greater"} corresponds to positive association, \code{"less"} to 
#' negative association. See \code{\link[stats]{cor.test}}.
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @return Data frame.
#' @seealso \code{\link[ezlimma]{limma_cor}}
#' @export

ezcor <- function(object, phenotype, method="pearson", reorder.rows=TRUE, 
                  prefix=NULL, adjust.method="BH", alternative="two.sided", check.names=TRUE){
  stopifnot(length(phenotype)==ncol(object))
  if (check.names){
    stopifnot(names(phenotype)==colnames(object))
  }
  
  #res.tmp$estimate has cor, res.tmp$statistic has t-stat
  res.mat <- t(apply(X=object, MARGIN=1, FUN=function(v){ 
    res.tmp <- stats::cor.test(v, phenotype, method=method, alternative=alternative)
    return(c(res.tmp$estimate, p=res.tmp$p.value))
  }))
  #adjust p-vals and name column by adjustment
  q <- stats::p.adjust(p=res.mat[,"p"], method=adjust.method)
  res.mat <- cbind(res.mat, adj.P.Val=q)
  if (adjust.method %in% c("BH", "fdr")){
    colnames(res.mat) <- sub("adj.P.Val", "FDR", colnames(res.mat))
  } else {
    colnames(res.mat) <- sub("adj.P.Val", adjust.method, colnames(res.mat))
  }
  if (reorder.rows){ res.mat <- res.mat[order(res.mat[,"p"]),] }
  if (!is.null(prefix)){ colnames(res.mat) <- paste(prefix, colnames(res.mat), sep=".") }
  return(res.mat)
}
