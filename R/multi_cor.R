#' Test correlation of each row of an object to each column of pheno.tab
#' 
#' Test correlation of each row of an object to each column of \code{pheno.tab} using 
#' one of Pearson's, Kendall's, or Spearman's correlation methods, or limma 
#' regression in \code{\link{limma_cor}}. See examples in vignette.
#' 
#' @param pheno.tab Matrix-like data object with columns as sample phenotypes, with \code{nrow(pheno.tab)==ncol(object)}.
#' @param method Character string indicating which association is to be used 
#'   for the test. One of \code{"pearson"}, \code{"spearman"}, \code{"kendall"}, 
#'   from \code{\link[stats]{cor.test}} or \code{"limma"} for \code{\link{limma_cor}}.
#' @param correlation Numeric vector of inter-duplicate or inter-technical replicate correlations. Must be given if 
#' \code{!is.null(block)}. Its length should be the same as the number of columns of \code{pheno.tab}.
#' @param covariates If \code{method="limma"}, numeric vector or matrix of covariates to include in 
#' \code{\link{limma_cor}} \code{design} matrix.
#' @param check.names Logical; should \code{rownames(pheno.tab)=colnames(object)} be checked?
#' @param limma.cols If \code{method="limma"}, \code{cols} from \code{\link{limma_cor}} to include.
#' @inheritParams limma_contrasts
#' @inheritParams limma_cor
#' @return Data frame with several statistical columns corresponding to each phenotype and one row per feature.
#' @details Each column of \code{pheno.tab} is tested independently. Arguments \code{covariates}, \code{block}, and
#' \code{correlation} only apply if \code{method="limma"}. When each individual \code{pheno.tab} column is tested, 
#' if some samples have \code{NA}s for that column, those samples are omitted for that column only. 
#' @export

multi_cor <- function(object, pheno.tab, method=c("pearson", "spearman", "kendall", "limma"), reorder.rows=TRUE, 
                      prefix=NULL, block=NULL, correlation=NULL, adjust.method="BH", covariates=NULL, check.names=TRUE,
                     limma.cols=c("AveExpr", "P.Value", "adj.P.Val", "logFC")){
  method <- match.arg(method)
  if (is.null(dim(pheno.tab))) stop("pheno.tab needs to have rows and columns.")
  stopifnot(ncol(object)==nrow(pheno.tab), is.null(covariates) || limma::isNumeric(covariates), colMeans(is.na(pheno.tab)) < 1,
            is.null(block) || length(correlation) == ncol(pheno.tab))
  if (check.names){
    stopifnot(rownames(pheno.tab)==colnames(object))
  }
  
  cor.mat <- NULL
  for (ind in 1:ncol(pheno.tab)){
    prefix.tmp <- ifelse(!is.null(prefix), paste(prefix, colnames(pheno.tab)[ind], sep="."), colnames(pheno.tab)[ind])
    if (method=="limma"){
      # na.omit() missing phenotypes in pheno.tab
      ph.idx <- which(!is.na(pheno.tab[,ind]))
      
      # block is on samples, so should be modified in case of NAs, & NULL[ph.idx] is still NULL
      if (is.null(covariates)){
        cor.tmp <- data.matrix(limma_cor(object[, ph.idx], pheno.tab[ph.idx, ind], reorder.rows=FALSE, prefix=prefix.tmp, block = block[ph.idx],
                                         correlation = correlation[ind], cols=limma.cols))
      } else {
        # model.matrix.lm, but not model.matrix, respects na.action
        # https://stackoverflow.com/questions/5616210/model-matrix-with-na-action-null
        des.tmp <- stats::model.matrix.lm(~1+pheno.tab[, ind]+covariates, na.action = stats::na.omit)
        cor.tmp <- data.matrix(limma_cor(object[, ph.idx], design=des.tmp, reorder.rows=FALSE, prefix=prefix.tmp, block = block[ph.idx],
                                         correlation = correlation[ind], cols=limma.cols))
      }
    } else {
      if (!is.null(covariates) || !is.null(block) || !is.null(correlation)){
        warning("method is not limma, so limma-specific arguments ignored.")
      }
      cor.tmp <- ezcor(object, pheno.tab[,ind], method=method, reorder.rows=FALSE, prefix=prefix.tmp, adjust.method=adjust.method)
    }
    if (!is.null(cor.mat)){ stopifnot(rownames(cor.mat)==rownames(cor.tmp)) }
    cor.mat <- cbind(cor.mat, cor.tmp)
  }
  if (reorder.rows){ cor.mat <- cor.mat[order(combine_pvalues(cor.mat)),] }
  # return data frame, for consistency with limma_*
  return(data.frame(cor.mat))
}